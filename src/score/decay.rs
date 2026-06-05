//! Global decay score aggregation.

use rayon::prelude::*;
use tracing::error;

use crate::config::weights::WeightsV1;
use crate::score::AxisScoresVec;
use crate::util::numeric::clamp01;

/// Decay score and robustness margin vectors.
#[derive(Debug, Clone, PartialEq)]
pub struct DecayScoreVec {
    pub decay: Vec<f32>,
    pub robustness_margin: Vec<f32>,
}

/// Compute decay score and robustness margin for v1.
///
/// Global weight sum is validated at config load time.
pub fn compute_decay_v1(axes: &AxisScoresVec, weights: &WeightsV1) -> DecayScoreVec {
    let len = axes.bioenergetics.len();
    assert_eq!(axes.ros.len(), len);
    assert_eq!(axes.dynamics.len(), len);
    assert_eq!(axes.regulation.len(), len);

    let w_bio = weights.global.bioenergetics;
    let w_ros = weights.global.ros;
    let w_dyn = weights.global.dynamics;
    let w_reg = weights.global.regulation;

    let mut decay = vec![0.0; len];
    let mut robustness_margin = vec![0.0; len];

    decay
        .par_iter_mut()
        .zip(robustness_margin.par_iter_mut())
        .enumerate()
        .with_min_len(1024)
        .for_each(|(i, (d, r))| {
            let value = w_bio * axes.bioenergetics[i]
                + w_ros * axes.ros[i]
                + w_dyn * axes.dynamics[i]
                + w_reg * axes.regulation[i];

            if value.is_nan() {
                error!(sample = i, "NaN encountered in decay score");
                *d = 0.0;
                *r = 0.0;
            } else {
                let clamped = clamp01(value);
                *d = clamped;
                *r = clamp01(1.0 - clamped);
            }
        });

    DecayScoreVec {
        decay,
        robustness_margin,
    }
}
