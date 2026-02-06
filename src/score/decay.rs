//! Global decay score aggregation.

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
pub fn compute_decay_v1(axes: &AxisScoresVec, weights: &WeightsV1) -> DecayScoreVec {
    debug_assert!(
        (weights.global.bioenergetics
            + weights.global.ros
            + weights.global.dynamics
            + weights.global.regulation
            - 1.0)
            .abs()
            <= 1e-6
    );

    let len = axes.bioenergetics.len();
    assert_eq!(axes.ros.len(), len);
    assert_eq!(axes.dynamics.len(), len);
    assert_eq!(axes.regulation.len(), len);

    let mut decay = vec![0.0; len];
    let mut robustness_margin = vec![0.0; len];

    for i in 0..len {
        let value = weights.global.bioenergetics * axes.bioenergetics[i]
            + weights.global.ros * axes.ros[i]
            + weights.global.dynamics * axes.dynamics[i]
            + weights.global.regulation * axes.regulation[i];

        if value.is_nan() {
            error!(sample = i, "NaN encountered in decay score");
            decay[i] = 0.0;
            robustness_margin[i] = 0.0;
        } else {
            let clamped = clamp01(value);
            decay[i] = clamped;
            robustness_margin[i] = clamp01(1.0 - clamped);
        }
    }

    DecayScoreVec {
        decay,
        robustness_margin,
    }
}
