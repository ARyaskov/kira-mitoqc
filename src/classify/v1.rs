//! v1 rule-based classification.

use rayon::prelude::*;
use tracing::error;

use crate::config::refs::RefsV1;
use crate::core::types::MitochondrialState;
use crate::redox::{RedoxMetrics, RedoxRegime};
use crate::score::{AxisScoresVec, DecayScoreVec};

/// Classify per-sample states for v1.
pub fn classify_v1(
    axes: &AxisScoresVec,
    decay: &DecayScoreVec,
    refs: &RefsV1,
) -> Vec<MitochondrialState> {
    classify_v1_with_redox(axes, decay, refs, None)
}

/// Redox-extended deterministic classification.
///
/// Priority order:
/// 1) explicit redox overload
/// 2) explicit unbuffered oxidative stress
/// 3) explicit compensated oxidative stress
/// 4) legacy v1 rules (backward-compatible)
pub fn classify_v1_with_redox(
    axes: &AxisScoresVec,
    _decay: &DecayScoreVec,
    refs: &RefsV1,
    redox: Option<&RedoxMetrics>,
) -> Vec<MitochondrialState> {
    let len = axes.bioenergetics.len();
    assert_eq!(axes.ros.len(), len);
    assert_eq!(axes.dynamics.len(), len);
    assert_eq!(axes.regulation.len(), len);

    if let Some(redox_metrics) = redox {
        assert_eq!(redox_metrics.redox_regime.len(), len);
    }

    let t = &refs.thresholds;

    (0..len)
        .into_par_iter()
        .with_min_len(1024)
        .map(|i| {
            if let Some(redox_metrics) = redox {
                match redox_metrics.redox_regime[i] {
                    RedoxRegime::RedoxOverload => return MitochondrialState::RedoxOverload,
                    RedoxRegime::UnbufferedOxidativeStress => {
                        return MitochondrialState::UnbufferedOxidativeStress;
                    }
                    RedoxRegime::CompensatedOxidativeStress => {
                        return MitochondrialState::CompensatedOxidativeStress;
                    }
                    RedoxRegime::Baseline => {}
                }
            }

            let bio = axes.bioenergetics[i];
            let ros = axes.ros[i];
            let dyns = axes.dynamics[i];
            let reg = axes.regulation[i];

            if bio.is_nan() || ros.is_nan() || dyns.is_nan() || reg.is_nan() {
                error!(sample = i, "NaN encountered in classification");
                return MitochondrialState::CompensatedButFragile;
            }

            if ros > t.ros_high && bio < t.bioenergetics_low {
                MitochondrialState::RosDominantDecay
            } else if bio > t.bioenergetics_high {
                MitochondrialState::BioenergeticCollapse
            } else if dyns > t.dynamics_high && reg < t.regulation_low {
                MitochondrialState::MitophagyLockedDepletion
            } else if bio > t.structural_bio_min && dyns > t.structural_dyn_min {
                MitochondrialState::StructuralFragmentation
            } else {
                MitochondrialState::CompensatedButFragile
            }
        })
        .collect()
}
