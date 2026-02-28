//! v1 rule-based classification.

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

    let mut states = Vec::with_capacity(len);
    for i in 0..len {
        if let Some(redox_metrics) = redox {
            match redox_metrics.redox_regime[i] {
                RedoxRegime::RedoxOverload => {
                    states.push(MitochondrialState::RedoxOverload);
                    continue;
                }
                RedoxRegime::UnbufferedOxidativeStress => {
                    states.push(MitochondrialState::UnbufferedOxidativeStress);
                    continue;
                }
                RedoxRegime::CompensatedOxidativeStress => {
                    states.push(MitochondrialState::CompensatedOxidativeStress);
                    continue;
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
            states.push(MitochondrialState::CompensatedButFragile);
            continue;
        }

        if ros > refs.thresholds.ros_high && bio < refs.thresholds.bioenergetics_low {
            states.push(MitochondrialState::RosDominantDecay);
        } else if bio > refs.thresholds.bioenergetics_high {
            states.push(MitochondrialState::BioenergeticCollapse);
        } else if dyns > refs.thresholds.dynamics_high && reg < refs.thresholds.regulation_low {
            states.push(MitochondrialState::MitophagyLockedDepletion);
        } else if bio > refs.thresholds.structural_bio_min
            && dyns > refs.thresholds.structural_dyn_min
        {
            states.push(MitochondrialState::StructuralFragmentation);
        } else {
            states.push(MitochondrialState::CompensatedButFragile);
        }
    }

    states
}
