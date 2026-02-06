//! v1 rule-based classification.

use tracing::error;

use crate::config::refs::RefsV1;
use crate::core::types::MitochondrialState;
use crate::score::{AxisScoresVec, DecayScoreVec};

/// Classify per-sample states for v1.
pub fn classify_v1(
    axes: &AxisScoresVec,
    _decay: &DecayScoreVec,
    refs: &RefsV1,
) -> Vec<MitochondrialState> {
    let len = axes.bioenergetics.len();
    assert_eq!(axes.ros.len(), len);
    assert_eq!(axes.dynamics.len(), len);
    assert_eq!(axes.regulation.len(), len);

    let mut states = Vec::with_capacity(len);
    for i in 0..len {
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
