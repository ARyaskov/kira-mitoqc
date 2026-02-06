//! Failure-mode classification for v1.

use crate::config::refs::RefsV1;
use crate::core::types::MitochondrialState;
use crate::score::{AxisScoresVec, DecayScoreVec};

pub mod v1;

/// Classified states for a profile.
#[derive(Debug, Clone, PartialEq)]
pub struct ClassifiedProfile {
    pub states: Vec<MitochondrialState>,
}

/// Classify per-sample states for v1.
pub fn classify_v1(
    axes: &AxisScoresVec,
    decay: &DecayScoreVec,
    refs: &RefsV1,
) -> Vec<MitochondrialState> {
    v1::classify_v1(axes, decay, refs)
}

/// Convenience helper returning a wrapped profile.
pub fn classify_profile_v1(
    axes: &AxisScoresVec,
    decay: &DecayScoreVec,
    refs: &RefsV1,
) -> ClassifiedProfile {
    ClassifiedProfile {
        states: classify_v1(axes, decay, refs),
    }
}
