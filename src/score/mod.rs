//! Score aggregation for axes and decay.

use crate::config::weights::WeightsV1;
use crate::core::types::ProxyScores;

pub mod axis;
pub mod decay;
pub mod v2;

pub use axis::{AxisScoresVec, compute_axes_v1};
pub use decay::{DecayScoreVec, compute_decay_v1};
pub use v2::{compute_axes_v2, compute_decay_v2};

/// Combined score output for a profile.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoredProfile {
    pub axes: AxisScoresVec,
    pub decay: DecayScoreVec,
}

/// Compute axes and decay scores for v1.
pub fn score_profile_v1(proxies: &ProxyScores, weights: &WeightsV1) -> ScoredProfile {
    let axes = compute_axes_v1(proxies, weights);
    let decay = compute_decay_v1(&axes, weights);
    ScoredProfile { axes, decay }
}
