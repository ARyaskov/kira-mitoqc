//! Explainability for v1 outputs.

use crate::config::weights::WeightsV1;
use crate::core::types::{MitochondrialState, ProxyKey, ProxyScores};
use crate::score::{AxisScoresVec, DecayScoreVec};

pub mod v1;

/// Axis grouping for explainability.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub enum AxisKind {
    Bioenergetics,
    Ros,
    Dynamics,
    Regulation,
}

/// A single driver contribution.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Driver {
    pub key: ProxyKey,
    pub axis: AxisKind,
    pub contribution: f32,
}

/// Explainability output per sample.
#[derive(Debug, Clone, PartialEq)]
pub struct Explainability {
    pub drivers: Vec<Driver>,
    pub interpretation: Vec<String>,
}

/// Compute explainability for v1.
pub fn explain_v1(
    proxies: &ProxyScores,
    axes: &AxisScoresVec,
    decay: &DecayScoreVec,
    states: &[MitochondrialState],
    weights: &WeightsV1,
) -> Vec<Explainability> {
    v1::explain_v1(proxies, axes, decay, states, weights)
}
