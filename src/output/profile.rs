//! Profile assembly for v1 output.

use serde::{Deserialize, Serialize};

use crate::core::types::{AxisScores, MitochondrialState, ProxyScores};
use crate::explain::{Driver, Explainability};
use crate::score::{AxisScoresVec, DecayScoreVec};

/// Serializable profile output for v1.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MitoProfileV1 {
    pub mitochondrial_state: String,
    pub decay_score: f32,
    pub robustness_margin: f32,
    pub axes: AxisScores,
    pub proxies: ProxyScores,
    pub drivers: Vec<Driver>,
    pub interpretation: Vec<String>,
}

/// Assemble per-sample profiles from computed components.
pub fn assemble_profiles_v1(
    states: &[MitochondrialState],
    decay: &DecayScoreVec,
    axes: &AxisScoresVec,
    proxies: &ProxyScores,
    explain: &[Explainability],
) -> Vec<MitoProfileV1> {
    let len = states.len();
    assert_eq!(decay.decay.len(), len);
    assert_eq!(decay.robustness_margin.len(), len);
    assert_eq!(axes.bioenergetics.len(), len);
    assert_eq!(axes.ros.len(), len);
    assert_eq!(axes.dynamics.len(), len);
    assert_eq!(axes.regulation.len(), len);
    assert_eq!(explain.len(), len);

    let mut profiles = Vec::with_capacity(len);
    for i in 0..len {
        let axes_scalar = AxisScores {
            bioenergetics: axes.bioenergetics[i],
            ros: axes.ros[i],
            dynamics: axes.dynamics[i],
            regulation: axes.regulation[i],
        };
        let proxies_single = slice_proxy_scores(proxies, i);
        let explain_sample = &explain[i];

        profiles.push(MitoProfileV1 {
            mitochondrial_state: states[i].as_str().to_string(),
            decay_score: decay.decay[i],
            robustness_margin: decay.robustness_margin[i],
            axes: axes_scalar,
            proxies: proxies_single,
            drivers: explain_sample.drivers.clone(),
            interpretation: explain_sample.interpretation.clone(),
        });
    }

    profiles
}

pub(crate) fn slice_proxy_scores(proxies: &ProxyScores, sample: usize) -> ProxyScores {
    let mut out = ProxyScores::default();

    for (key, values) in proxies.normalized.iter() {
        out.set(*key, vec![values[sample]]);
    }
    for (key, values) in proxies.raw.iter() {
        out.set_raw(*key, vec![values[sample]]);
    }

    out
}
