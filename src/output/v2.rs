//! v2 output assembly.

use serde::{Deserialize, Serialize};

use crate::core::types::AxisScores;
use crate::output::profile::MitoProfileV1;
use crate::proxy::ProxyScoresV2;
use crate::score::{AxisScoresVec, DecayScoreVec};

/// v2 output bundle for a single sample.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct V2Profile {
    pub axes: AxisScores,
    pub decay_score: f32,
    pub robustness_margin: f32,
    pub proxies: ProxyScoresV2,
    pub refs_version: String,
}

/// Combined v1 + v2 output per sample.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MitoProfileBundleV2 {
    pub v1: MitoProfileV1,
    pub v2: V2Profile,
}

/// Assemble v2 bundle output.
pub fn assemble_profiles_v2(
    v1_profiles: &[MitoProfileV1],
    axes: &AxisScoresVec,
    decay: &DecayScoreVec,
    proxies: &ProxyScoresV2,
) -> Vec<MitoProfileBundleV2> {
    let len = v1_profiles.len();
    assert_eq!(axes.bioenergetics.len(), len);
    assert_eq!(decay.decay.len(), len);

    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        let axes_scalar = AxisScores {
            bioenergetics: axes.bioenergetics[i],
            ros: axes.ros[i],
            dynamics: axes.dynamics[i],
            regulation: axes.regulation[i],
        };
        let v2_profile = V2Profile {
            axes: axes_scalar,
            decay_score: decay.decay[i],
            robustness_margin: decay.robustness_margin[i],
            proxies: slice_proxy_scores_v2(proxies, i),
            refs_version: "v2".to_string(),
        };
        out.push(MitoProfileBundleV2 {
            v1: v1_profiles[i].clone(),
            v2: v2_profile,
        });
    }

    out
}

fn slice_proxy_scores_v2(proxies: &ProxyScoresV2, sample: usize) -> ProxyScoresV2 {
    let mut out = ProxyScoresV2::default();
    out.v1 = super::profile::slice_proxy_scores(&proxies.v1, sample);

    for (key, values) in proxies.v2_normalized.iter() {
        out.v2_normalized.insert(*key, vec![values[sample]]);
    }
    for (key, values) in proxies.v2_raw.iter() {
        out.v2_raw.insert(*key, vec![values[sample]]);
    }

    out
}
