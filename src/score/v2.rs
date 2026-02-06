//! v2 axis aggregation with multi-omics corrections.

use tracing::error;

use crate::config::refs_v2::RefsV2;
use crate::config::weights_v2::WeightsV2;
use crate::core::types::{ProxyKey, ProxyScores};
use crate::proxy::{ProxyKeyV2, ProxyScoresV2};
use crate::score::{AxisScoresVec, DecayScoreVec};
use crate::util::numeric::clamp01;

/// Compute v2 axis scores.
pub fn compute_axes_v2(
    proxies: &ProxyScoresV2,
    weights: &WeightsV2,
    refs_v2: &RefsV2,
) -> AxisScoresVec {
    let v1 = &proxies.v1;

    let etc_rna = get_vec(v1, ProxyKey::ETCStoichiometryLoss);
    let uncoupling = get_vec(v1, ProxyKey::MtdnaExpressionUncoupling);
    let atp_rna = get_vec(v1, ProxyKey::AtpCouplingLoss);

    let ros = get_vec(v1, ProxyKey::RosResponseOverdrive);
    let nadh = get_vec(v1, ProxyKey::NadhImbalance);

    let dynamics = get_vec(v1, ProxyKey::DynamicsImbalance);
    let mitophagy = get_vec(v1, ProxyKey::MitophagyExcess);

    let biogenesis = get_vec(v1, ProxyKey::BiogenesisFailure);

    let samples = etc_rna.len();

    let cn = get_v2_or_zero(proxies, ProxyKeyV2::MtDnaCopyNumberInstability, samples);
    let het = get_v2_or_zero(proxies, ProxyKeyV2::MtDnaHeteroplasmyBurden, samples);

    let etc_corrected = mix_proxy(
        etc_rna,
        proxies
            .v2_normalized
            .get(&ProxyKeyV2::ProteomicsEtcStoichiometryLoss),
        refs_v2.mixing.alpha_rna_protein_stoich,
        samples,
    );
    let atp_corrected = mix_proxy(
        atp_rna,
        proxies
            .v2_normalized
            .get(&ProxyKeyV2::ProteomicsAtpCouplingLoss),
        refs_v2.mixing.alpha_rna_protein_atp,
        samples,
    );

    let mut bioenergetics = vec![0.0; samples];
    let mut ros_axis = vec![0.0; samples];
    let mut dynamics_axis = vec![0.0; samples];
    let mut regulation_axis = vec![0.0; samples];

    for i in 0..samples {
        let bio = weights.axis.bioenergetics.etc_stoichiometry_loss * etc_corrected[i]
            + weights.axis.bioenergetics.mtdna_expression_uncoupling * uncoupling[i]
            + weights.axis.bioenergetics.atp_coupling_loss * atp_corrected[i]
            + weights.axis.bioenergetics.mtdna_copy_number_instability * cn[i]
            + weights.axis.bioenergetics.mtdna_heteroplasmy_burden * het[i];
        let ros_val = weights.axis.ros.ros_response_overdrive * ros[i]
            + weights.axis.ros.nadh_imbalance * nadh[i];
        let dyn_val = weights.axis.dynamics.dynamics_imbalance * dynamics[i]
            + weights.axis.dynamics.mitophagy_excess * mitophagy[i];
        let reg_val = weights.axis.regulation.biogenesis_failure * biogenesis[i]
            + weights.axis.regulation.mtdna_copy_number_instability * cn[i]
            + weights.axis.regulation.mtdna_heteroplasmy_burden * het[i];

        if bio.is_nan() || ros_val.is_nan() || dyn_val.is_nan() || reg_val.is_nan() {
            error!(sample = i, "NaN encountered in v2 axis aggregation");
            bioenergetics[i] = 0.0;
            ros_axis[i] = 0.0;
            dynamics_axis[i] = 0.0;
            regulation_axis[i] = 0.0;
        } else {
            bioenergetics[i] = bio;
            ros_axis[i] = ros_val;
            dynamics_axis[i] = dyn_val;
            regulation_axis[i] = reg_val;
        }
    }

    AxisScoresVec {
        bioenergetics,
        ros: ros_axis,
        dynamics: dynamics_axis,
        regulation: regulation_axis,
    }
}

/// Compute decay scores for v2 axes.
pub fn compute_decay_v2(axes: &AxisScoresVec, weights: &WeightsV2) -> DecayScoreVec {
    let len = axes.bioenergetics.len();
    let mut decay = vec![0.0; len];
    let mut robustness_margin = vec![0.0; len];

    for i in 0..len {
        let value = weights.global.bioenergetics * axes.bioenergetics[i]
            + weights.global.ros * axes.ros[i]
            + weights.global.dynamics * axes.dynamics[i]
            + weights.global.regulation * axes.regulation[i];

        if value.is_nan() {
            error!(sample = i, "NaN encountered in v2 decay score");
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

fn get_vec<'a>(proxies: &'a ProxyScores, key: ProxyKey) -> &'a [f32] {
    proxies
        .normalized
        .get(&key)
        .expect("missing proxy values")
        .as_slice()
}

fn get_v2_or_zero(proxies: &ProxyScoresV2, key: ProxyKeyV2, samples: usize) -> Vec<f32> {
    proxies
        .v2_normalized
        .get(&key)
        .cloned()
        .unwrap_or_else(|| vec![0.0; samples])
}

fn mix_proxy(rna: &[f32], protein: Option<&Vec<f32>>, alpha: f32, samples: usize) -> Vec<f32> {
    match protein {
        Some(values) => {
            assert_eq!(values.len(), samples);
            let mut out = Vec::with_capacity(samples);
            for i in 0..samples {
                out.push(alpha * rna[i] + (1.0 - alpha) * values[i]);
            }
            out
        }
        None => rna.to_vec(),
    }
}
