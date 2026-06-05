//! v2 axis aggregation with multi-omics corrections.

use rayon::prelude::*;
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

    let w_bio = &weights.axis.bioenergetics;
    let w_ros = &weights.axis.ros;
    let w_dyn = &weights.axis.dynamics;
    let w_reg = &weights.axis.regulation;

    let (w_bio_etc, w_bio_unc, w_bio_atp, w_bio_cn, w_bio_het) = (
        w_bio.etc_stoichiometry_loss,
        w_bio.mtdna_expression_uncoupling,
        w_bio.atp_coupling_loss,
        w_bio.mtdna_copy_number_instability,
        w_bio.mtdna_heteroplasmy_burden,
    );
    let (w_ros_over, w_ros_nadh) = (w_ros.ros_response_overdrive, w_ros.nadh_imbalance);
    let (w_dyn_imb, w_dyn_mito) = (w_dyn.dynamics_imbalance, w_dyn.mitophagy_excess);
    let (w_reg_bio, w_reg_cn, w_reg_het) = (
        w_reg.biogenesis_failure,
        w_reg.mtdna_copy_number_instability,
        w_reg.mtdna_heteroplasmy_burden,
    );

    let mut bioenergetics = vec![0.0; samples];
    let mut ros_axis = vec![0.0; samples];
    let mut dynamics_axis = vec![0.0; samples];
    let mut regulation_axis = vec![0.0; samples];

    bioenergetics
        .par_iter_mut()
        .zip(ros_axis.par_iter_mut())
        .zip(dynamics_axis.par_iter_mut())
        .zip(regulation_axis.par_iter_mut())
        .enumerate()
        .with_min_len(1024)
        .for_each(|(i, (((bio_o, ros_o), dyn_o), reg_o))| {
            let bio = w_bio_etc * etc_corrected[i]
                + w_bio_unc * uncoupling[i]
                + w_bio_atp * atp_corrected[i]
                + w_bio_cn * cn[i]
                + w_bio_het * het[i];
            let ros_val = w_ros_over * ros[i] + w_ros_nadh * nadh[i];
            let dyn_val = w_dyn_imb * dynamics[i] + w_dyn_mito * mitophagy[i];
            let reg_val = w_reg_bio * biogenesis[i] + w_reg_cn * cn[i] + w_reg_het * het[i];

            if bio.is_nan() || ros_val.is_nan() || dyn_val.is_nan() || reg_val.is_nan() {
                error!(sample = i, "NaN encountered in v2 axis aggregation");
                *bio_o = 0.0;
                *ros_o = 0.0;
                *dyn_o = 0.0;
                *reg_o = 0.0;
            } else {
                *bio_o = bio;
                *ros_o = ros_val;
                *dyn_o = dyn_val;
                *reg_o = reg_val;
            }
        });

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
    let w_bio = weights.global.bioenergetics;
    let w_ros = weights.global.ros;
    let w_dyn = weights.global.dynamics;
    let w_reg = weights.global.regulation;

    let mut decay = vec![0.0; len];
    let mut robustness_margin = vec![0.0; len];

    decay
        .par_iter_mut()
        .zip(robustness_margin.par_iter_mut())
        .enumerate()
        .with_min_len(1024)
        .for_each(|(i, (d, r))| {
            let value = w_bio * axes.bioenergetics[i]
                + w_ros * axes.ros[i]
                + w_dyn * axes.dynamics[i]
                + w_reg * axes.regulation[i];

            if value.is_nan() {
                error!(sample = i, "NaN encountered in v2 decay score");
                *d = 0.0;
                *r = 0.0;
            } else {
                let clamped = clamp01(value);
                *d = clamped;
                *r = clamp01(1.0 - clamped);
            }
        });

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
