//! Axis aggregation from proxy scores.

use rayon::prelude::*;
use tracing::error;

use crate::config::weights::WeightsV1;
use crate::core::types::{ProxyKey, ProxyScores};

/// Axis score vectors.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisScoresVec {
    pub bioenergetics: Vec<f32>,
    pub ros: Vec<f32>,
    pub dynamics: Vec<f32>,
    pub regulation: Vec<f32>,
}

/// Compute axis scores for v1.
///
/// Weight sums are validated at config load time (`WeightsV1::validate`); we
/// don't re-check here on every call.
pub fn compute_axes_v1(proxies: &ProxyScores, weights: &WeightsV1) -> AxisScoresVec {
    let etc = proxies
        .normalized
        .get(&ProxyKey::ETCStoichiometryLoss)
        .expect("missing ETC_stoichiometry_loss");
    let uncoupling = proxies
        .normalized
        .get(&ProxyKey::MtdnaExpressionUncoupling)
        .expect("missing mtDNA_expression_uncoupling");
    let atp = proxies
        .normalized
        .get(&ProxyKey::AtpCouplingLoss)
        .expect("missing ATP_coupling_loss");

    let ros = proxies
        .normalized
        .get(&ProxyKey::RosResponseOverdrive)
        .expect("missing ROS_response_overdrive");
    let nadh = proxies
        .normalized
        .get(&ProxyKey::NadhImbalance)
        .expect("missing NADH_imbalance");

    let dynamics = proxies
        .normalized
        .get(&ProxyKey::DynamicsImbalance)
        .expect("missing dynamics_imbalance");
    let mitophagy = proxies
        .normalized
        .get(&ProxyKey::MitophagyExcess)
        .expect("missing mitophagy_excess");

    let biogenesis = proxies
        .normalized
        .get(&ProxyKey::BiogenesisFailure)
        .expect("missing biogenesis_failure");

    let len = etc.len();
    assert_eq!(uncoupling.len(), len);
    assert_eq!(atp.len(), len);
    assert_eq!(ros.len(), len);
    assert_eq!(nadh.len(), len);
    assert_eq!(dynamics.len(), len);
    assert_eq!(mitophagy.len(), len);
    assert_eq!(biogenesis.len(), len);

    // Hoist weights into Send+Sync locals for rayon.
    let w_bio_etc = weights.axis.bioenergetics.etc_stoichiometry_loss;
    let w_bio_unc = weights.axis.bioenergetics.mtdna_expression_uncoupling;
    let w_bio_atp = weights.axis.bioenergetics.atp_coupling_loss;
    let w_ros_over = weights.axis.ros.ros_response_overdrive;
    let w_ros_nadh = weights.axis.ros.nadh_imbalance;
    let w_dyn_imb = weights.axis.dynamics.dynamics_imbalance;
    let w_dyn_mito = weights.axis.dynamics.mitophagy_excess;
    let w_reg_bio = weights.axis.regulation.biogenesis_failure;

    let mut bioenergetics = vec![0.0; len];
    let mut ros_axis = vec![0.0; len];
    let mut dynamics_axis = vec![0.0; len];
    let mut regulation_axis = vec![0.0; len];

    bioenergetics
        .par_iter_mut()
        .zip(ros_axis.par_iter_mut())
        .zip(dynamics_axis.par_iter_mut())
        .zip(regulation_axis.par_iter_mut())
        .enumerate()
        .with_min_len(1024)
        .for_each(|(i, (((bio_out, ros_out), dyn_out), reg_out))| {
            let bio = w_bio_etc * etc[i] + w_bio_unc * uncoupling[i] + w_bio_atp * atp[i];
            let ros_val = w_ros_over * ros[i] + w_ros_nadh * nadh[i];
            let dyn_val = w_dyn_imb * dynamics[i] + w_dyn_mito * mitophagy[i];
            let reg_val = w_reg_bio * biogenesis[i];

            if bio.is_nan() || ros_val.is_nan() || dyn_val.is_nan() || reg_val.is_nan() {
                error!(sample = i, "NaN encountered in axis aggregation");
                *bio_out = 0.0;
                *ros_out = 0.0;
                *dyn_out = 0.0;
                *reg_out = 0.0;
            } else {
                *bio_out = bio;
                *ros_out = ros_val;
                *dyn_out = dyn_val;
                *reg_out = reg_val;
            }
        });

    AxisScoresVec {
        bioenergetics,
        ros: ros_axis,
        dynamics: dynamics_axis,
        regulation: regulation_axis,
    }
}
