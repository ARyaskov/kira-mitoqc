//! Axis aggregation from proxy scores.

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
pub fn compute_axes_v1(proxies: &ProxyScores, weights: &WeightsV1) -> AxisScoresVec {
    debug_assert!(
        (weights.axis.bioenergetics.etc_stoichiometry_loss
            + weights.axis.bioenergetics.mtdna_expression_uncoupling
            + weights.axis.bioenergetics.atp_coupling_loss
            - 1.0)
            .abs()
            <= 1e-6
    );
    debug_assert!(
        (weights.axis.ros.ros_response_overdrive + weights.axis.ros.nadh_imbalance - 1.0).abs()
            <= 1e-6
    );
    debug_assert!(
        (weights.axis.dynamics.dynamics_imbalance + weights.axis.dynamics.mitophagy_excess - 1.0)
            .abs()
            <= 1e-6
    );
    debug_assert!((weights.axis.regulation.biogenesis_failure - 1.0).abs() <= 1e-6);

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

    let mut bioenergetics = vec![0.0; len];
    let mut ros_axis = vec![0.0; len];
    let mut dynamics_axis = vec![0.0; len];
    let mut regulation_axis = vec![0.0; len];

    for i in 0..len {
        let bio = weights.axis.bioenergetics.etc_stoichiometry_loss * etc[i]
            + weights.axis.bioenergetics.mtdna_expression_uncoupling * uncoupling[i]
            + weights.axis.bioenergetics.atp_coupling_loss * atp[i];
        let ros_val = weights.axis.ros.ros_response_overdrive * ros[i]
            + weights.axis.ros.nadh_imbalance * nadh[i];
        let dyn_val = weights.axis.dynamics.dynamics_imbalance * dynamics[i]
            + weights.axis.dynamics.mitophagy_excess * mitophagy[i];
        let reg_val = weights.axis.regulation.biogenesis_failure * biogenesis[i];

        if bio.is_nan() || ros_val.is_nan() || dyn_val.is_nan() || reg_val.is_nan() {
            error!(sample = i, "NaN encountered in axis aggregation");
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
