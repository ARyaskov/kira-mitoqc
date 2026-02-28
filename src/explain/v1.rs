//! v1 explainability implementation.

use tracing::error;

use crate::config::weights::WeightsV1;
use crate::core::types::{MitochondrialState, ProxyKey, ProxyScores};
use crate::explain::{AxisKind, Driver, Explainability};
use crate::score::{AxisScoresVec, DecayScoreVec};

const AXIS_THRESHOLD: f32 = 0.5;

/// Compute explainability for v1 outputs.
pub fn explain_v1(
    proxies: &ProxyScores,
    axes: &AxisScoresVec,
    _decay: &DecayScoreVec,
    states: &[MitochondrialState],
    weights: &WeightsV1,
) -> Vec<Explainability> {
    let len = axes.bioenergetics.len();
    assert_eq!(axes.ros.len(), len);
    assert_eq!(axes.dynamics.len(), len);
    assert_eq!(axes.regulation.len(), len);
    assert_eq!(states.len(), len);

    let proxy_order = [
        ProxyKey::ETCStoichiometryLoss,
        ProxyKey::MtdnaExpressionUncoupling,
        ProxyKey::AtpCouplingLoss,
        ProxyKey::RosResponseOverdrive,
        ProxyKey::NadhImbalance,
        ProxyKey::DynamicsImbalance,
        ProxyKey::MitophagyExcess,
        ProxyKey::BiogenesisFailure,
    ];

    let mut outputs = Vec::with_capacity(len);
    for i in 0..len {
        let drivers = compute_drivers_for_sample(i, proxies, weights, &proxy_order);
        let interpretation = build_interpretation(i, axes, states[i]);
        outputs.push(Explainability {
            drivers,
            interpretation,
        });
    }

    outputs
}

fn compute_drivers_for_sample(
    sample: usize,
    proxies: &ProxyScores,
    weights: &WeightsV1,
    proxy_order: &[ProxyKey],
) -> Vec<Driver> {
    let mut drivers = Vec::with_capacity(proxy_order.len());

    for key in proxy_order {
        let values = proxies.normalized.get(key).expect("missing proxy values");
        let value = values[sample];
        let (axis, axis_weight, global_weight) = proxy_weights(*key, weights);
        let mut contribution = value * axis_weight * global_weight;
        if contribution.is_nan() {
            error!(sample, proxy = key.as_str(), "NaN contribution");
            contribution = 0.0;
        }
        drivers.push(Driver {
            key: *key,
            axis,
            contribution,
        });
    }

    drivers.sort_by(|a, b| {
        let abs_a = a.contribution.abs();
        let abs_b = b.contribution.abs();
        abs_b
            .partial_cmp(&abs_a)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.key.cmp(&b.key))
    });

    let max_drivers = weights.explainability.max_drivers;
    let min_abs = weights.explainability.min_abs_contribution;

    drivers
        .into_iter()
        .take(max_drivers)
        .filter(|d| d.contribution.abs() >= min_abs)
        .collect()
}

fn proxy_weights(key: ProxyKey, weights: &WeightsV1) -> (AxisKind, f32, f32) {
    match key {
        ProxyKey::ETCStoichiometryLoss => (
            AxisKind::Bioenergetics,
            weights.axis.bioenergetics.etc_stoichiometry_loss,
            weights.global.bioenergetics,
        ),
        ProxyKey::MtdnaExpressionUncoupling => (
            AxisKind::Bioenergetics,
            weights.axis.bioenergetics.mtdna_expression_uncoupling,
            weights.global.bioenergetics,
        ),
        ProxyKey::AtpCouplingLoss => (
            AxisKind::Bioenergetics,
            weights.axis.bioenergetics.atp_coupling_loss,
            weights.global.bioenergetics,
        ),
        ProxyKey::RosResponseOverdrive => (
            AxisKind::Ros,
            weights.axis.ros.ros_response_overdrive,
            weights.global.ros,
        ),
        ProxyKey::NadhImbalance => (
            AxisKind::Ros,
            weights.axis.ros.nadh_imbalance,
            weights.global.ros,
        ),
        ProxyKey::DynamicsImbalance => (
            AxisKind::Dynamics,
            weights.axis.dynamics.dynamics_imbalance,
            weights.global.dynamics,
        ),
        ProxyKey::MitophagyExcess => (
            AxisKind::Dynamics,
            weights.axis.dynamics.mitophagy_excess,
            weights.global.dynamics,
        ),
        ProxyKey::BiogenesisFailure => (
            AxisKind::Regulation,
            weights.axis.regulation.biogenesis_failure,
            weights.global.regulation,
        ),
    }
}

fn build_interpretation(
    sample: usize,
    axes: &AxisScoresVec,
    state: MitochondrialState,
) -> Vec<String> {
    let mut lines = Vec::new();

    let bio = axes.bioenergetics[sample];
    let ros = axes.ros[sample];
    let dyns = axes.dynamics[sample];
    let reg = axes.regulation[sample];

    if bio > AXIS_THRESHOLD {
        lines
            .push("Bioenergetic instability driven by ETC imbalance and coupling loss".to_string());
    }
    if ros > AXIS_THRESHOLD {
        lines.push("Oxidative stress response exceeds compensatory capacity".to_string());
    }
    if dyns > AXIS_THRESHOLD {
        lines.push(
            "Mitochondrial dynamics shifted toward fragmentation or excessive turnover".to_string(),
        );
    }
    if reg > AXIS_THRESHOLD {
        lines.push("Insufficient mitochondrial biogenesis response".to_string());
    }

    let state_line = match state {
        MitochondrialState::RosDominantDecay => {
            "ROS-driven mitochondrial decay with preserved core bioenergetics"
        }
        MitochondrialState::BioenergeticCollapse => "Primary failure of energy production capacity",
        MitochondrialState::MitophagyLockedDepletion => {
            "Excessive mitophagy without sufficient biogenesis compensation"
        }
        MitochondrialState::StructuralFragmentation => {
            "Concurrent bioenergetic stress and structural instability"
        }
        MitochondrialState::CompensatedButFragile => {
            "Mitochondrial function maintained by compensatory mechanisms"
        }
        MitochondrialState::CompensatedOxidativeStress => {
            "Oxidative stress proxy is elevated but buffering programs remain active"
        }
        MitochondrialState::UnbufferedOxidativeStress => {
            "Oxidative stress proxy exceeds current redox buffering capacity"
        }
        MitochondrialState::RedoxOverload => {
            "Persistent oxidative stress proxy overload with maladaptive buffering"
        }
    };
    lines.push(state_line.to_string());

    lines
}
