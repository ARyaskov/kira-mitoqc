use kira_mitoqc::config::weights::{
    AxisGroupWeights, AxisWeights, DynamicsWeights, Explainability, GlobalWeights, Metadata,
    RegulationWeights, RosWeights, WeightsV1,
};
use kira_mitoqc::core::types::{ProxyKey, ProxyScores};
use kira_mitoqc::score::{compute_axes_v1, compute_decay_v1};

fn weights_fixture() -> WeightsV1 {
    WeightsV1 {
        metadata: Metadata {
            version: "v1.0".to_string(),
            description: "test".to_string(),
            notes: "test".to_string(),
        },
        axis: AxisGroupWeights {
            bioenergetics: AxisWeights {
                etc_stoichiometry_loss: 0.4,
                mtdna_expression_uncoupling: 0.35,
                atp_coupling_loss: 0.25,
            },
            ros: RosWeights {
                ros_response_overdrive: 0.6,
                nadh_imbalance: 0.4,
            },
            dynamics: DynamicsWeights {
                dynamics_imbalance: 0.6,
                mitophagy_excess: 0.4,
            },
            regulation: RegulationWeights {
                biogenesis_failure: 1.0,
            },
        },
        global: GlobalWeights {
            bioenergetics: 0.35,
            ros: 0.3,
            dynamics: 0.2,
            regulation: 0.15,
        },
        explainability: Explainability {
            max_drivers: 8,
            min_abs_contribution: 0.03,
        },
    }
}

fn proxy_scores_fixture() -> ProxyScores {
    let mut scores = ProxyScores::default();
    scores.set(ProxyKey::ETCStoichiometryLoss, vec![0.2, 0.4]);
    scores.set(ProxyKey::MtdnaExpressionUncoupling, vec![0.1, 0.3]);
    scores.set(ProxyKey::AtpCouplingLoss, vec![0.5, 0.7]);

    scores.set(ProxyKey::RosResponseOverdrive, vec![0.2, 0.6]);
    scores.set(ProxyKey::NadhImbalance, vec![0.4, 0.2]);

    scores.set(ProxyKey::DynamicsImbalance, vec![0.3, 0.8]);
    scores.set(ProxyKey::MitophagyExcess, vec![0.5, 0.1]);

    scores.set(ProxyKey::BiogenesisFailure, vec![0.6, 0.2]);

    scores
}

#[test]
fn axis_aggregation_correctness() {
    let weights = weights_fixture();
    let proxies = proxy_scores_fixture();
    let axes = compute_axes_v1(&proxies, &weights);

    let expected_bio = 0.4 * 0.2 + 0.35 * 0.1 + 0.25 * 0.5;
    let expected_ros = 0.6 * 0.2 + 0.4 * 0.4;
    let expected_dyn = 0.6 * 0.3 + 0.4 * 0.5;
    let expected_reg = 1.0 * 0.6;

    assert!((axes.bioenergetics[0] - expected_bio).abs() < 1e-6);
    assert!((axes.ros[0] - expected_ros).abs() < 1e-6);
    assert!((axes.dynamics[0] - expected_dyn).abs() < 1e-6);
    assert!((axes.regulation[0] - expected_reg).abs() < 1e-6);
}

#[test]
fn decay_score_correctness() {
    let weights = weights_fixture();
    let proxies = proxy_scores_fixture();
    let axes = compute_axes_v1(&proxies, &weights);
    let decay = compute_decay_v1(&axes, &weights);

    let expected = 0.35 * axes.bioenergetics[0]
        + 0.3 * axes.ros[0]
        + 0.2 * axes.dynamics[0]
        + 0.15 * axes.regulation[0];
    assert!((decay.decay[0] - expected).abs() < 1e-6);
}

#[test]
fn robustness_margin_bounds() {
    let weights = weights_fixture();
    let mut proxies = proxy_scores_fixture();
    proxies.set(ProxyKey::ETCStoichiometryLoss, vec![1.0, 1.0]);
    proxies.set(ProxyKey::MtdnaExpressionUncoupling, vec![1.0, 1.0]);
    proxies.set(ProxyKey::AtpCouplingLoss, vec![1.0, 1.0]);
    proxies.set(ProxyKey::RosResponseOverdrive, vec![1.0, 1.0]);
    proxies.set(ProxyKey::NadhImbalance, vec![1.0, 1.0]);
    proxies.set(ProxyKey::DynamicsImbalance, vec![1.0, 1.0]);
    proxies.set(ProxyKey::MitophagyExcess, vec![1.0, 1.0]);
    proxies.set(ProxyKey::BiogenesisFailure, vec![1.0, 1.0]);

    let axes = compute_axes_v1(&proxies, &weights);
    let decay = compute_decay_v1(&axes, &weights);

    assert!(decay.decay[0] <= 1.0);
    assert!(decay.robustness_margin[0] >= 0.0);
}

#[test]
fn deterministic_output() {
    let weights = weights_fixture();
    let proxies = proxy_scores_fixture();
    let axes1 = compute_axes_v1(&proxies, &weights);
    let axes2 = compute_axes_v1(&proxies, &weights);
    assert_eq!(axes1, axes2);
}
