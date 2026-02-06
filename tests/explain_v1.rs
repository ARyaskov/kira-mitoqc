use kira_mitoqc::config::weights::{
    AxisGroupWeights, AxisWeights, DynamicsWeights, Explainability, GlobalWeights, Metadata,
    RegulationWeights, RosWeights, WeightsV1,
};
use kira_mitoqc::core::types::{MitochondrialState, ProxyKey, ProxyScores};
use kira_mitoqc::explain::{AxisKind, explain_v1};
use kira_mitoqc::score::{AxisScoresVec, DecayScoreVec};

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
            max_drivers: 3,
            min_abs_contribution: 0.05,
        },
    }
}

fn proxies_fixture() -> ProxyScores {
    let mut scores = ProxyScores::default();
    scores.set(ProxyKey::ETCStoichiometryLoss, vec![0.9]);
    scores.set(ProxyKey::MtdnaExpressionUncoupling, vec![0.2]);
    scores.set(ProxyKey::AtpCouplingLoss, vec![0.1]);
    scores.set(ProxyKey::RosResponseOverdrive, vec![0.3]);
    scores.set(ProxyKey::NadhImbalance, vec![0.7]);
    scores.set(ProxyKey::DynamicsImbalance, vec![0.4]);
    scores.set(ProxyKey::MitophagyExcess, vec![0.05]);
    scores.set(ProxyKey::BiogenesisFailure, vec![0.8]);
    scores
}

fn axes_fixture() -> AxisScoresVec {
    AxisScoresVec {
        bioenergetics: vec![0.6],
        ros: vec![0.7],
        dynamics: vec![0.2],
        regulation: vec![0.6],
    }
}

#[test]
fn contribution_computation_and_driver_sorting() {
    let weights = weights_fixture();
    let proxies = proxies_fixture();
    let axes = axes_fixture();
    let decay = DecayScoreVec {
        decay: vec![0.0],
        robustness_margin: vec![1.0],
    };
    let states = vec![MitochondrialState::BioenergeticCollapse];

    let expl = explain_v1(&proxies, &axes, &decay, &states, &weights);
    let drivers = &expl[0].drivers;

    assert!(drivers.len() <= 3);
    assert_eq!(drivers[0].key, ProxyKey::ETCStoichiometryLoss);
    assert_eq!(drivers[0].axis, AxisKind::Bioenergetics);
}

#[test]
fn min_abs_contribution_filters() {
    let weights = weights_fixture();
    let mut proxies = proxies_fixture();
    proxies.set(ProxyKey::MitophagyExcess, vec![0.01]);
    let axes = axes_fixture();
    let decay = DecayScoreVec {
        decay: vec![0.0],
        robustness_margin: vec![1.0],
    };
    let states = vec![MitochondrialState::CompensatedButFragile];

    let expl = explain_v1(&proxies, &axes, &decay, &states, &weights);
    assert!(expl[0].drivers.iter().all(|d| d.contribution.abs() >= 0.05));
}

#[test]
fn interpretation_strings_include_axes_and_state() {
    let weights = weights_fixture();
    let proxies = proxies_fixture();
    let axes = axes_fixture();
    let decay = DecayScoreVec {
        decay: vec![0.0],
        robustness_margin: vec![1.0],
    };
    let states = vec![MitochondrialState::RosDominantDecay];

    let expl = explain_v1(&proxies, &axes, &decay, &states, &weights);
    let lines = &expl[0].interpretation;
    assert!(lines.iter().any(|l| l.contains("Bioenergetic instability")));
    assert!(
        lines
            .iter()
            .any(|l| l.contains("Oxidative stress response"))
    );
    assert!(
        lines
            .iter()
            .any(|l| l.contains("Insufficient mitochondrial biogenesis"))
    );
    assert!(
        lines
            .iter()
            .any(|l| l.contains("ROS-driven mitochondrial decay"))
    );
}
