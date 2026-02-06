use std::fs;

use kira_mitoqc::core::types::MitochondrialState;
use kira_mitoqc::core::types::{AxisScores, ProxyKey, ProxyScores};
use kira_mitoqc::explain::{AxisKind, Driver, Explainability};
use kira_mitoqc::output::profile::{MitoProfileV1, assemble_profiles_v1};
use kira_mitoqc::output::{write_axes_tsv, write_decay_tsv, write_json, write_proxies_tsv};
use kira_mitoqc::score::{AxisScoresVec, DecayScoreVec};

fn temp_dir(name: &str) -> std::path::PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_mitoqc_out_{name}_{nanos}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn proxy_scores_fixture() -> ProxyScores {
    let mut scores = ProxyScores::default();
    scores.set(ProxyKey::ETCStoichiometryLoss, vec![0.1, 0.2]);
    scores.set(ProxyKey::MtdnaExpressionUncoupling, vec![0.3, 0.4]);
    scores.set(ProxyKey::AtpCouplingLoss, vec![0.5, 0.6]);
    scores.set(ProxyKey::RosResponseOverdrive, vec![0.2, 0.3]);
    scores.set(ProxyKey::NadhImbalance, vec![0.4, 0.1]);
    scores.set(ProxyKey::DynamicsImbalance, vec![0.6, 0.2]);
    scores.set(ProxyKey::MitophagyExcess, vec![0.1, 0.2]);
    scores.set(ProxyKey::BiogenesisFailure, vec![0.7, 0.8]);

    scores.set_raw(ProxyKey::ETCStoichiometryLoss, vec![0.01, 0.02]);
    scores.set_raw(ProxyKey::MtdnaExpressionUncoupling, vec![0.03, 0.04]);
    scores
}

#[test]
fn assemble_profiles_wires_fields() {
    let states = vec![
        MitochondrialState::BioenergeticCollapse,
        MitochondrialState::CompensatedButFragile,
    ];
    let decay = DecayScoreVec {
        decay: vec![0.7, 0.2],
        robustness_margin: vec![0.3, 0.8],
    };
    let axes = AxisScoresVec {
        bioenergetics: vec![0.6, 0.1],
        ros: vec![0.4, 0.2],
        dynamics: vec![0.3, 0.2],
        regulation: vec![0.2, 0.1],
    };
    let proxies = proxy_scores_fixture();
    let explain = vec![
        Explainability {
            drivers: vec![Driver {
                key: ProxyKey::ETCStoichiometryLoss,
                axis: AxisKind::Bioenergetics,
                contribution: 0.1,
            }],
            interpretation: vec!["A".to_string()],
        },
        Explainability {
            drivers: vec![],
            interpretation: vec!["B".to_string()],
        },
    ];

    let profiles = assemble_profiles_v1(&states, &decay, &axes, &proxies, &explain);

    assert_eq!(profiles.len(), 2);
    assert_eq!(profiles[0].mitochondrial_state, "Bioenergetic collapse");
    assert_eq!(profiles[0].decay_score, 0.7);
    assert_eq!(profiles[0].axes.bioenergetics, 0.6);
}

#[test]
fn json_roundtrip() {
    let profile = MitoProfileV1 {
        mitochondrial_state: "Compensated but fragile".to_string(),
        decay_score: 0.2,
        robustness_margin: 0.8,
        axes: AxisScores {
            bioenergetics: 0.1,
            ros: 0.2,
            dynamics: 0.3,
            regulation: 0.4,
        },
        proxies: proxy_scores_fixture(),
        drivers: vec![],
        interpretation: vec!["ok".to_string()],
    };
    let dir = temp_dir("json");
    write_json(&dir, &[profile.clone()]).expect("write json");

    let data = fs::read_to_string(dir.join("mitochondrial_profile.json")).unwrap();
    let decoded: Vec<MitoProfileV1> = serde_json::from_str(&data).unwrap();
    assert_eq!(decoded[0], profile);
}

#[test]
fn tsv_headers_are_stable() {
    let dir = temp_dir("tsv");
    let axes = AxisScoresVec {
        bioenergetics: vec![0.1],
        ros: vec![0.2],
        dynamics: vec![0.3],
        regulation: vec![0.4],
    };
    let decay = DecayScoreVec {
        decay: vec![0.5],
        robustness_margin: vec![0.5],
    };
    let proxies = proxy_scores_fixture();

    write_axes_tsv(&dir, &axes).expect("axes tsv");
    write_decay_tsv(&dir, &decay).expect("decay tsv");
    write_proxies_tsv(&dir, &proxies).expect("proxies tsv");

    let axes_header = fs::read_to_string(dir.join("axes.tsv")).unwrap();
    let decay_header = fs::read_to_string(dir.join("decay.tsv")).unwrap();
    let proxies_header = fs::read_to_string(dir.join("proxies.tsv")).unwrap();

    assert!(
        axes_header
            .lines()
            .next()
            .unwrap()
            .starts_with("sample\tbioenergetics\tros\tdynamics\tregulation")
    );
    assert!(
        decay_header
            .lines()
            .next()
            .unwrap()
            .starts_with("sample\tdecay_score\trobustness_margin")
    );
    assert!(
        proxies_header
            .lines()
            .next()
            .unwrap()
            .starts_with("sample\tETC_stoichiometry_loss")
    );
}
