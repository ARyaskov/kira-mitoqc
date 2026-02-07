use std::fs;

use kira_mitoqc::core::types::{AxisScores, ProxyScores};
use kira_mitoqc::output::pipeline_contract::{
    write_mito_metrics_tsv, write_pipeline_step_json, write_summary_json,
};
use kira_mitoqc::output::profile::MitoProfileV1;

fn temp_dir(name: &str) -> std::path::PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_mitoqc_contract_{name}_{nanos}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn sample_profiles() -> Vec<MitoProfileV1> {
    vec![
        MitoProfileV1 {
            mitochondrial_state: "Bioenergetic collapse".to_string(),
            decay_score: 0.9,
            robustness_margin: 0.1,
            axes: AxisScores {
                bioenergetics: 0.8,
                ros: 0.5,
                dynamics: 0.4,
                regulation: 0.3,
            },
            proxies: ProxyScores::default(),
            drivers: vec![],
            interpretation: vec![],
        },
        MitoProfileV1 {
            mitochondrial_state: "Compensated but fragile".to_string(),
            decay_score: 0.4,
            robustness_margin: 0.6,
            axes: AxisScores {
                bioenergetics: 0.2,
                ros: 0.3,
                dynamics: 0.7,
                regulation: 0.5,
            },
            proxies: ProxyScores::default(),
            drivers: vec![],
            interpretation: vec![],
        },
    ]
}

#[test]
fn writes_summary_with_required_fields() {
    let dir = temp_dir("summary");
    let profiles = sample_profiles();

    write_summary_json(&dir, &profiles).expect("summary");
    let value: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("summary.json")).unwrap()).unwrap();

    assert_eq!(value["tool"], "kira-mitoqc");
    assert_eq!(value["input"]["mode"], "pipeline");
    assert_eq!(value["input"]["n_samples"], 2);
    assert!(value["mitochondrial_state_distribution"].is_object());
    assert!(value["decay"]["decay_score_median"].is_number());
    assert!(value["decay"]["robustness_margin_median"].is_number());
    assert!(value["axes_median"]["bioenergetics"].is_number());
    assert!(value["axes_median"]["ros"].is_number());
    assert!(value["axes_median"]["dynamics"].is_number());
    assert!(value["axes_median"]["regulation"].is_number());
}

#[test]
fn writes_metrics_tsv_with_stable_columns() {
    let dir = temp_dir("metrics");
    let profiles = sample_profiles();

    let barcodes = vec!["cellA".to_string(), "cellB".to_string()];
    write_mito_metrics_tsv(&dir, &barcodes, &profiles).expect("metrics");
    let contents = fs::read_to_string(dir.join("mito_metrics.tsv")).unwrap();
    let mut lines = contents.lines();
    assert_eq!(
        lines.next().unwrap(),
        "cell_id\tmitochondrial_state\tdecay_score\trobustness_margin\tbioenergetics\tros\tdynamics\tregulation"
    );
    assert_eq!(
        lines.next().unwrap(),
        "cellA\tBioenergetic collapse\t0.900000\t0.100000\t0.800000\t0.500000\t0.400000\t0.300000"
    );
    assert_eq!(
        lines.next().unwrap(),
        "cellB\tCompensated but fragile\t0.400000\t0.600000\t0.200000\t0.300000\t0.700000\t0.500000"
    );
}

#[test]
fn writes_pipeline_step_manifest() {
    let dir = temp_dir("manifest");

    write_pipeline_step_json(&dir, "kira-organelle.bin").expect("manifest");
    let value: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("pipeline_step.json")).unwrap()).unwrap();

    assert_eq!(value["tool"], "kira-mitoqc");
    assert_eq!(value["mode"], "pipeline");
    assert_eq!(value["artifacts"]["summary"], "summary.json");
    assert_eq!(value["artifacts"]["primary_metrics"], "mito_metrics.tsv");
    assert_eq!(value["artifacts"]["shared_cache"], "kira-organelle.bin");
    assert_eq!(value["sample_metrics"]["file"], "mito_metrics.tsv");
    assert_eq!(value["sample_metrics"]["id_column"], "cell_id");
    assert_eq!(
        value["sample_metrics"]["state_column"],
        "mitochondrial_state"
    );
}

#[test]
fn summary_output_is_stable_on_repeat_write() {
    let dir = temp_dir("stable");
    let profiles = sample_profiles();

    write_summary_json(&dir, &profiles).expect("first write");
    let first = fs::read(dir.join("summary.json")).unwrap();
    write_summary_json(&dir, &profiles).expect("second write");
    let second = fs::read(dir.join("summary.json")).unwrap();
    assert_eq!(first, second);
}
