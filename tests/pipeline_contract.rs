use std::fs;

use kira_mitoqc::core::types::{AxisScores, ProxyScores};
use kira_mitoqc::input::ExpressionSource;
use kira_mitoqc::metrics::metabolic_extension::aggregate::build_summary as build_metabolic_summary;
use kira_mitoqc::metrics::metabolic_extension::scores::{
    MetabolicMetrics, MetabolicMissingness, MetabolicThresholds,
};
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

fn sample_metabolic_metrics() -> MetabolicMetrics {
    MetabolicMetrics {
        oxphos_core: vec![0.11, 0.22],
        gly_core: vec![0.33, 0.44],
        fao_core: vec![0.55, 0.66],
        ros_core: vec![0.77, 0.88],
        bio_core: vec![0.99, 0.12],
        mri: vec![1.1, 2.2],
        osl: vec![1.3, 2.4],
        ess: vec![0.5, 1.6],
        mcb: vec![-0.5, -1.6],
        ogi: vec![1.7, -1.8],
        metabolic_rigid_high: vec![false, true],
        ros_high: vec![false, true],
        energetic_strain_high: vec![false, true],
        compensation_failure: vec![false, true],
        oxphos_dominant: vec![true, false],
        glycolysis_dominant: vec![false, true],
        thresholds: MetabolicThresholds::default(),
        missingness: MetabolicMissingness {
            oxphos_found: 10,
            oxphos_total: 10,
            glycolysis_found: 10,
            glycolysis_total: 10,
            fao_found: 6,
            fao_total: 6,
            ros_found: 8,
            ros_total: 13,
            biogenesis_found: 6,
            biogenesis_total: 7,
        },
    }
}

#[test]
fn writes_summary_with_required_fields() {
    let dir = temp_dir("summary");
    let profiles = sample_profiles();

    write_summary_json(&dir, &profiles, "10x", ExpressionSource::RawUmiCounts).expect("summary");
    let value: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("summary.json")).unwrap()).unwrap();

    assert_eq!(value["tool"], "kira-mitoqc");
    assert_eq!(value["input"]["mode"], "pipeline");
    assert_eq!(value["input"]["n_samples"], 2);
    assert_eq!(value["input"]["input_format"], "10x");
    assert_eq!(value["input"]["expression_type"], "raw_umi_counts");
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
    let metabolic = sample_metabolic_metrics();

    let barcodes = vec!["cellA".to_string(), "cellB".to_string()];
    write_mito_metrics_tsv(&dir, &barcodes, &profiles, &metabolic).expect("metrics");
    let contents = fs::read_to_string(dir.join("mito_metrics.tsv")).unwrap();
    let mut lines = contents.lines();
    assert_eq!(
        lines.next().unwrap(),
        "cell_id\tmitochondrial_state\tdecay_score\trobustness_margin\tbioenergetics\tros\tdynamics\tregulation\toxphos_core\tgly_core\tfao_core\tros_core\tbio_core\tMRI\tOSL\tESS\tMCB\tOGI\tmetabolic_rigid_high\tros_high\tenergetic_strain_high\tcompensation_failure"
    );
    assert_eq!(
        lines.next().unwrap(),
        "cellA\tBioenergetic collapse\t0.900000\t0.100000\t0.800000\t0.500000\t0.400000\t0.300000\t0.110000\t0.330000\t0.550000\t0.770000\t0.990000\t1.100000\t1.300000\t0.500000\t-0.500000\t1.700000\tfalse\tfalse\tfalse\tfalse"
    );
    assert_eq!(
        lines.next().unwrap(),
        "cellB\tCompensated but fragile\t0.400000\t0.600000\t0.200000\t0.300000\t0.700000\t0.500000\t0.220000\t0.440000\t0.660000\t0.880000\t0.120000\t2.200000\t2.400000\t1.600000\t-1.600000\t-1.800000\ttrue\ttrue\ttrue\ttrue"
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

    write_summary_json(&dir, &profiles, "10x", ExpressionSource::RawUmiCounts)
        .expect("first write");
    let first = fs::read(dir.join("summary.json")).unwrap();
    write_summary_json(&dir, &profiles, "10x", ExpressionSource::RawUmiCounts)
        .expect("second write");
    let second = fs::read(dir.join("summary.json")).unwrap();
    assert_eq!(first, second);
}

#[test]
fn summary_unchanged_when_redox_disabled() {
    let dir = temp_dir("disabled_redox");
    let profiles = sample_profiles();

    write_summary_json(&dir, &profiles, "10x", ExpressionSource::RawUmiCounts)
        .expect("base summary");
    let base = fs::read(dir.join("summary.json")).expect("read base");

    kira_mitoqc::output::pipeline_contract::write_summary_json_with_redox(
        &dir,
        &profiles,
        "10x",
        ExpressionSource::RawUmiCounts,
        None,
        None,
    )
    .expect("summary with disabled redox");
    let disabled = fs::read(dir.join("summary.json")).expect("read disabled");

    assert_eq!(base, disabled);
}

#[test]
fn summary_can_include_mitochondrial_metabolic_block() {
    let dir = temp_dir("metabolic_summary");
    let profiles = sample_profiles();
    let metrics = sample_metabolic_metrics();
    let sample_ids = vec!["clusterA".to_string(), "clusterB".to_string()];
    let metabolic_summary = build_metabolic_summary(&metrics, &sample_ids, true);

    kira_mitoqc::output::pipeline_contract::write_summary_json_with_redox(
        &dir,
        &profiles,
        "10x",
        ExpressionSource::RawUmiCounts,
        None,
        Some(&metabolic_summary),
    )
    .expect("summary with metabolic block");

    let value: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(dir.join("summary.json")).unwrap()).unwrap();
    assert_eq!(
        value["mitochondrial_metabolic"]["panel_version"],
        "MITO_METABOLIC_PANEL_V1"
    );
    assert!(value["mitochondrial_metabolic"]["thresholds"]["metabolic_rigid_high"].is_number());
    assert!(value["mitochondrial_metabolic"]["global_stats"]["mri_median"].is_number());
    assert!(value["mitochondrial_metabolic"]["cluster_stats"].is_array());
}
