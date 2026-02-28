use std::fs;
use std::path::PathBuf;

use kira_mitoqc::input::{DetectedInputFormat, InputFormat, detect_input_format};

fn temp_path(name: &str, content: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("kira_mitoqc_{name}_{nanos}.txt"));
    fs::write(&path, content).expect("write");
    path
}

#[test]
fn auto_detects_bd_rhapsody_for_comment_and_float_shape() {
    let file = temp_path(
        "bd_detect",
        "#meta\n\
         cellA\tcellB\n\
         GENE1\t1.0\t2.5\n",
    );
    let detected = detect_input_format(&file, InputFormat::Auto).expect("detect");
    assert_eq!(detected, DetectedInputFormat::BDRhapsodyDense);
}

#[test]
fn auto_defaults_to_tenx_when_not_matching_bd_signature() {
    let file = temp_path(
        "not_bd",
        "cellA\tcellB\n\
         GENE1\t1\t2\n",
    );
    let detected = detect_input_format(&file, InputFormat::Auto).expect("detect");
    assert_eq!(detected, DetectedInputFormat::Tenx);
}

#[test]
fn explicit_override_bd_rhapsody_is_honored() {
    let file = temp_path(
        "bd_override",
        "cellA\tcellB\n\
         GENE1\t1\t2\n",
    );
    let detected = detect_input_format(&file, InputFormat::BdRhapsody).expect("detect");
    assert_eq!(detected, DetectedInputFormat::BDRhapsodyDense);
}

#[test]
fn auto_detects_bd_rhapsody_for_raw_counts_filename() {
    let file = temp_path("sample_raw_counts.tsv", "gene\tcellA\tcellB\nGENE1\t1\t2\n");
    let detected = detect_input_format(&file, InputFormat::Auto).expect("detect");
    assert_eq!(detected, DetectedInputFormat::BDRhapsodyDense);
}

#[test]
fn auto_detects_bd_rhapsody_for_directory_with_prefixed_raw_counts() {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_mitoqc_bd_dir_{nanos}"));
    fs::create_dir_all(&dir).expect("mkdir");
    fs::write(
        dir.join("GSM3635278_human_p1t1_raw_counts.tsv"),
        "gene\tcellA\tcellB\nGENE1\t1\t2\n",
    )
    .expect("write");

    let detected = detect_input_format(&dir, InputFormat::Auto).expect("detect");
    assert_eq!(detected, DetectedInputFormat::BDRhapsodyDense);
}
