use std::fs;
use std::path::PathBuf;
use std::process::Command;

use kira_mitoqc::data::aggregate::AggregationMode;
use kira_mitoqc::input::bd_rhapsody::compute_mito_fraction_from_file;

fn temp_dir(name: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_mitoqc_{name}_{nanos}"));
    fs::create_dir_all(&dir).expect("mkdir");
    dir
}

#[test]
fn run_accepts_bd_rhapsody_and_computes_mito_fraction() {
    let out = temp_dir("bd_run_out");
    let input = out.join("sample_bd_rhapsody.tsv");
    fs::write(
        &input,
        "#Sequencing: BD Rhapsody Whole Transcriptome 3-end sequencing\n\
         #Species: Homo sapiens\n\
         cellA\tcellB\tcellC\n\
         MT-ND1\t2.0\t0.5\t1.0\n\
         MT-ND2\t1.0\t0.5\t1.0\n\
         MT-ND3\t1.0\t0.5\t1.0\n\
         MT-CYB\t2.0\t0.5\t1.0\n\
         MT-CO1\t2.0\t0.5\t1.0\n\
         MT-CO2\t2.0\t0.5\t1.0\n\
         MT-ATP6\t2.0\t0.5\t1.0\n\
         MT-ATP8\t2.0\t0.5\t1.0\n\
         NDUFS1\t10.0\t10.0\t10.0\n\
         SDHA\t10.0\t10.0\t10.0\n\
         UQCRC1\t10.0\t10.0\t10.0\n\
         COX4I1\t10.0\t10.0\t10.0\n\
         COX5A\t10.0\t10.0\t10.0\n\
         ATP5F1A\t10.0\t10.0\t10.0\n\
         SOD2\t5.0\t5.0\t5.0\n\
         GPX1\t5.0\t5.0\t5.0\n\
         GPX4\t5.0\t5.0\t5.0\n\
         PINK1\t5.0\t5.0\t5.0\n\
         MFN1\t5.0\t5.0\t5.0\n\
         DNM1L\t5.0\t5.0\t5.0\n\
         TFAM\t5.0\t5.0\t5.0\n",
    )
    .expect("write input");

    let exe = std::env::var("CARGO_BIN_EXE_kira-mitoqc").expect("kira-mitoqc bin path");
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let run = Command::new(exe)
        .arg("run")
        .arg("--input")
        .arg(&input)
        .arg("--out")
        .arg(&out)
        .arg("--assets")
        .arg(manifest_dir.join("assets"))
        .output()
        .expect("run kira-mitoqc");
    assert!(
        run.status.success(),
        "stdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&run.stdout),
        String::from_utf8_lossy(&run.stderr)
    );

    let _logs = format!(
        "{}\n{}",
        String::from_utf8_lossy(&run.stdout),
        String::from_utf8_lossy(&run.stderr)
    );
    assert!(out.join("mitochondrial_profile.json").exists());

    let mito = std::collections::BTreeSet::from([
        "MT-ND1".to_string(),
        "MT-ND2".to_string(),
        "MT-ND3".to_string(),
        "MT-CYB".to_string(),
        "MT-CO1".to_string(),
        "MT-CO2".to_string(),
        "MT-ATP6".to_string(),
        "MT-ATP8".to_string(),
    ]);
    let fractions =
        compute_mito_fraction_from_file(&input, &mito, AggregationMode::Cell, None).expect("frac");
    assert_eq!(fractions.len(), 3);
    assert!((fractions[0] - 14.0 / 109.0).abs() < 1e-6);
    assert!((fractions[1] - 4.0 / 99.0).abs() < 1e-6);
    assert!((fractions[2] - 8.0 / 103.0).abs() < 1e-6);
}
