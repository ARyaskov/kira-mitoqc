use flate2::Compression;
use flate2::write::GzEncoder;
use kira_mitoqc::io::mtx::{load_mtx_dir, load_mtx_metadata};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn temp_dir(name: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_mitoqc_legacy_{name}_{nanos}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_gz(path: &std::path::Path, contents: &str) {
    let file = File::create(path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(contents.as_bytes()).unwrap();
    encoder.finish().unwrap();
}

fn write_basic_matrix(dir: &PathBuf, rows: usize, cols: usize) {
    fs::write(
        dir.join("matrix.mtx"),
        format!("%%MatrixMarket matrix coordinate real general\n{rows} {cols} 1\n1 1 1.0\n"),
    )
    .unwrap();
    fs::write(dir.join("barcodes.tsv"), "C1\n").unwrap();
}

#[test]
fn legacy_genes_tsv_used_when_features_missing() {
    let dir = temp_dir("genes_only");
    write_basic_matrix(&dir, 2, 1);
    fs::write(dir.join("genes.tsv"), "ENSG1\tMT-ND1\nENSG2\tNDUFS1\n").unwrap();

    let (features, _) = load_mtx_metadata(&dir, None).expect("metadata");
    assert_eq!(features, vec!["MT-ND1".to_string(), "NDUFS1".to_string()]);

    let mtx = load_mtx_dir(&dir, None).expect("mtx load");
    assert_eq!(mtx.features.len(), 2);
}

#[test]
fn features_tsv_takes_priority_over_genes_tsv() {
    let dir = temp_dir("mixed");
    write_basic_matrix(&dir, 1, 1);
    fs::write(dir.join("features.tsv"), "SYM1\tX\n").unwrap();
    fs::write(dir.join("genes.tsv"), "ENSG1\tSYM2\n").unwrap();

    let (features, _) = load_mtx_metadata(&dir, Some(0)).expect("metadata");
    assert_eq!(features, vec!["SYM1".to_string()]);
}

#[test]
fn legacy_genes_tsv_gz_supported() {
    let dir = temp_dir("genes_gz");
    write_basic_matrix(&dir, 2, 1);
    write_gz(&dir.join("genes.tsv.gz"), "ENSG1\tMT-CO1\nENSG2\tATP5F1A\n");

    let (features, _) = load_mtx_metadata(&dir, None).expect("metadata");
    assert_eq!(features, vec!["MT-CO1".to_string(), "ATP5F1A".to_string()]);
}

#[test]
fn invalid_genes_tsv_errors() {
    let dir = temp_dir("genes_invalid");
    write_basic_matrix(&dir, 1, 1);
    fs::write(dir.join("genes.tsv"), "ONLYONECOLUMN\n").unwrap();

    let err = load_mtx_metadata(&dir, None).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("genes.tsv"));
}
