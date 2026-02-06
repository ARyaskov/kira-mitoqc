use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

use kira_mitoqc::core::types::GeneSet;
use kira_mitoqc::data::aggregate::{AggregationMode, ClusterMap, aggregate};
use kira_mitoqc::data::prepare_expression_with_clusters;
use kira_mitoqc::input::{GeneIndex, resolve_all_genesets};
use kira_mitoqc::io::mtx::{MtxInput, load_mtx_dir};
use sprs::CsMat;

fn temp_dir() -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir =
        std::env::temp_dir().join(format!("kira_mitoqc_test_{}_{}", std::process::id(), nanos));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_file(path: &PathBuf, content: &str) {
    let mut file = File::create(path).unwrap();
    file.write_all(content.as_bytes()).unwrap();
}

#[test]
fn load_mtx_dir_reads_matrix_and_metadata() {
    let dir = temp_dir();
    write_file(
        &dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate real general\n3 3 3\n1 1 1.0\n2 2 2.0\n3 1 3.0\n",
    );
    write_file(&dir.join("features.tsv"), "G1\tX\nG2\tY\nG3\tZ\n");
    write_file(&dir.join("barcodes.tsv"), "C1\nC2\nC3\n");

    let mtx = load_mtx_dir(&dir, Some(0)).expect("mtx load");
    assert_eq!(mtx.features.len(), 3);
    assert_eq!(mtx.barcodes.len(), 3);
    assert_eq!(mtx.matrix.rows(), 3);
    assert_eq!(mtx.matrix.cols(), 3);
}

fn sample_matrix() -> CsMat<f32> {
    let rows = 2;
    let cols = 3;
    let indptr = vec![0, 2, 3, 4];
    let indices = vec![0, 1, 0, 1];
    let data = vec![1.0, 2.0, 3.0, 4.0];
    CsMat::new_csc((rows, cols), indptr, indices, data)
}

#[test]
fn aggregate_sample_mode() {
    let matrix = sample_matrix();
    let aggregated = aggregate(&matrix, AggregationMode::Sample, None);
    assert_eq!(aggregated.genes, 2);
    assert_eq!(aggregated.samples, 1);
    let values = aggregated.values;
    assert!((values[0] - (1.0 + 3.0) / 3.0).abs() < 1e-6);
    assert!((values[1] - (2.0 + 4.0) / 3.0).abs() < 1e-6);
}

#[test]
fn aggregate_cell_mode() {
    let matrix = sample_matrix();
    let aggregated = aggregate(&matrix, AggregationMode::Cell, None);
    assert_eq!(aggregated.genes, 2);
    assert_eq!(aggregated.samples, 3);
    assert_eq!(aggregated.values, vec![1.0, 3.0, 0.0, 2.0, 0.0, 4.0]);
}

#[test]
fn aggregate_cluster_mode() {
    let matrix = sample_matrix();
    let clusters = ClusterMap {
        cluster_ids: vec!["A".to_string(), "B".to_string()],
        cell_to_cluster: vec![0, 1, 0],
    };
    let aggregated = aggregate(&matrix, AggregationMode::Cluster, Some(&clusters));
    assert_eq!(aggregated.genes, 2);
    assert_eq!(aggregated.samples, 2);
    assert_eq!(aggregated.values, vec![0.5, 3.0, 3.0, 0.0]);
}

#[test]
fn soa_preparation_fills_missing_genes_with_zero() {
    let features = vec!["MT1", "N1", "R1"]
        .into_iter()
        .map(|s| s.to_string())
        .collect();
    let barcodes = vec!["C1", "C2"]
        .into_iter()
        .map(|s| s.to_string())
        .collect();

    let rows = 3;
    let cols = 2;
    let indptr = vec![0, 3, 6];
    let indices = vec![0, 1, 2, 0, 1, 2];
    let data = vec![1.0, 3.0, 5.0, 2.0, 4.0, 6.0];
    let matrix = CsMat::new_csc((rows, cols), indptr, indices, data);

    let mtx = MtxInput {
        matrix,
        features,
        barcodes,
    };

    let geneset = GeneSet {
        mtdna_complex_i: vec!["MT1".to_string()],
        mtdna_complex_iii: vec!["MT3".to_string()],
        mtdna_complex_iv: vec!["MT4".to_string()],
        mtdna_complex_v: vec!["MT5".to_string()],
        nuclear_oxphos_complex_i: vec!["N1".to_string()],
        nuclear_oxphos_complex_ii: vec!["N2".to_string()],
        nuclear_oxphos_complex_iii: vec!["N3".to_string()],
        nuclear_oxphos_complex_iv: vec!["N4".to_string()],
        nuclear_oxphos_complex_v: vec!["N5".to_string()],
        ros_detox_genes: vec!["R1".to_string(), "R2".to_string()],
        mitophagy_genes: vec!["M1".to_string()],
        dynamics_fusion: vec!["F1".to_string()],
        dynamics_fission: vec!["F2".to_string()],
        biogenesis_genes: vec!["B1".to_string()],
    };

    let gene_index = GeneIndex::from_feature_list(&mtx.features);
    let resolved = resolve_all_genesets(&gene_index, &geneset);

    let prepared = prepare_expression_with_clusters(&mtx, &resolved, AggregationMode::Cell, None)
        .expect("prepare expression");

    assert_eq!(prepared.soa.samples, 2);
    assert_eq!(prepared.soa.get(0, 0), 1.0);
    assert_eq!(prepared.soa.get(0, 1), 2.0);
    assert_eq!(prepared.soa.get(1, 0), 0.0);
    assert_eq!(prepared.soa.get(1, 1), 0.0);
}
