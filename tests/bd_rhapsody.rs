use std::fs;
use std::path::PathBuf;

use kira_mitoqc::core::types::GeneSet;
use kira_mitoqc::data::aggregate::AggregationMode;
use kira_mitoqc::data::prepare_expression_with_clusters;
use kira_mitoqc::input::bd_rhapsody::{
    compute_mito_fraction_from_file, load_bd_rhapsody, load_bd_rhapsody_metadata,
    resolve_bd_input_path,
};
use kira_mitoqc::input::{GeneIndex, resolve_all_genesets};
use kira_mitoqc::io::mtx::MtxInput;

fn temp_file(name: &str, content: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!("kira_mitoqc_{name}_{nanos}.tsv"));
    fs::write(&path, content).expect("write temp file");
    path
}

fn tiny_geneset() -> GeneSet {
    GeneSet {
        mtdna_complex_i: vec!["MT-ND1".to_string()],
        mtdna_complex_iii: vec!["MT-CYB".to_string()],
        mtdna_complex_iv: vec!["MT-CO1".to_string()],
        mtdna_complex_v: vec!["MT-ATP6".to_string(), "MT-ATP8".to_string()],
        nuclear_oxphos_complex_i: vec!["NDUFS1".to_string()],
        nuclear_oxphos_complex_ii: vec!["SDHA".to_string()],
        nuclear_oxphos_complex_iii: vec!["UQCRC1".to_string()],
        nuclear_oxphos_complex_iv: vec!["COX4I1".to_string()],
        nuclear_oxphos_complex_v: vec!["ATP5F1A".to_string()],
        ros_detox_genes: vec!["SOD2".to_string()],
        mitophagy_genes: vec!["PINK1".to_string()],
        dynamics_fusion: vec!["MFN1".to_string()],
        dynamics_fission: vec!["DNM1L".to_string()],
        biogenesis_genes: vec!["TFAM".to_string()],
    }
}

#[test]
fn loads_minimal_bd_rhapsody_matrix() {
    let path = temp_file(
        "bd_minimal",
        "#Sequencing: BD Rhapsody Whole Transcriptome 3-end sequencing\n\
         cellA\tcellB\tcellC\n\
         MT-ND1\t1.0\t2.0\t3.0\n\
         NDUFS1\t4.0\t5.0\t6.0\n\
         SOD2\t7.0\t8.0\t9.0\n",
    );
    let input = load_bd_rhapsody(&path).expect("load bd");
    assert_eq!(input.features.len(), 3);
    assert_eq!(input.barcodes.len(), 3);
    assert_eq!(input.matrix.rows(), 3);
    assert_eq!(input.matrix.cols(), 3);
    assert_eq!(input.matrix.get(0, 0).copied().unwrap_or(0.0), 1.0);
    assert_eq!(input.matrix.get(1, 2).copied().unwrap_or(0.0), 6.0);
    assert_eq!(input.matrix.get(2, 1).copied().unwrap_or(0.0), 8.0);

    let (features, barcodes) = load_bd_rhapsody_metadata(&path).expect("meta");
    assert_eq!(features, input.features);
    assert_eq!(barcodes, input.barcodes);
}

#[test]
fn duplicate_gene_symbols_are_loaded_deterministically() {
    let path = temp_file(
        "bd_dups",
        "#meta\n\
         cellA\tcellB\n\
         MT-ND1\t1.0\t2.0\n\
         MT-ND1\t3.0\t4.0\n\
         NDUFS1\t5.0\t6.0\n",
    );
    let input = load_bd_rhapsody(&path).expect("load bd");
    assert_eq!(input.features, vec!["MT-ND1", "MT-ND1", "NDUFS1"]);
    assert_eq!(input.matrix.rows(), 3);
    assert_eq!(input.matrix.cols(), 2);
}

#[test]
fn mixed_integer_and_float_values_parse() {
    let path = temp_file(
        "bd_mixed",
        "#meta\n\
         cellA\tcellB\n\
         MT-ND1\t1\t2.5\n\
         NDUFS1\t3.0\t4\n",
    );
    let input = load_bd_rhapsody(&path).expect("load bd");
    assert_eq!(input.matrix.get(0, 0).copied().unwrap_or(0.0), 1.0);
    assert_eq!(input.matrix.get(0, 1).copied().unwrap_or(0.0), 2.5);
    assert_eq!(input.matrix.get(1, 0).copied().unwrap_or(0.0), 3.0);
    assert_eq!(input.matrix.get(1, 1).copied().unwrap_or(0.0), 4.0);
}

#[test]
fn prepared_expression_is_stable_under_gene_row_shuffle() {
    let path_a = temp_file(
        "bd_order_a",
        "#meta\n\
         cellA\tcellB\n\
         MT-ND1\t1.0\t2.0\n\
         NDUFS1\t3.0\t4.0\n\
         SOD2\t5.0\t6.0\n",
    );
    let path_b = temp_file(
        "bd_order_b",
        "#meta\n\
         cellA\tcellB\n\
         SOD2\t5.0\t6.0\n\
         NDUFS1\t3.0\t4.0\n\
         MT-ND1\t1.0\t2.0\n",
    );
    let a = load_bd_rhapsody(&path_a).expect("load a");
    let b = load_bd_rhapsody(&path_b).expect("load b");
    let geneset = tiny_geneset();

    let idx_a = GeneIndex::from_feature_list(&a.features);
    let idx_b = GeneIndex::from_feature_list(&b.features);
    let resolved_a = resolve_all_genesets(&idx_a, &geneset);
    let resolved_b = resolve_all_genesets(&idx_b, &geneset);

    let prep_a = prepare_expression_with_clusters(
        &MtxInput {
            matrix: a.matrix,
            features: a.features,
            barcodes: a.barcodes,
        },
        &resolved_a,
        AggregationMode::Cell,
        None,
    )
    .expect("prepare a");
    let prep_b = prepare_expression_with_clusters(
        &MtxInput {
            matrix: b.matrix,
            features: b.features,
            barcodes: b.barcodes,
        },
        &resolved_b,
        AggregationMode::Cell,
        None,
    )
    .expect("prepare b");

    assert_eq!(prep_a.soa.values, prep_b.soa.values);
}

#[test]
fn mito_fraction_uses_sum_mito_over_sum_all() {
    let path = temp_file(
        "bd_fraction",
        "#meta\n\
         cellA\tcellB\tcellC\n\
         MT-ND1\t2.0\t0.0\t1.0\n\
         NDUFS1\t2.0\t4.0\t1.0\n\
         SOD2\t6.0\t6.0\t8.0\n",
    );
    let mito = std::collections::BTreeSet::from(["MT-ND1".to_string()]);
    let fraction =
        compute_mito_fraction_from_file(&path, &mito, AggregationMode::Cell, None).expect("frac");
    assert_eq!(fraction.len(), 3);
    assert!((fraction[0] - 0.2).abs() < 1e-6);
    assert!((fraction[1] - 0.0).abs() < 1e-6);
    assert!((fraction[2] - 0.1).abs() < 1e-6);
}

#[test]
fn supports_gene_header_column_in_raw_counts_format() {
    let path = temp_file(
        "bd_gene_header",
        "gene\tcellA\tcellB\n\
         MT-ND1\t1\t2\n\
         NDUFS1\t3\t4\n",
    );
    let input = load_bd_rhapsody(&path).expect("load");
    assert_eq!(input.features, vec!["MT-ND1", "NDUFS1"]);
    assert_eq!(input.barcodes, vec!["cellA", "cellB"]);
    assert_eq!(input.matrix.get(0, 0).copied().unwrap_or(0.0), 1.0);
    assert_eq!(input.matrix.get(1, 1).copied().unwrap_or(0.0), 4.0);
}

#[test]
fn resolve_prefixed_raw_counts_in_directory() {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_mitoqc_bd_resolve_{nanos}"));
    fs::create_dir_all(&dir).expect("mkdir");
    let file = dir.join("GSM3635278_human_p1t1_raw_counts.tsv");
    fs::write(&file, "gene\tcellA\tcellB\nGENE1\t1\t2\n").expect("write");

    let resolved = resolve_bd_input_path(&dir).expect("resolve");
    assert_eq!(resolved, file);
}

#[test]
fn supports_empty_first_header_cell_in_raw_counts() {
    let path = temp_file(
        "bd_empty_first_header",
        "\tcellA\tcellB\n\
         MT-ND1\t1\t2\n\
         NDUFS1\t3\t4\n",
    );
    let input = load_bd_rhapsody(&path).expect("load");
    assert_eq!(input.barcodes, vec!["cellA", "cellB"]);
    assert_eq!(input.features, vec!["MT-ND1", "NDUFS1"]);
}

#[test]
fn keeps_trailing_empty_columns_in_rows() {
    let path = temp_file(
        "bd_trailing_empty",
        "gene\tcellA\tcellB\tcellC\n\
         MT-ND1\t1\t2\t\n\
         NDUFS1\t3\t4\t5\n",
    );
    let input = load_bd_rhapsody(&path).expect("load");
    assert_eq!(input.barcodes, vec!["cellA", "cellB", "cellC"]);
    assert_eq!(input.matrix.cols(), 3);
    assert_eq!(input.matrix.get(0, 2).copied().unwrap_or(0.0), 0.0);
}

#[test]
fn supports_cell_major_raw_counts_with_barcode_header() {
    let path = temp_file(
        "bd_cell_major",
        "barcode\tMT-ND1\tNDUFS1\n\
         CELL_A\t1\t3\n\
         CELL_B\t2\t4\n",
    );
    let input = load_bd_rhapsody(&path).expect("load");
    assert_eq!(input.features, vec!["MT-ND1", "NDUFS1"]);
    assert_eq!(input.barcodes, vec!["CELL_A", "CELL_B"]);
    assert_eq!(input.matrix.get(0, 0).copied().unwrap_or(0.0), 1.0);
    assert_eq!(input.matrix.get(1, 1).copied().unwrap_or(0.0), 4.0);
}
