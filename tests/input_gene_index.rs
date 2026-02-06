use kira_mitoqc::core::types::GeneSet;
use kira_mitoqc::input::{GeneIndex, resolve_all_genesets};

fn vecs(values: &[&str]) -> Vec<String> {
    values.iter().map(|v| v.to_string()).collect()
}

#[test]
fn gene_index_preserves_first_occurrence_order() {
    let features = vecs(&["A", "B", "A", "C"]);
    let index = GeneIndex::from_feature_list(&features);
    let resolution = index.resolve(&vecs(&["A", "B", "C"]));

    assert_eq!(
        resolution.genes,
        vec!["A".to_string(), "B".to_string(), "C".to_string()]
    );
    assert_eq!(resolution.found, vec![0, 1, 3]);
    assert!(resolution.missing.is_empty());
}

#[test]
fn gene_index_resolve_reports_missing() {
    let features = vecs(&["G1", "G2"]);
    let index = GeneIndex::from_feature_list(&features);
    let resolution = index.resolve(&vecs(&["G2", "G3", "G1"]));

    assert_eq!(
        resolution.genes,
        vec!["G2".to_string(), "G3".to_string(), "G1".to_string()]
    );
    assert_eq!(resolution.found, vec![1, 0]);
    assert_eq!(resolution.missing, vec!["G3".to_string()]);
}

#[test]
fn resolve_all_genesets_resolves_every_section() {
    let geneset = GeneSet {
        mtdna_complex_i: vecs(&["MT1"]),
        mtdna_complex_iii: vecs(&["MT3"]),
        mtdna_complex_iv: vecs(&["MT4"]),
        mtdna_complex_v: vecs(&["MT5"]),
        nuclear_oxphos_complex_i: vecs(&["N1"]),
        nuclear_oxphos_complex_ii: vecs(&["N2"]),
        nuclear_oxphos_complex_iii: vecs(&["N3"]),
        nuclear_oxphos_complex_iv: vecs(&["N4"]),
        nuclear_oxphos_complex_v: vecs(&["N5"]),
        ros_detox_genes: vecs(&["R1"]),
        mitophagy_genes: vecs(&["M1"]),
        dynamics_fusion: vecs(&["F1"]),
        dynamics_fission: vecs(&["F2"]),
        biogenesis_genes: vecs(&["B1"]),
    };

    let features = vecs(&[
        "MT1", "MT3", "MT4", "MT5", "N1", "N2", "N3", "N4", "N5", "R1", "M1", "F1", "F2", "B1",
    ]);
    let index = GeneIndex::from_feature_list(&features);

    let resolved = resolve_all_genesets(&index, &geneset);

    assert_eq!(resolved.mtdna_complex_i.found, vec![0]);
    assert_eq!(resolved.mtdna_complex_iii.found, vec![1]);
    assert_eq!(resolved.mtdna_complex_iv.found, vec![2]);
    assert_eq!(resolved.mtdna_complex_v.found, vec![3]);
    assert_eq!(resolved.mtdna_all.found, vec![0, 1, 2, 3]);
    assert_eq!(resolved.nuclear_oxphos_complex_i.found, vec![4]);
    assert_eq!(resolved.nuclear_oxphos_complex_ii.found, vec![5]);
    assert_eq!(resolved.nuclear_oxphos_complex_iii.found, vec![6]);
    assert_eq!(resolved.nuclear_oxphos_complex_iv.found, vec![7]);
    assert_eq!(resolved.nuclear_oxphos_complex_v.found, vec![8]);
    assert_eq!(resolved.nuclear_oxphos_all.found, vec![4, 5, 6, 7, 8]);
    assert_eq!(resolved.ros.found, vec![9]);
    assert_eq!(resolved.mitophagy.found, vec![10]);
    assert_eq!(resolved.fusion.found, vec![11]);
    assert_eq!(resolved.fission.found, vec![12]);
    assert_eq!(resolved.biogenesis.found, vec![13]);
}
