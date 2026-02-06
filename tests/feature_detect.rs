use kira_mitoqc::io::feature_detect::detect_gene_symbol_column;
use kira_mitoqc::io::mtx::load_mtx_metadata;
use std::io::Cursor;

#[test]
fn detect_10x_v3_symbols_col2() {
    let data = "ENSG0001\tMT-ND1\tGene Expression\nENSG0002\tCOX4I1\tGene Expression\n";
    let mut cursor = Cursor::new(data.as_bytes());
    let col = detect_gene_symbol_column(&mut cursor, 100).expect("detect");
    assert_eq!(col, 1);
}

#[test]
fn detect_10x_v2_symbols_col1() {
    let data = "MT-ND1\tENSG0001\nCOX4I1\tENSG0002\n";
    let mut cursor = Cursor::new(data.as_bytes());
    let col = detect_gene_symbol_column(&mut cursor, 100).expect("detect");
    assert_eq!(col, 0);
}

#[test]
fn detect_symbols_only() {
    let data = "MT-CO1\nMT-ND2\nATP5F1A\n";
    let mut cursor = Cursor::new(data.as_bytes());
    let col = detect_gene_symbol_column(&mut cursor, 100).expect("detect");
    assert_eq!(col, 0);
}

#[test]
fn detect_ensembl_only_fails() {
    let data = "ENSG000001\nENSG000002\n";
    let mut cursor = Cursor::new(data.as_bytes());
    let err = detect_gene_symbol_column(&mut cursor, 100).unwrap_err();
    match err {
        kira_mitoqc::input::InputError::GeneSymbolNotDetected => {}
        _ => panic!("unexpected error: {err:?}"),
    }
}

#[test]
fn manual_override_precedence() {
    let dir = std::env::temp_dir().join("kira_mitoqc_detect_override");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(
        dir.join("features.tsv"),
        "ENSG0001\tMT-ND1\nENSG0002\tCOX4I1\n",
    )
    .unwrap();
    std::fs::write(dir.join("barcodes.tsv"), "C1\nC2\n").unwrap();

    let (features, _) = load_mtx_metadata(&dir, Some(0)).expect("meta");
    assert_eq!(features[0], "ENSG0001");
}
