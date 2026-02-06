#[cfg(feature = "h5ad")]
mod h5ad_tests {
    use hdf5::File;
    use kira_mitoqc::io::h5ad::{load_h5ad, load_h5ad_clusters, load_h5ad_metadata};
    use std::fs;

    fn temp_path(name: &str) -> std::path::PathBuf {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("kira_mitoqc_h5ad_{name}_{nanos}.h5ad"))
    }

    fn write_basic_h5ad(path: &std::path::Path) {
        let file = File::create(path).unwrap();

        let obs = file.create_group("obs").unwrap();
        obs.new_dataset_builder()
            .with_data(&["cell1", "cell2"])
            .create("_index")
            .unwrap();
        obs.new_dataset_builder()
            .with_data(&["A", "B"])
            .create("cluster")
            .unwrap();

        let var = file.create_group("var").unwrap();
        var.new_dataset_builder()
            .with_data(&["gene1", "gene2"])
            .create("_index")
            .unwrap();
        var.new_dataset_builder()
            .with_data(&["gA", "gB"])
            .create("gene_symbols")
            .unwrap();

        let x = file
            .new_dataset_builder()
            .with_data(&[0.1_f32, 0.2_f32, 0.3_f32, 0.4_f32])
            .shape([2, 2])
            .create("X")
            .unwrap();
        let _ = x;
    }

    #[test]
    fn load_dense_h5ad() {
        let path = temp_path("dense");
        write_basic_h5ad(&path);

        let input = load_h5ad(&path, None).expect("load h5ad");
        assert_eq!(input.features, vec!["gA".to_string(), "gB".to_string()]);
        assert_eq!(
            input.barcodes,
            vec!["cell1".to_string(), "cell2".to_string()]
        );
        assert_eq!(input.matrix.rows(), 2);
        assert_eq!(input.matrix.cols(), 2);

        fs::remove_file(path).ok();
    }

    #[test]
    fn gene_symbol_key_precedence() {
        let path = temp_path("symbols");
        write_basic_h5ad(&path);

        let (features, _) = load_h5ad_metadata(&path, Some("_index")).expect("meta");
        assert_eq!(features, vec!["gene1".to_string(), "gene2".to_string()]);

        fs::remove_file(path).ok();
    }

    #[test]
    fn cluster_labels_from_obs() {
        let path = temp_path("clusters");
        write_basic_h5ad(&path);

        let cluster_map = load_h5ad_clusters(&path, "cluster").expect("clusters");
        assert_eq!(
            cluster_map.cluster_ids,
            vec!["A".to_string(), "B".to_string()]
        );
        assert_eq!(cluster_map.cell_to_cluster, vec![0, 1]);

        fs::remove_file(path).ok();
    }
}
