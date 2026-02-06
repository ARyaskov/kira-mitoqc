mod gzip_tests {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use kira_mitoqc::io::mtx::{load_mtx_dir, load_mtx_metadata};
    use std::fs::{self, File};
    use std::io::Write;

    fn temp_dir(name: &str) -> std::path::PathBuf {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!("kira_mitoqc_gzip_{name}_{nanos}"));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_gz(path: &std::path::Path, contents: &str) {
        let file = File::create(path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(contents.as_bytes()).unwrap();
        encoder.finish().unwrap();
    }

    #[test]
    fn features_gz_only() {
        let dir = temp_dir("features");
        write_gz(&dir.join("features.tsv.gz"), "G1\tX\nG2\tY\n");
        fs::write(dir.join("barcodes.tsv"), "C1\nC2\n").unwrap();
        fs::write(
            dir.join("matrix.mtx"),
            "%%MatrixMarket matrix coordinate real general\n2 2 2\n1 1 1.0\n2 2 2.0\n",
        )
        .unwrap();

        let (features, _) = load_mtx_metadata(&dir, Some(0)).expect("metadata");
        assert_eq!(features, vec!["G1".to_string(), "G2".to_string()]);
    }

    #[test]
    fn barcodes_gz_only() {
        let dir = temp_dir("barcodes");
        fs::write(dir.join("features.tsv"), "G1\tX\nG2\tY\n").unwrap();
        write_gz(&dir.join("barcodes.tsv.gz"), "C1\nC2\n");
        fs::write(
            dir.join("matrix.mtx"),
            "%%MatrixMarket matrix coordinate real general\n2 2 2\n1 1 1.0\n2 2 2.0\n",
        )
        .unwrap();

        let (_, barcodes) = load_mtx_metadata(&dir, Some(0)).expect("metadata");
        assert_eq!(barcodes, vec!["C1".to_string(), "C2".to_string()]);
    }

    #[test]
    fn mixed_plain_and_gz() {
        let dir = temp_dir("mixed");
        write_gz(&dir.join("features.tsv.gz"), "G1\tX\nG2\tY\n");
        fs::write(dir.join("barcodes.tsv"), "C1\nC2\n").unwrap();
        write_gz(
            &dir.join("matrix.mtx.gz"),
            "%%MatrixMarket matrix coordinate real general\n2 2 2\n1 1 1.0\n2 2 2.0\n",
        );

        let input = load_mtx_dir(&dir, Some(0)).expect("mtx");
        assert_eq!(input.features.len(), 2);
        assert_eq!(input.barcodes.len(), 2);
        assert_eq!(input.matrix.rows(), 2);
        assert_eq!(input.matrix.cols(), 2);
    }
}
