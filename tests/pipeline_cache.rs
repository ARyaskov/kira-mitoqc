use std::fs;
use std::path::PathBuf;

use kira_mitoqc::cache::{mmap_organelle_bin, write_organelle_bin_from_mtx};
use kira_mitoqc::io::mtx::{discover_dataset_files, resolve_shared_cache_filename};

fn temp_dir(name: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!("kira_mitoqc_pipeline_{name}_{nanos}"));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn write_basic_dataset(dir: &PathBuf, prefix: Option<&str>) {
    let stem = |name: &str| match prefix {
        Some(p) => format!("{p}_{name}"),
        None => name.to_string(),
    };

    fs::write(
        dir.join(stem("matrix.mtx")),
        "%%MatrixMarket matrix coordinate integer general\n3 2 3\n1 1 5\n3 1 1\n2 2 7\n",
    )
    .unwrap();
    fs::write(
        dir.join(stem("features.tsv")),
        "ENSG1\tMT-ND1\tGene Expression\nENSG2\tATP5F1A\tGene Expression\nENSG3\tSOD2\tGene Expression\n",
    )
    .unwrap();
    fs::write(dir.join(stem("barcodes.tsv")), "C1\nC2\n").unwrap();
}

#[test]
fn prefix_detection_non_prefixed() {
    let dir = temp_dir("noprefix");
    write_basic_dataset(&dir, None);
    let d = discover_dataset_files(&dir).expect("discover");
    assert_eq!(d.prefix, None);
}

#[test]
fn prefix_detection_prefixed() {
    let dir = temp_dir("prefixed");
    write_basic_dataset(&dir, Some("GSM123"));
    let d = discover_dataset_files(&dir).expect("discover");
    assert_eq!(d.prefix.as_deref(), Some("GSM123"));
}

#[test]
fn shared_bin_filename_resolution() {
    assert_eq!(resolve_shared_cache_filename(None), "kira-organelle.bin");
    assert_eq!(
        resolve_shared_cache_filename(Some("GSM123")),
        "GSM123.kira-organelle.bin"
    );
}

#[test]
fn organelle_cache_roundtrip_and_csc_equivalence() {
    let dir = temp_dir("roundtrip");
    write_basic_dataset(&dir, None);
    let out = dir.join("kira-organelle.bin");

    write_organelle_bin_from_mtx(&out, &dir, Some(1)).expect("write");
    let view = mmap_organelle_bin(&out).expect("mmap");

    assert_eq!(view.n_genes, 3);
    assert_eq!(view.n_cells, 2);
    assert_eq!(view.nnz, 3);
    assert_eq!(view.genes, vec!["MT-ND1", "ATP5F1A", "SOD2"]);
    assert_eq!(view.barcodes, vec!["C1", "C2"]);
    assert_eq!(view.col_ptr(), &[0, 2, 3]);
    assert_eq!(view.row_idx(), &[0, 2, 1]);
    assert_eq!(view.values_u32(), &[5, 1, 7]);
}

#[test]
fn header_and_crc_are_valid() {
    let dir = temp_dir("crc");
    write_basic_dataset(&dir, None);
    let out = dir.join("kira-organelle.bin");

    write_organelle_bin_from_mtx(&out, &dir, Some(1)).expect("write");
    let bytes = fs::read(&out).expect("read bin");
    assert!(bytes.len() >= 256);
    let header = &bytes[..256];
    assert_eq!(&header[0..4], b"KORG");
    assert_eq!(
        u16::from_le_bytes(*header[4..6].as_array::<2>().unwrap()),
        1
    );
    assert_eq!(
        u16::from_le_bytes(*header[6..8].as_array::<2>().unwrap()),
        0
    );
    assert_eq!(
        u32::from_le_bytes(*header[8..12].as_array::<4>().unwrap()),
        0x1234_5678
    );
    assert_eq!(
        u32::from_le_bytes(*header[12..16].as_array::<4>().unwrap()),
        256
    );
    let file_bytes = u64::from_le_bytes(*header[112..120].as_array::<8>().unwrap()) as usize;
    assert_eq!(file_bytes, bytes.len());
    let stored_crc = u64::from_le_bytes(*header[120..128].as_array::<8>().unwrap());
    let mut crc_header = [0u8; 256];
    crc_header.copy_from_slice(header);
    crc_header[120..128].fill(0);
    let computed_crc = crc64_ecma(&crc_header);
    assert_eq!(stored_crc, computed_crc);
}

#[test]
fn pipeline_like_path_uses_written_cache() {
    let dir = temp_dir("pipeline_path");
    write_basic_dataset(&dir, Some("XYZ"));
    let discovery = discover_dataset_files(&dir).expect("discover");
    let shared = discovery
        .input_dir
        .join(resolve_shared_cache_filename(discovery.prefix.as_deref()));

    write_organelle_bin_from_mtx(&shared, &dir, Some(1)).expect("write");
    let view = mmap_organelle_bin(&shared).expect("mmap");
    let input = view.to_mtx_input();

    assert_eq!(input.features.len(), 3);
    assert_eq!(input.barcodes.len(), 2);
    assert_eq!(input.matrix.nnz(), 3);
}

fn crc64_ecma(bytes: &[u8]) -> u64 {
    let poly: u64 = 0x42F0_E1EB_A9EA_3693;
    let mut crc: u64 = 0;
    for byte in bytes {
        crc ^= (*byte as u64) << 56;
        for _ in 0..8 {
            if (crc & 0x8000_0000_0000_0000) != 0 {
                crc = (crc << 1) ^ poly;
            } else {
                crc <<= 1;
            }
        }
    }
    crc
}
