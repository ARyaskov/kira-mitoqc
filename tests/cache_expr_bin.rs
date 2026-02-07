use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use kira_mitoqc::cache::{
    CacheError, ExprCacheMode, mmap_expr_bin, write_expr_bin, write_expr_bin_with_mode,
};
use kira_mitoqc::data::ExpressionSoA;

fn temp_path(name: &str) -> PathBuf {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("kira_mitoqc_cache_{name}_{nanos}"))
}

#[test]
fn roundtrip_write_and_mmap() {
    let soa = ExpressionSoA {
        values: vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        genes: 2,
        samples: 3,
    };
    let path = temp_path("roundtrip");
    write_expr_bin(&path, &soa).expect("write cache");

    let view = mmap_expr_bin(&path).expect("mmap cache");
    assert_eq!(view.genes, soa.genes);
    assert_eq!(view.samples, soa.samples);
    assert_eq!(view.values, soa.values);
    assert_eq!(view.get(1, 2), 6.0);
}

#[test]
fn invalid_magic_is_rejected() {
    let path = temp_path("bad_magic");
    let mut file = File::create(&path).unwrap();
    file.write_all(b"BADMAGIC").unwrap();
    file.write_all(&1u32.to_le_bytes()).unwrap();
    file.write_all(&1u32.to_le_bytes()).unwrap();
    file.write_all(&1u32.to_le_bytes()).unwrap();
    file.write_all(&0u32.to_le_bytes()).unwrap();
    file.write_all(&[0u8; 4]).unwrap();

    let err = mmap_expr_bin(&path).unwrap_err();
    match err {
        CacheError::InvalidMagic { .. } => {}
        _ => panic!("unexpected error: {err:?}"),
    }
}

#[test]
fn size_mismatch_is_detected() {
    let path = temp_path("size_mismatch");
    let mut file = File::create(&path).unwrap();
    file.write_all(b"KIRAMTX\0").unwrap();
    file.write_all(&1u32.to_le_bytes()).unwrap();
    file.write_all(&2u32.to_le_bytes()).unwrap();
    file.write_all(&3u32.to_le_bytes()).unwrap();
    file.write_all(&0u32.to_le_bytes()).unwrap();
    file.write_all(&[0u8; 4]).unwrap();

    let err = mmap_expr_bin(&path).unwrap_err();
    match err {
        CacheError::SizeMismatch { .. } => {}
        _ => panic!("unexpected error: {err:?}"),
    }
}

#[test]
fn mode_metadata_roundtrip() {
    let soa = ExpressionSoA {
        values: vec![1.0, 2.0, 3.0, 4.0],
        genes: 2,
        samples: 2,
    };
    let path = temp_path("mode_roundtrip");
    write_expr_bin_with_mode(&path, &soa, ExprCacheMode::Cell).expect("write cache");

    let view = mmap_expr_bin(&path).expect("mmap cache");
    assert_eq!(view.mode, ExprCacheMode::Cell);
}

#[test]
fn legacy_writer_keeps_unknown_mode() {
    let soa = ExpressionSoA {
        values: vec![1.0, 2.0, 3.0, 4.0],
        genes: 2,
        samples: 2,
    };
    let path = temp_path("unknown_mode");
    write_expr_bin(&path, &soa).expect("write cache");

    let view = mmap_expr_bin(&path).expect("mmap cache");
    assert_eq!(view.mode, ExprCacheMode::Unknown);
}
