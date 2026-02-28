use kira_mitoqc::cache::{mmap_expr_bin, write_expr_bin};
use kira_mitoqc::compute::dispatch::compute_primitives;
use kira_mitoqc::compute::{GeneOffsets, SIMD_EPS, scalar};
use kira_mitoqc::core::types::GeneSet;
use kira_mitoqc::data::ExpressionSoA;
use kira_mitoqc::input::{GeneIndex, resolve_all_genesets};

fn vecs(values: &[&str]) -> Vec<String> {
    values.iter().map(|v| v.to_string()).collect()
}

fn make_fixture(samples: usize) -> (ExpressionSoA, GeneSet, Vec<String>) {
    let genes = vec![
        "MT-ND1", "MT-CYB", "MT-CO1", "MT-ATP6", "MT-ATP8", "NDUFS1", "SDHA", "UQCRC1", "COX4I1",
        "ATP5F1A", "SOD2", "PINK1", "MFN1", "DNM1L", "TFAM",
    ];
    let mut values = Vec::new();
    for g in 0..genes.len() {
        for s in 0..samples {
            values.push((g as f32) + (s as f32) * 0.1);
        }
    }
    let soa = ExpressionSoA {
        values,
        genes: genes.len(),
        samples,
    };

    let geneset = GeneSet {
        mtdna_complex_i: vecs(&["MT-ND1"]),
        mtdna_complex_iii: vecs(&["MT-CYB"]),
        mtdna_complex_iv: vecs(&["MT-CO1"]),
        mtdna_complex_v: vecs(&["MT-ATP6", "MT-ATP8"]),
        nuclear_oxphos_complex_i: vecs(&["NDUFS1"]),
        nuclear_oxphos_complex_ii: vecs(&["SDHA"]),
        nuclear_oxphos_complex_iii: vecs(&["UQCRC1"]),
        nuclear_oxphos_complex_iv: vecs(&["COX4I1"]),
        nuclear_oxphos_complex_v: vecs(&["ATP5F1A"]),
        ros_detox_genes: vecs(&["SOD2"]),
        mitophagy_genes: vecs(&["PINK1"]),
        dynamics_fusion: vecs(&["MFN1"]),
        dynamics_fission: vecs(&["DNM1L"]),
        biogenesis_genes: vecs(&["TFAM"]),
    };

    (
        soa,
        geneset,
        genes.into_iter().map(|s| s.to_string()).collect(),
    )
}

fn assert_close(a: &[f32], b: &[f32]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b.iter()) {
        let diff = (x - y).abs();
        assert!(diff <= SIMD_EPS, "diff {diff} exceeded {SIMD_EPS}");
    }
}

#[test]
fn compute_primitives_matches_scalar() {
    let (soa, geneset, features) = make_fixture(7);
    let gene_index = GeneIndex::from_feature_list(&features);
    let resolved = resolve_all_genesets(&gene_index, &geneset);

    let path = std::env::temp_dir().join("kira_mitoqc_primitives.bin");
    write_expr_bin(&path, &soa).expect("write cache");
    let view = mmap_expr_bin(&path).expect("mmap cache");

    let offsets = GeneOffsets::from_resolved(&resolved);
    let scalar_out = scalar::compute_primitives_scalar(&view, &offsets);
    let dispatch_out = compute_primitives(&view, &resolved);

    assert_close(&scalar_out.mtdna_mean, &dispatch_out.mtdna_mean);
    assert_close(&scalar_out.nuclear_mean, &dispatch_out.nuclear_mean);
    assert_close(&scalar_out.c_i, &dispatch_out.c_i);
    assert_close(&scalar_out.c_iii, &dispatch_out.c_iii);
    assert_close(&scalar_out.c_iv, &dispatch_out.c_iv);
    assert_close(&scalar_out.c_v, &dispatch_out.c_v);
    assert_close(&scalar_out.ros_mean, &dispatch_out.ros_mean);
    assert_close(&scalar_out.mitophagy_mean, &dispatch_out.mitophagy_mean);
    assert_close(&scalar_out.fusion_mean, &dispatch_out.fusion_mean);
    assert_close(&scalar_out.fission_mean, &dispatch_out.fission_mean);
    assert_close(&scalar_out.biogenesis_mean, &dispatch_out.biogenesis_mean);
    assert_close(&scalar_out.atp_mt, &dispatch_out.atp_mt);
    assert_close(&scalar_out.atp_nu, &dispatch_out.atp_nu);
    assert_close(&scalar_out.stoich_variance, &dispatch_out.stoich_variance);
}
