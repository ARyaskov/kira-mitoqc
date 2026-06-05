use kira_mitoqc::cache::{ExprCacheMode, mmap_expr_bin, write_expr_bin_with_mode};
use kira_mitoqc::classify::classify_v1_with_redox;
use kira_mitoqc::compute::{GeneOffsets, scalar};
use kira_mitoqc::config::refs::{Eps, Metadata, Normalization, Qc, Refs, RefsV1, Thresholds};
use kira_mitoqc::core::types::{GeneSet, ProxyKey};
use kira_mitoqc::data::{ExpressionSoA, SoaIndex};
use kira_mitoqc::input::{GeneIndex, resolve_all_genesets};
use kira_mitoqc::proxy::compute_proxies_v1;
use kira_mitoqc::redox::{RedoxRegime, compute_redox_metrics};
use kira_mitoqc::score::{DecayScoreVec, compute_axes_v1};

fn vecs(values: &[&str]) -> Vec<String> {
    values.iter().map(|v| v.to_string()).collect()
}

fn refs_fixture() -> RefsV1 {
    RefsV1 {
        metadata: Metadata {
            version: "v1.0".to_string(),
            description: "test".to_string(),
            notes: "test".to_string(),
        },
        eps: Eps { value: 1e-6 },
        refs: Refs {
            stoich_ref: 0.25,
            uncoupling_ref: 1.0,
            ros_ref: 1.0,
            redox_ref: 1.0,
            atp_ref: 1.0,
            mito_ref: 1.0,
            dyn_ref: 1.0,
            bio_ref: 1.0,
        },
        thresholds: Thresholds {
            ros_high: 0.75,
            bioenergetics_low: 0.5,
            bioenergetics_high: 0.8,
            dynamics_high: 0.7,
            regulation_low: 0.4,
            structural_bio_min: 0.6,
            structural_dyn_min: 0.6,
        },
        normalization: Normalization {
            log1p: true,
            expression_unit: "log1p".to_string(),
        },
        qc: Qc {
            min_mtdna_genes_found: 1,
            min_nuclear_oxphos_found: 1,
            min_ros_genes_found: 1,
        },
    }
}

fn weights_fixture() -> kira_mitoqc::config::weights::WeightsV1 {
    kira_mitoqc::config::weights::WeightsV1 {
        metadata: kira_mitoqc::config::weights::Metadata {
            version: "v1.0".to_string(),
            description: "test".to_string(),
            notes: "test".to_string(),
        },
        axis: kira_mitoqc::config::weights::AxisGroupWeights {
            bioenergetics: kira_mitoqc::config::weights::AxisWeights {
                etc_stoichiometry_loss: 0.4,
                mtdna_expression_uncoupling: 0.35,
                atp_coupling_loss: 0.25,
            },
            ros: kira_mitoqc::config::weights::RosWeights {
                ros_response_overdrive: 0.6,
                nadh_imbalance: 0.4,
            },
            dynamics: kira_mitoqc::config::weights::DynamicsWeights {
                dynamics_imbalance: 0.6,
                mitophagy_excess: 0.4,
            },
            regulation: kira_mitoqc::config::weights::RegulationWeights {
                biogenesis_failure: 1.0,
            },
        },
        global: kira_mitoqc::config::weights::GlobalWeights {
            bioenergetics: 0.35,
            ros: 0.3,
            dynamics: 0.2,
            regulation: 0.15,
        },
        explainability: kira_mitoqc::config::weights::Explainability {
            max_drivers: 8,
            min_abs_contribution: 0.03,
        },
    }
}

fn fixture() -> (ExpressionSoA, GeneSet, Vec<String>) {
    let genes = vec![
        "MT-ND1", "MT-CYB", "MT-CO1", "MT-ATP6", "MT-ATP8", "NDUFS1", "SDHA", "UQCRC1", "COX4I1",
        "ATP5F1A", "SOD2", "PRDX3", "TXN2", "TXN", "GSR", "PINK1", "MFN1", "DNM1L", "TFAM",
    ];
    let mut values = Vec::new();
    let samples = 6;
    for g in 0..genes.len() {
        for s in 0..samples {
            values.push((g as f32) * 0.02 + (s as f32) * 0.05);
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
        ros_detox_genes: vecs(&["SOD2", "PRDX3"]),
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

#[test]
fn redox_metrics_deterministic() {
    let (soa, geneset, features) = fixture();
    let path = std::env::temp_dir().join("kira_mitoqc_redox_metrics_det.bin");
    write_expr_bin_with_mode(&path, &soa, ExprCacheMode::Cell).expect("write");
    let view = mmap_expr_bin(&path).expect("mmap");

    let index = GeneIndex::from_feature_list(&features);
    let resolved = resolve_all_genesets(&index, &geneset);
    let soa_index = SoaIndex::from_resolved(&resolved);
    let refs = refs_fixture();
    let weights = weights_fixture();
    let primitives =
        scalar::compute_primitives_scalar(&view, &GeneOffsets::from_resolved(&resolved));
    let proxies = compute_proxies_v1(&primitives, &resolved, &refs).expect("proxies");
    let axes = compute_axes_v1(&proxies, &weights);

    let a = compute_redox_metrics(&view, &soa_index, &axes, &proxies).expect("redox");
    let b = compute_redox_metrics(&view, &soa_index, &axes, &proxies).expect("redox");
    assert_eq!(a, b);
}

#[test]
fn redox_regime_thresholds_are_deterministic() {
    let decay = DecayScoreVec {
        decay: vec![0.0, 0.0, 0.0],
        robustness_margin: vec![1.0, 1.0, 1.0],
    };
    let axes = kira_mitoqc::score::AxisScoresVec {
        bioenergetics: vec![0.9, 0.7, 0.2],
        ros: vec![0.9, 0.7, 0.2],
        dynamics: vec![0.1, 0.1, 0.1],
        regulation: vec![0.9, 0.7, 0.2],
    };
    let redox = kira_mitoqc::redox::RedoxMetrics {
        mito_oxidative_stress_index: vec![0.9, 0.7, 0.3],
        redox_buffering_capacity: vec![0.2, 0.3, 0.7],
        mito_redox_mismatch: vec![0.7, 0.4, -0.4],
        mitochondrial_stress_adaptation_score: vec![0.9, 0.7, 0.2],
        redox_regime: vec![
            RedoxRegime::RedoxOverload,
            RedoxRegime::UnbufferedOxidativeStress,
            RedoxRegime::CompensatedOxidativeStress,
        ],
        low_confidence: vec![false, false, false],
    };
    let states = classify_v1_with_redox(&axes, &decay, &refs_fixture(), Some(&redox));
    assert_eq!(states[0].as_str(), "RedoxOverload");
    assert_eq!(states[1].as_str(), "UnbufferedOxidativeStress");
    assert_eq!(states[2].as_str(), "CompensatedOxidativeStress");
}

#[test]
fn redox_pipeline_is_scalar_dispatch_equivalent() {
    let (soa, geneset, features) = fixture();
    let path = std::env::temp_dir().join("kira_mitoqc_redox_equiv.bin");
    write_expr_bin_with_mode(&path, &soa, ExprCacheMode::Cell).expect("write");
    let view = mmap_expr_bin(&path).expect("mmap");

    let index = GeneIndex::from_feature_list(&features);
    let resolved = resolve_all_genesets(&index, &geneset);
    let soa_index = SoaIndex::from_resolved(&resolved);
    let offsets = GeneOffsets::from_resolved(&resolved);

    let refs = refs_fixture();
    let weights = weights_fixture();

    let primitives_scalar = scalar::compute_primitives_scalar(&view, &offsets);
    let primitives_dispatch = kira_mitoqc::compute::dispatch::compute_primitives(&view, &resolved);

    let proxies_scalar = compute_proxies_v1(&primitives_scalar, &resolved, &refs).expect("proxies");
    let proxies_dispatch =
        compute_proxies_v1(&primitives_dispatch, &resolved, &refs).expect("proxies");

    let axes_scalar = compute_axes_v1(&proxies_scalar, &weights);
    let axes_dispatch = compute_axes_v1(&proxies_dispatch, &weights);

    let redox_scalar =
        compute_redox_metrics(&view, &soa_index, &axes_scalar, &proxies_scalar).expect("redox");
    let redox_dispatch =
        compute_redox_metrics(&view, &soa_index, &axes_dispatch, &proxies_dispatch)
            .expect("redox");

    assert_eq!(redox_scalar, redox_dispatch);
    assert!(
        redox_scalar
            .mito_oxidative_stress_index
            .iter()
            .all(|v| *v >= 0.0 && *v <= 1.0)
    );
    assert!(
        redox_scalar
            .redox_buffering_capacity
            .iter()
            .all(|v| *v >= 0.0 && *v <= 1.0)
    );

    // Ensure key proxy path is actually in use.
    assert!(
        proxies_scalar
            .normalized
            .contains_key(&ProxyKey::RosResponseOverdrive)
    );
}
