use kira_mitoqc::compute::PrimitiveSignals;
use kira_mitoqc::config::refs::{Eps, Metadata, Normalization, Qc, Refs, RefsV1, Thresholds};
use kira_mitoqc::core::types::{ProxyKey, ProxyScores};
use kira_mitoqc::input::{GeneResolution, ResolvedGeneSets};
use kira_mitoqc::proxy::{ProxyError, compute_proxies_v1};

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
            expression_unit: "log1p_cpm_or_tpm".to_string(),
        },
        qc: Qc {
            min_mtdna_genes_found: 2,
            min_nuclear_oxphos_found: 2,
            min_ros_genes_found: 1,
        },
    }
}

fn resolution(found: usize) -> GeneResolution {
    GeneResolution {
        genes: vec!["G".to_string(); found],
        found: (0..found).collect(),
        missing: Vec::new(),
    }
}

fn resolved_fixture(mtdna: usize, nuc: usize, ros: usize) -> ResolvedGeneSets {
    ResolvedGeneSets {
        mtdna_complex_i: resolution(1),
        mtdna_complex_iii: resolution(1),
        mtdna_complex_iv: resolution(0),
        mtdna_complex_v: resolution(0),
        mtdna_all: resolution(mtdna),
        nuclear_oxphos_complex_i: resolution(1),
        nuclear_oxphos_complex_ii: resolution(1),
        nuclear_oxphos_complex_iii: resolution(0),
        nuclear_oxphos_complex_iv: resolution(0),
        nuclear_oxphos_complex_v: resolution(0),
        nuclear_oxphos_all: resolution(nuc),
        ros: resolution(ros),
        mitophagy: resolution(1),
        fusion: resolution(1),
        fission: resolution(1),
        biogenesis: resolution(1),
    }
}

fn primitives_fixture() -> PrimitiveSignals {
    PrimitiveSignals {
        mtdna_mean: vec![0.2],
        nuclear_mean: vec![0.6],
        c_i: vec![0.5],
        c_iii: vec![0.5],
        c_iv: vec![0.5],
        c_v: vec![0.5],
        ros_mean: vec![0.4],
        mitophagy_mean: vec![0.3],
        fusion_mean: vec![0.8],
        fission_mean: vec![0.4],
        biogenesis_mean: vec![0.2],
        atp_mt: vec![0.7],
        atp_nu: vec![0.5],
        stoich_variance: vec![0.05],
    }
}

fn get(scores: &ProxyScores, key: ProxyKey) -> f32 {
    scores.normalized.get(&key).expect("missing key")[0]
}

#[test]
fn proxy_formulas_match_spec() {
    let refs = refs_fixture();
    let resolved = resolved_fixture(2, 2, 1);
    let primitives = primitives_fixture();

    let scores = compute_proxies_v1(&primitives, &resolved, &refs).expect("compute");

    assert!((get(&scores, ProxyKey::ETCStoichiometryLoss) - 0.2).abs() < 1e-6);
    assert!((get(&scores, ProxyKey::MtdnaExpressionUncoupling) - 0.4).abs() < 1e-6);
    assert!((get(&scores, ProxyKey::RosResponseOverdrive) - 0.4).abs() < 1e-6);
    let redox_proxy = primitives.c_i[0] / (primitives.ros_mean[0] + refs.eps.value);
    let redox_norm = (1.0 - redox_proxy).clamp(0.0, 1.0);
    assert!((get(&scores, ProxyKey::NadhImbalance) - redox_norm).abs() < 1e-6);
    let atp_raw = 1.0 - (primitives.atp_mt[0] / (primitives.atp_nu[0] + refs.eps.value));
    let atp_norm = atp_raw.clamp(0.0, 1.0);
    assert!((get(&scores, ProxyKey::AtpCouplingLoss) - atp_norm).abs() < 1e-6);
    let dynamics_raw = (primitives.fusion_mean[0] / (primitives.fission_mean[0] + refs.eps.value))
        .log2()
        .abs();
    let dynamics_norm = dynamics_raw.clamp(0.0, 1.0);
    assert!((get(&scores, ProxyKey::DynamicsImbalance) - dynamics_norm).abs() < 1e-6);
    let bio_norm = (1.0 - primitives.biogenesis_mean[0]).clamp(0.0, 1.0);
    assert!((get(&scores, ProxyKey::BiogenesisFailure) - bio_norm).abs() < 1e-6);
}

#[test]
fn qc_failure_is_error_for_critical_sets() {
    let refs = refs_fixture();
    let resolved = resolved_fixture(1, 2, 1);
    let primitives = primitives_fixture();

    let err = compute_proxies_v1(&primitives, &resolved, &refs).unwrap_err();
    match err {
        ProxyError::InsufficientGenes { set, .. } => assert_eq!(set, "mtDNA"),
        _ => panic!("unexpected error: {err:?}"),
    }
}

#[test]
fn clamp_behavior_bounds() {
    let refs = refs_fixture();
    let resolved = resolved_fixture(2, 2, 1);
    let mut primitives = primitives_fixture();
    primitives.ros_mean = vec![2.5];

    let scores = compute_proxies_v1(&primitives, &resolved, &refs).expect("compute");
    assert_eq!(get(&scores, ProxyKey::RosResponseOverdrive), 1.0);
}
