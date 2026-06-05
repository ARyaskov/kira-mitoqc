use kira_mitoqc::compute::PrimitiveSignals;
use kira_mitoqc::config::refs::{Eps, Metadata, Normalization, Qc, Refs, RefsV1, Thresholds};
use kira_mitoqc::config::refs_v2::{Mixing, RefsV2, RefsV2Consts};
use kira_mitoqc::config::weights::{DynamicsWeights, Explainability, GlobalWeights, RosWeights};
use kira_mitoqc::config::weights_v2::{
    AxisBioenergeticsV2, AxisGroupWeightsV2, AxisRegulationV2, WeightsV2,
};
use kira_mitoqc::input::{GeneResolution, ResolvedGeneSets};
use kira_mitoqc::proxy::{OptionalOmicsInputs, ProxyKeyV2, compute_proxies_v1, compute_proxies_v2};
use kira_mitoqc::score::{compute_axes_v2, compute_decay_v2};

fn refs_v1_fixture() -> RefsV1 {
    RefsV1 {
        metadata: Metadata {
            version: "v1".to_string(),
            description: "".to_string(),
            notes: "".to_string(),
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
            min_mtdna_genes_found: 1,
            min_nuclear_oxphos_found: 1,
            min_ros_genes_found: 1,
        },
    }
}

fn dummy_resolved() -> ResolvedGeneSets {
    let res = |n| GeneResolution {
        genes: vec!["G".to_string(); n],
        found: (0..n).collect(),
        missing: vec![],
    };
    ResolvedGeneSets {
        mtdna_complex_i: res(1),
        mtdna_complex_iii: res(1),
        mtdna_complex_iv: res(1),
        mtdna_complex_v: res(1),
        mtdna_all: res(2),
        nuclear_oxphos_complex_i: res(1),
        nuclear_oxphos_complex_ii: res(1),
        nuclear_oxphos_complex_iii: res(1),
        nuclear_oxphos_complex_iv: res(1),
        nuclear_oxphos_complex_v: res(1),
        nuclear_oxphos_all: res(2),
        ros: res(1),
        mitophagy: res(1),
        fusion: res(1),
        fission: res(1),
        biogenesis: res(1),
    }
}

fn primitives_fixture() -> PrimitiveSignals {
    PrimitiveSignals {
        mtdna_mean: vec![0.2, 0.3],
        nuclear_mean: vec![0.2, 0.3],
        c_i: vec![0.2, 0.3],
        c_iii: vec![0.2, 0.3],
        c_iv: vec![0.2, 0.3],
        c_v: vec![0.2, 0.3],
        ros_mean: vec![0.2, 0.3],
        mitophagy_mean: vec![0.2, 0.3],
        fusion_mean: vec![0.2, 0.3],
        fission_mean: vec![0.2, 0.3],
        biogenesis_mean: vec![0.2, 0.3],
        atp_mt: vec![0.2, 0.3],
        atp_nu: vec![0.2, 0.3],
        stoich_variance: vec![0.1, 0.2],
    }
}

fn refs_v2_fixture() -> RefsV2 {
    RefsV2 {
        refs: RefsV2Consts {
            cn_ref: 1.0,
            het_ref: 1.0,
            del_ref: 1.0,
            prot_stoich_ref: 1.0,
            prot_atp_ref: 1.0,
        },
        mixing: Mixing {
            alpha_rna_protein_stoich: 0.6,
            alpha_rna_protein_atp: 0.6,
        },
    }
}

fn weights_v2_fixture() -> WeightsV2 {
    WeightsV2 {
        axis: AxisGroupWeightsV2 {
            bioenergetics: AxisBioenergeticsV2 {
                etc_stoichiometry_loss: 0.32,
                mtdna_expression_uncoupling: 0.28,
                atp_coupling_loss: 0.20,
                mtdna_copy_number_instability: 0.10,
                mtdna_heteroplasmy_burden: 0.10,
            },
            ros: RosWeights {
                ros_response_overdrive: 0.6,
                nadh_imbalance: 0.4,
            },
            dynamics: DynamicsWeights {
                dynamics_imbalance: 0.6,
                mitophagy_excess: 0.4,
            },
            regulation: AxisRegulationV2 {
                biogenesis_failure: 0.70,
                mtdna_copy_number_instability: 0.15,
                mtdna_heteroplasmy_burden: 0.15,
            },
        },
        global: GlobalWeights {
            bioenergetics: 0.35,
            ros: 0.3,
            dynamics: 0.2,
            regulation: 0.15,
        },
        explainability: Explainability {
            max_drivers: 8,
            min_abs_contribution: 0.03,
        },
    }
}

#[test]
fn v2_disabled_matches_v1_proxies() {
    let refs = refs_v1_fixture();
    let resolved = dummy_resolved();
    let primitives = primitives_fixture();
    let v1 = compute_proxies_v1(&primitives, &resolved, &refs).expect("v1");
    let v1_clone = v1.clone();
    let v2 = compute_proxies_v2(
        &primitives,
        v1,
        &refs_v2_fixture(),
        &OptionalOmicsInputs::default(),
    );
    assert_eq!(v2.v1, v1_clone);
    assert!(v2.v2_normalized.is_empty());
}

#[test]
fn mixing_corrections_apply() {
    let refs = refs_v1_fixture();
    let resolved = dummy_resolved();
    let primitives = primitives_fixture();
    let v1 = compute_proxies_v1(&primitives, &resolved, &refs).expect("v1");
    let extra = OptionalOmicsInputs {
        mt_dna_copy_number: None,
        heteroplasmy: None,
        mt_dna_deletions: None,
        proteomics_etc: Some(vec![0.9, 0.9]),
        proteomics_atp: Some(vec![0.8, 0.8]),
    };
    let refs_v2 = refs_v2_fixture();
    let v2 = compute_proxies_v2(&primitives, v1, &refs_v2, &extra);

    let axes_v2 = compute_axes_v2(&v2, &weights_v2_fixture(), &refs_v2);
    let decay_v2 = compute_decay_v2(&axes_v2, &weights_v2_fixture());

    assert_eq!(axes_v2.bioenergetics.len(), 2);
    assert_eq!(decay_v2.decay.len(), 2);

    let prot = v2
        .v2_normalized
        .get(&ProxyKeyV2::ProteomicsEtcStoichiometryLoss)
        .unwrap();
    assert!((prot[0] - 0.9).abs() < 1e-6);
}
