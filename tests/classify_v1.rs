use kira_mitoqc::classify::{ClassifiedProfile, classify_v1};
use kira_mitoqc::config::refs::{Eps, Metadata, Normalization, Qc, Refs, RefsV1, Thresholds};
use kira_mitoqc::core::types::MitochondrialState;
use kira_mitoqc::score::{AxisScoresVec, DecayScoreVec};

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
            min_mtdna_genes_found: 1,
            min_nuclear_oxphos_found: 1,
            min_ros_genes_found: 1,
        },
    }
}

fn axes_for(samples: &[([f32; 4])]) -> AxisScoresVec {
    let mut bio = Vec::new();
    let mut ros = Vec::new();
    let mut dyns = Vec::new();
    let mut reg = Vec::new();
    for [b, r, d, g] in samples {
        bio.push(*b);
        ros.push(*r);
        dyns.push(*d);
        reg.push(*g);
    }
    AxisScoresVec {
        bioenergetics: bio,
        ros,
        dynamics: dyns,
        regulation: reg,
    }
}

#[test]
fn rule_ordering_matches_spec() {
    let refs = refs_fixture();
    let axes = axes_for(&[
        [0.4, 0.8, 0.2, 0.5],
        [0.9, 0.2, 0.2, 0.2],
        [0.5, 0.2, 0.8, 0.2],
        [0.7, 0.2, 0.7, 0.7],
        [0.2, 0.2, 0.2, 0.2],
    ]);
    let decay = DecayScoreVec {
        decay: vec![0.0; 5],
        robustness_margin: vec![1.0; 5],
    };

    let states = classify_v1(&axes, &decay, &refs);

    assert_eq!(states[0], MitochondrialState::RosDominantDecay);
    assert_eq!(states[1], MitochondrialState::BioenergeticCollapse);
    assert_eq!(states[2], MitochondrialState::MitophagyLockedDepletion);
    assert_eq!(states[3], MitochondrialState::StructuralFragmentation);
    assert_eq!(states[4], MitochondrialState::CompensatedButFragile);
}

#[test]
fn threshold_edges_are_exclusive() {
    let refs = refs_fixture();
    let axes = axes_for(&[[0.5, 0.75, 0.7, 0.4]]);
    let decay = DecayScoreVec {
        decay: vec![0.0],
        robustness_margin: vec![1.0],
    };

    let states = classify_v1(&axes, &decay, &refs);
    assert_eq!(states[0], MitochondrialState::CompensatedButFragile);
}

#[test]
fn classified_profile_wrapper() {
    let refs = refs_fixture();
    let axes = axes_for(&[[0.9, 0.2, 0.2, 0.2]]);
    let decay = DecayScoreVec {
        decay: vec![0.0],
        robustness_margin: vec![1.0],
    };

    let profile = ClassifiedProfile {
        states: classify_v1(&axes, &decay, &refs),
    };
    assert_eq!(profile.states[0], MitochondrialState::BioenergeticCollapse);
}
