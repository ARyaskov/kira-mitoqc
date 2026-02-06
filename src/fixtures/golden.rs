//! Golden regression expectations for synthetic archetypes.

use std::path::Path;

use crate::classify::classify_v1;
use crate::config::ConfigV1;
use crate::core::types::{MitochondrialState, ProxyKey};
use crate::explain::explain_v1;
use crate::fixtures::archetypes::{Archetype, ArchetypeSpec, all_archetypes, archetype_spec};
use crate::fixtures::generator::{SyntheticFactors, generate_primitives_v1};
use crate::input::GeneResolution;
use crate::proxy::compute_proxies_v1;
use crate::score::score_profile_v1;

/// Run the full pipeline for an archetype spec.
pub fn run_archetype(spec: &ArchetypeSpec) -> Vec<MitochondrialState> {
    let config = ConfigV1::load_from_assets_dir(Path::new("assets")).expect("load config");
    let primitives = combine_primitives(&spec.factors);
    let resolved = dummy_resolved(&config.refs);
    let proxies = compute_proxies_v1(&primitives, &resolved, &config.refs).expect("proxies");
    let scored = score_profile_v1(&proxies, &config.weights);
    let states = classify_v1(&scored.axes, &scored.decay, &config.refs);
    let _explain = explain_v1(
        &proxies,
        &scored.axes,
        &scored.decay,
        &states,
        &config.weights,
    );
    states
}

/// Verify golden expectations for all archetypes.
pub fn verify_all() {
    let config = ConfigV1::load_from_assets_dir(Path::new("assets")).expect("load config");

    for archetype in all_archetypes() {
        let spec = archetype_spec(archetype);
        let primitives = combine_primitives(&spec.factors);
        let resolved = dummy_resolved(&config.refs);
        let proxies = compute_proxies_v1(&primitives, &resolved, &config.refs).expect("proxies");
        let scored = score_profile_v1(&proxies, &config.weights);
        let states = classify_v1(&scored.axes, &scored.decay, &config.refs);
        let explain = explain_v1(
            &proxies,
            &scored.axes,
            &scored.decay,
            &states,
            &config.weights,
        );

        for state in &states {
            assert_eq!(*state, spec.expected_state, "archetype {}", spec.name);
        }

        if archetype == Archetype::AgingGradualDecline {
            assert_monotonic_increasing(&scored.decay.decay);
        }

        assert_driver_matches_intent(archetype, &explain);
    }
}

fn combine_primitives(factors: &[SyntheticFactors]) -> crate::compute::PrimitiveSignals {
    let mut combined: Option<crate::compute::PrimitiveSignals> = None;
    for factor in factors {
        let single = generate_primitives_v1(factor);
        combined = Some(match combined {
            None => single,
            Some(mut agg) => {
                agg.mtdna_mean.extend(single.mtdna_mean);
                agg.nuclear_mean.extend(single.nuclear_mean);
                agg.c_i.extend(single.c_i);
                agg.c_iii.extend(single.c_iii);
                agg.c_iv.extend(single.c_iv);
                agg.c_v.extend(single.c_v);
                agg.ros_mean.extend(single.ros_mean);
                agg.mitophagy_mean.extend(single.mitophagy_mean);
                agg.fusion_mean.extend(single.fusion_mean);
                agg.fission_mean.extend(single.fission_mean);
                agg.biogenesis_mean.extend(single.biogenesis_mean);
                agg.atp_mt.extend(single.atp_mt);
                agg.atp_nu.extend(single.atp_nu);
                agg.stoich_variance.extend(single.stoich_variance);
                agg
            }
        });
    }
    combined.expect("non-empty factors")
}

fn dummy_resolved(refs: &crate::config::refs::RefsV1) -> crate::input::ResolvedGeneSets {
    let mtdna = resolution(refs.qc.min_mtdna_genes_found + 1);
    let nuc = resolution(refs.qc.min_nuclear_oxphos_found + 1);
    let ros = resolution(refs.qc.min_ros_genes_found + 1);

    crate::input::ResolvedGeneSets {
        mtdna_complex_i: resolution(2),
        mtdna_complex_iii: resolution(2),
        mtdna_complex_iv: resolution(2),
        mtdna_complex_v: resolution(2),
        mtdna_all: mtdna,
        nuclear_oxphos_complex_i: resolution(2),
        nuclear_oxphos_complex_ii: resolution(2),
        nuclear_oxphos_complex_iii: resolution(2),
        nuclear_oxphos_complex_iv: resolution(2),
        nuclear_oxphos_complex_v: resolution(2),
        nuclear_oxphos_all: nuc,
        ros,
        mitophagy: resolution(2),
        fusion: resolution(2),
        fission: resolution(2),
        biogenesis: resolution(2),
    }
}

fn resolution(count: usize) -> GeneResolution {
    GeneResolution {
        genes: (0..count).map(|i| format!("G{i}")).collect(),
        found: (0..count).collect(),
        missing: Vec::new(),
    }
}

fn assert_monotonic_increasing(values: &[f32]) {
    assert!(!values.is_empty());
    assert!(values[values.len() - 1] >= values[0]);
}

fn assert_driver_matches_intent(archetype: Archetype, explain: &[crate::explain::Explainability]) {
    let expected = match archetype {
        Archetype::HealthyControl => None,
        Archetype::AgingGradualDecline => Some(ProxyKey::BiogenesisFailure),
        Archetype::RosDominantDecay => Some(ProxyKey::RosResponseOverdrive),
        Archetype::BioenergeticCollapse => Some(ProxyKey::BiogenesisFailure),
        Archetype::MitophagyLockedDepletion => Some(ProxyKey::BiogenesisFailure),
        Archetype::StructuralFragmentation => Some(ProxyKey::BiogenesisFailure),
        Archetype::CompensatedButFragile => Some(ProxyKey::BiogenesisFailure),
    };

    if let Some(expected) = expected {
        let idx = if archetype == Archetype::AgingGradualDecline {
            explain.len() - 1
        } else {
            0
        };
        let top = explain[idx].drivers.first().expect("no drivers").key;
        assert_eq!(top, expected, "archetype intent mismatch");
    }
}
