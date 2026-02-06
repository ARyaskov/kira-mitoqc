use kira_mitoqc::fixtures::archetypes::{Archetype, all_archetypes, archetype_spec};
use kira_mitoqc::fixtures::golden::{run_archetype, verify_all};

#[test]
fn archetype_healthy_control() {
    let spec = archetype_spec(Archetype::HealthyControl);
    let states = run_archetype(&spec);
    for state in states {
        assert_eq!(state, spec.expected_state);
    }
}

#[test]
fn archetype_aging_gradual_decline() {
    let spec = archetype_spec(Archetype::AgingGradualDecline);
    let states = run_archetype(&spec);
    for state in states {
        assert_eq!(state, spec.expected_state);
    }
}

#[test]
fn archetype_ros_dominant_decay() {
    let spec = archetype_spec(Archetype::RosDominantDecay);
    let states = run_archetype(&spec);
    for state in states {
        assert_eq!(state, spec.expected_state);
    }
}

#[test]
fn archetype_bioenergetic_collapse() {
    let spec = archetype_spec(Archetype::BioenergeticCollapse);
    let states = run_archetype(&spec);
    for state in states {
        assert_eq!(state, spec.expected_state);
    }
}

#[test]
fn archetype_mitophagy_locked_depletion() {
    let spec = archetype_spec(Archetype::MitophagyLockedDepletion);
    let states = run_archetype(&spec);
    for state in states {
        assert_eq!(state, spec.expected_state);
    }
}

#[test]
fn archetype_structural_fragmentation() {
    let spec = archetype_spec(Archetype::StructuralFragmentation);
    let states = run_archetype(&spec);
    for state in states {
        assert_eq!(state, spec.expected_state);
    }
}

#[test]
fn archetype_compensated_but_fragile() {
    let spec = archetype_spec(Archetype::CompensatedButFragile);
    let states = run_archetype(&spec);
    for state in states {
        assert_eq!(state, spec.expected_state);
    }
}

#[test]
fn regression_all_archetypes() {
    for archetype in all_archetypes() {
        let spec = archetype_spec(archetype);
        let states = run_archetype(&spec);
        for state in states {
            assert_eq!(state, spec.expected_state);
        }
    }
    verify_all();
}
