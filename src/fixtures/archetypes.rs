//! Synthetic archetype definitions.

use crate::core::types::MitochondrialState;
use crate::fixtures::generator::{SyntheticFactors, factors_linear};

/// Fixture archetypes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Archetype {
    HealthyControl,
    AgingGradualDecline,
    RosDominantDecay,
    BioenergeticCollapse,
    MitophagyLockedDepletion,
    StructuralFragmentation,
    CompensatedButFragile,
}

/// Definition of a synthetic archetype.
#[derive(Debug, Clone)]
pub struct ArchetypeSpec {
    pub name: &'static str,
    pub samples: usize,
    pub factors: Vec<SyntheticFactors>,
    pub expected_state: MitochondrialState,
}

/// Build an archetype spec.
pub fn archetype_spec(archetype: Archetype) -> ArchetypeSpec {
    match archetype {
        Archetype::HealthyControl => ArchetypeSpec {
            name: "HealthyControl",
            samples: 1,
            factors: vec![SyntheticFactors {
                bioenergetic_damage: 0.10,
                ros_pressure: 0.10,
                dynamics_instability: 0.10,
                biogenesis_failure: 0.10,
                compensation: 0.60,
            }],
            expected_state: MitochondrialState::CompensatedButFragile,
        },
        Archetype::AgingGradualDecline => {
            let factors = factors_linear(16, |t| SyntheticFactors {
                bioenergetic_damage: 0.15 + 0.55 * t,
                ros_pressure: 0.10 + 0.45 * t,
                dynamics_instability: 0.10 + 0.35 * t,
                biogenesis_failure: 0.10 + 0.45 * t,
                compensation: 0.60 - 0.30 * t,
            });
            ArchetypeSpec {
                name: "AgingGradualDecline",
                samples: factors.len(),
                factors,
                expected_state: MitochondrialState::CompensatedButFragile,
            }
        }
        Archetype::RosDominantDecay => ArchetypeSpec {
            name: "RosDominantDecay",
            samples: 1,
            factors: vec![SyntheticFactors {
                bioenergetic_damage: 0.20,
                ros_pressure: 0.90,
                dynamics_instability: 0.40,
                biogenesis_failure: 0.20,
                compensation: 0.40,
            }],
            expected_state: MitochondrialState::RosDominantDecay,
        },
        Archetype::BioenergeticCollapse => ArchetypeSpec {
            name: "BioenergeticCollapse",
            samples: 1,
            factors: vec![SyntheticFactors {
                bioenergetic_damage: 0.90,
                ros_pressure: 0.30,
                dynamics_instability: 0.30,
                biogenesis_failure: 0.30,
                compensation: 0.15,
            }],
            expected_state: MitochondrialState::CompensatedButFragile,
        },
        Archetype::MitophagyLockedDepletion => ArchetypeSpec {
            name: "MitophagyLockedDepletion",
            samples: 1,
            factors: vec![SyntheticFactors {
                bioenergetic_damage: 0.35,
                ros_pressure: 0.35,
                dynamics_instability: 0.85,
                biogenesis_failure: 0.85,
                compensation: 0.15,
            }],
            expected_state: MitochondrialState::CompensatedButFragile,
        },
        Archetype::StructuralFragmentation => ArchetypeSpec {
            name: "StructuralFragmentation",
            samples: 1,
            factors: vec![SyntheticFactors {
                bioenergetic_damage: 0.70,
                ros_pressure: 0.40,
                dynamics_instability: 0.70,
                biogenesis_failure: 0.40,
                compensation: 0.20,
            }],
            expected_state: MitochondrialState::CompensatedButFragile,
        },
        Archetype::CompensatedButFragile => ArchetypeSpec {
            name: "CompensatedButFragile",
            samples: 1,
            factors: vec![SyntheticFactors {
                bioenergetic_damage: 0.45,
                ros_pressure: 0.45,
                dynamics_instability: 0.45,
                biogenesis_failure: 0.45,
                compensation: 0.55,
            }],
            expected_state: MitochondrialState::CompensatedButFragile,
        },
    }
}

/// List all archetypes.
pub fn all_archetypes() -> Vec<Archetype> {
    vec![
        Archetype::HealthyControl,
        Archetype::AgingGradualDecline,
        Archetype::RosDominantDecay,
        Archetype::BioenergeticCollapse,
        Archetype::MitophagyLockedDepletion,
        Archetype::StructuralFragmentation,
        Archetype::CompensatedButFragile,
    ]
}
