//! Core domain types for mitochondrial QC scoring.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// Proxy identifiers for normative scoring.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ProxyKey {
    ETCStoichiometryLoss,
    MtdnaExpressionUncoupling,
    RosResponseOverdrive,
    NadhImbalance,
    AtpCouplingLoss,
    MitophagyExcess,
    DynamicsImbalance,
    BiogenesisFailure,
}

impl ProxyKey {
    /// Stable string identifier for configuration and output.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ETCStoichiometryLoss => "ETC_stoichiometry_loss",
            Self::MtdnaExpressionUncoupling => "mtDNA_expression_uncoupling",
            Self::RosResponseOverdrive => "ROS_response_overdrive",
            Self::NadhImbalance => "NADH_imbalance",
            Self::AtpCouplingLoss => "ATP_coupling_loss",
            Self::MitophagyExcess => "mitophagy_excess",
            Self::DynamicsImbalance => "dynamics_imbalance",
            Self::BiogenesisFailure => "biogenesis_failure",
        }
    }
}

/// Axis-level aggregate scores.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AxisScores {
    pub bioenergetics: f32,
    pub ros: f32,
    pub dynamics: f32,
    pub regulation: f32,
}

/// Raw and normalized proxy scores (per-sample vectors).
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ProxyScores {
    pub normalized: BTreeMap<ProxyKey, Vec<f32>>,
    pub raw: BTreeMap<ProxyKey, Vec<f32>>,
}

impl ProxyScores {
    /// Get a normalized proxy score vector.
    pub fn get(&self, key: ProxyKey) -> Option<&[f32]> {
        self.normalized.get(&key).map(|v| v.as_slice())
    }

    /// Set a normalized proxy score vector.
    pub fn set(&mut self, key: ProxyKey, values: Vec<f32>) {
        self.normalized.insert(key, values);
    }

    /// Set a raw proxy score vector.
    pub fn set_raw(&mut self, key: ProxyKey, values: Vec<f32>) {
        self.raw.insert(key, values);
    }
}

impl Serialize for ProxyKey {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(self.as_str())
    }
}

impl<'de> Deserialize<'de> for ProxyKey {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        match s.as_str() {
            "ETC_stoichiometry_loss" => Ok(Self::ETCStoichiometryLoss),
            "mtDNA_expression_uncoupling" => Ok(Self::MtdnaExpressionUncoupling),
            "ROS_response_overdrive" => Ok(Self::RosResponseOverdrive),
            "NADH_imbalance" => Ok(Self::NadhImbalance),
            "ATP_coupling_loss" => Ok(Self::AtpCouplingLoss),
            "mitophagy_excess" => Ok(Self::MitophagyExcess),
            "dynamics_imbalance" => Ok(Self::DynamicsImbalance),
            "biogenesis_failure" => Ok(Self::BiogenesisFailure),
            _ => Err(serde::de::Error::custom(format!("unknown proxy key: {s}"))),
        }
    }
}

/// Deterministic failure-mode classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MitochondrialState {
    RosDominantDecay,
    BioenergeticCollapse,
    StructuralFragmentation,
    MitophagyLockedDepletion,
    CompensatedButFragile,
    CompensatedOxidativeStress,
    UnbufferedOxidativeStress,
    RedoxOverload,
}

impl MitochondrialState {
    /// Stable label for output.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::RosDominantDecay => "ROS-dominant decay",
            Self::BioenergeticCollapse => "Bioenergetic collapse",
            Self::StructuralFragmentation => "Structural fragmentation",
            Self::MitophagyLockedDepletion => "Mitophagy-locked depletion",
            Self::CompensatedButFragile => "Compensated but fragile",
            Self::CompensatedOxidativeStress => "CompensatedOxidativeStress",
            Self::UnbufferedOxidativeStress => "UnbufferedOxidativeStress",
            Self::RedoxOverload => "RedoxOverload",
        }
    }
}

/// Gene set collection for v1 scoring.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneSet {
    pub mtdna_complex_i: Vec<String>,
    pub mtdna_complex_iii: Vec<String>,
    pub mtdna_complex_iv: Vec<String>,
    pub mtdna_complex_v: Vec<String>,
    pub nuclear_oxphos_complex_i: Vec<String>,
    pub nuclear_oxphos_complex_ii: Vec<String>,
    pub nuclear_oxphos_complex_iii: Vec<String>,
    pub nuclear_oxphos_complex_iv: Vec<String>,
    pub nuclear_oxphos_complex_v: Vec<String>,
    pub ros_detox_genes: Vec<String>,
    pub mitophagy_genes: Vec<String>,
    pub dynamics_fusion: Vec<String>,
    pub dynamics_fission: Vec<String>,
    pub biogenesis_genes: Vec<String>,
}

impl GeneSet {
    /// Return all mtDNA genes across ETC complexes.
    pub fn all_mtdna(&self) -> Vec<&str> {
        self.mtdna_complex_i
            .iter()
            .chain(self.mtdna_complex_iii.iter())
            .chain(self.mtdna_complex_iv.iter())
            .chain(self.mtdna_complex_v.iter())
            .map(|s| s.as_str())
            .collect()
    }

    /// Return all nuclear OXPHOS anchors across ETC complexes.
    pub fn all_nuclear_oxphos(&self) -> Vec<&str> {
        self.nuclear_oxphos_complex_i
            .iter()
            .chain(self.nuclear_oxphos_complex_ii.iter())
            .chain(self.nuclear_oxphos_complex_iii.iter())
            .chain(self.nuclear_oxphos_complex_iv.iter())
            .chain(self.nuclear_oxphos_complex_v.iter())
            .map(|s| s.as_str())
            .collect()
    }
}

/// Full profile output for a sample or aggregate.
#[derive(Debug, Clone, PartialEq)]
pub struct MitoProfile {
    pub state: MitochondrialState,
    pub decay_score: f32,
    pub axis: AxisScores,
    pub robustness_margin: f32,
    pub proxies: ProxyScores,
    pub interpretation: Vec<String>,
}
