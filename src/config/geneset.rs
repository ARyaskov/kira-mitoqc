//! Geneset configuration parsing.

use serde::Deserialize;

use crate::core::types::GeneSet;

#[derive(Debug, Clone, Deserialize)]
pub(crate) struct Metadata {
    #[allow(dead_code)]
    pub version: String,
    #[allow(dead_code)]
    pub description: String,
    #[allow(dead_code)]
    pub source: String,
    #[allow(dead_code)]
    pub notes: String,
}

#[derive(Debug, Clone, Deserialize)]
pub(crate) struct Mtdna {
    #[serde(rename = "complex_I")]
    pub complex_i: Vec<String>,
    #[serde(rename = "complex_III")]
    pub complex_iii: Vec<String>,
    #[serde(rename = "complex_IV")]
    pub complex_iv: Vec<String>,
    #[serde(rename = "complex_V")]
    pub complex_v: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
pub(crate) struct NuclearOxphos {
    #[serde(rename = "complex_I")]
    pub complex_i: Vec<String>,
    #[serde(rename = "complex_II")]
    pub complex_ii: Vec<String>,
    #[serde(rename = "complex_III")]
    pub complex_iii: Vec<String>,
    #[serde(rename = "complex_IV")]
    pub complex_iv: Vec<String>,
    #[serde(rename = "complex_V")]
    pub complex_v: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
pub(crate) struct GenesOnly {
    pub genes: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
pub(crate) struct Dynamics {
    pub fusion: Vec<String>,
    pub fission: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct GenesetV1 {
    #[allow(dead_code)]
    pub(crate) metadata: Metadata,
    #[serde(rename = "mtDNA")]
    pub(crate) mtdna: Mtdna,
    pub(crate) nuclear_oxphos: NuclearOxphos,
    pub(crate) ros_detox: GenesOnly,
    pub(crate) mitophagy: GenesOnly,
    pub(crate) dynamics: Dynamics,
    pub(crate) biogenesis: GenesOnly,
}

impl GenesetV1 {
    /// Convert to the core GeneSet struct.
    pub fn into_geneset(self) -> GeneSet {
        GeneSet {
            mtdna_complex_i: self.mtdna.complex_i,
            mtdna_complex_iii: self.mtdna.complex_iii,
            mtdna_complex_iv: self.mtdna.complex_iv,
            mtdna_complex_v: self.mtdna.complex_v,
            nuclear_oxphos_complex_i: self.nuclear_oxphos.complex_i,
            nuclear_oxphos_complex_ii: self.nuclear_oxphos.complex_ii,
            nuclear_oxphos_complex_iii: self.nuclear_oxphos.complex_iii,
            nuclear_oxphos_complex_iv: self.nuclear_oxphos.complex_iv,
            nuclear_oxphos_complex_v: self.nuclear_oxphos.complex_v,
            ros_detox_genes: self.ros_detox.genes,
            mitophagy_genes: self.mitophagy.genes,
            dynamics_fusion: self.dynamics.fusion,
            dynamics_fission: self.dynamics.fission,
            biogenesis_genes: self.biogenesis.genes,
        }
    }
}
