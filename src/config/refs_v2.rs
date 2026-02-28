//! v2 reference constants for multi-omics.

use std::fs;
use std::path::Path;

use serde::Deserialize;

use crate::config::ConfigError;

const EMBEDDED_REFS_V2: &str =
    include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/assets/refs_v2.toml"));

/// v2 reference constants.
#[derive(Debug, Clone, Deserialize)]
pub struct RefsV2Consts {
    pub cn_ref: f32,
    pub het_ref: f32,
    pub del_ref: f32,
    pub prot_stoich_ref: f32,
    pub prot_atp_ref: f32,
}

/// Mixing coefficients for RNA-protein correction.
#[derive(Debug, Clone, Deserialize)]
pub struct Mixing {
    pub alpha_rna_protein_stoich: f32,
    pub alpha_rna_protein_atp: f32,
}

/// v2 references bundle.
#[derive(Debug, Clone, Deserialize)]
pub struct RefsV2 {
    #[serde(rename = "refs_v2")]
    pub refs: RefsV2Consts,
    pub mixing: Mixing,
}

/// Load refs_v2.toml from a path.
pub fn load_refs_v2(path: &Path) -> Result<RefsV2, ConfigError> {
    let contents = fs::read_to_string(path).map_err(|source| ConfigError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    toml::from_str(&contents).map_err(|source| ConfigError::Toml {
        path: path.to_path_buf(),
        source,
    })
}

pub fn load_refs_v2_embedded() -> Result<RefsV2, ConfigError> {
    toml::from_str(EMBEDDED_REFS_V2).map_err(|source| ConfigError::Toml {
        path: std::path::PathBuf::from("embedded://refs_v2.toml"),
        source,
    })
}
