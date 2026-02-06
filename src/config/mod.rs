//! Configuration loading and validation for v1.

use std::fs;
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::core::types::GeneSet;

pub mod geneset;
pub mod refs;
pub mod refs_v2;
pub mod weights;
pub mod weights_v2;

use geneset::GenesetV1;
use refs::RefsV1;
use weights::WeightsV1;

/// Fully loaded configuration bundle for v1.
#[derive(Debug, Clone)]
pub struct ConfigV1 {
    pub geneset: GeneSet,
    pub weights: WeightsV1,
    pub refs: RefsV1,
}

impl ConfigV1 {
    /// Load geneset, weights, and refs from a directory containing TOML assets.
    pub fn load_from_assets_dir(path: &Path) -> Result<Self, ConfigError> {
        let geneset_path = path.join("geneset_v1.toml");
        let weights_path = path.join("weights_v1.toml");
        let refs_path = path.join("refs_v1.toml");

        let geneset_toml = read_to_string(&geneset_path)?;
        let weights_toml = read_to_string(&weights_path)?;
        let refs_toml = read_to_string(&refs_path)?;

        let geneset_v1: GenesetV1 = parse_toml(&geneset_path, &geneset_toml)?;
        let weights_v1: WeightsV1 = parse_toml(&weights_path, &weights_toml)?;
        let refs_v1: RefsV1 = parse_toml(&refs_path, &refs_toml)?;

        let geneset = geneset_v1.into_geneset();
        validate_geneset(&geneset)?;
        weights_v1.validate(1e-6).map_err(ConfigError::Validation)?;

        Ok(Self {
            geneset,
            weights: weights_v1,
            refs: refs_v1,
        })
    }
}

fn read_to_string(path: &Path) -> Result<String, ConfigError> {
    fs::read_to_string(path).map_err(|source| ConfigError::Io {
        path: path.to_path_buf(),
        source,
    })
}

fn parse_toml<T>(path: &Path, contents: &str) -> Result<T, ConfigError>
where
    T: serde::de::DeserializeOwned,
{
    toml::from_str(contents).map_err(|source| ConfigError::Toml {
        path: path.to_path_buf(),
        source,
    })
}

fn validate_geneset(geneset: &GeneSet) -> Result<(), ConfigError> {
    let checks: [(&str, usize); 14] = [
        ("mtdna_complex_i", geneset.mtdna_complex_i.len()),
        ("mtdna_complex_iii", geneset.mtdna_complex_iii.len()),
        ("mtdna_complex_iv", geneset.mtdna_complex_iv.len()),
        ("mtdna_complex_v", geneset.mtdna_complex_v.len()),
        (
            "nuclear_oxphos_complex_i",
            geneset.nuclear_oxphos_complex_i.len(),
        ),
        (
            "nuclear_oxphos_complex_ii",
            geneset.nuclear_oxphos_complex_ii.len(),
        ),
        (
            "nuclear_oxphos_complex_iii",
            geneset.nuclear_oxphos_complex_iii.len(),
        ),
        (
            "nuclear_oxphos_complex_iv",
            geneset.nuclear_oxphos_complex_iv.len(),
        ),
        (
            "nuclear_oxphos_complex_v",
            geneset.nuclear_oxphos_complex_v.len(),
        ),
        ("ros_detox_genes", geneset.ros_detox_genes.len()),
        ("mitophagy_genes", geneset.mitophagy_genes.len()),
        ("dynamics_fusion", geneset.dynamics_fusion.len()),
        ("dynamics_fission", geneset.dynamics_fission.len()),
        ("biogenesis_genes", geneset.biogenesis_genes.len()),
    ];

    for (name, len) in checks {
        if len == 0 {
            return Err(ConfigError::Validation(format!(
                "geneset list {name} must be non-empty"
            )));
        }
    }

    Ok(())
}

/// Configuration loading errors.
#[derive(Debug, Error)]
pub enum ConfigError {
    #[error("I/O error reading {path:?}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("TOML parse error in {path:?}: {source}")]
    Toml {
        path: PathBuf,
        source: toml::de::Error,
    },

    #[error("Validation error: {0}")]
    Validation(String),
}
