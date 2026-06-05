//! Reference constants and thresholds for v1.

use serde::Deserialize;

/// Metadata for configuration provenance.
#[derive(Debug, Clone, Deserialize)]
pub struct Metadata {
    #[allow(dead_code)]
    pub version: String,
    #[allow(dead_code)]
    pub description: String,
    #[allow(dead_code)]
    pub notes: String,
}

/// Epsilon configuration.
#[derive(Debug, Clone, Deserialize)]
pub struct Eps {
    pub value: f32,
}

/// Reference constants for normalization.
#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct Refs {
    pub stoich_ref: f32,
    pub uncoupling_ref: f32,
    pub ros_ref: f32,
    pub redox_ref: f32,
    pub atp_ref: f32,
    pub mito_ref: f32,
    pub dyn_ref: f32,
    pub bio_ref: f32,
}

/// Thresholds for classification.
#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub struct Thresholds {
    pub ros_high: f32,
    pub bioenergetics_low: f32,
    pub bioenergetics_high: f32,
    pub dynamics_high: f32,
    pub regulation_low: f32,
    pub structural_bio_min: f32,
    pub structural_dyn_min: f32,
}

/// QC thresholds.
#[derive(Debug, Clone, Deserialize)]
pub struct Qc {
    pub min_mtdna_genes_found: usize,
    pub min_nuclear_oxphos_found: usize,
    pub min_ros_genes_found: usize,
}

/// Normalization assumptions for input data.
#[derive(Debug, Clone, Deserialize)]
pub struct Normalization {
    pub log1p: bool,
    pub expression_unit: String,
}

/// References configuration schema for v1.
#[derive(Debug, Clone, Deserialize)]
pub struct RefsV1 {
    #[allow(dead_code)]
    pub metadata: Metadata,
    pub eps: Eps,
    pub refs: Refs,
    pub thresholds: Thresholds,
    pub normalization: Normalization,
    pub qc: Qc,
}
