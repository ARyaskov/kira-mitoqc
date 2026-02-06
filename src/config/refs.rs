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
pub struct Refs {
    #[serde(rename = "STOICH_REF")]
    pub stoich_ref: f32,
    #[serde(rename = "UNCOUPLING_REF")]
    pub uncoupling_ref: f32,
    #[serde(rename = "ROS_REF")]
    pub ros_ref: f32,
    #[serde(rename = "REDOX_REF")]
    pub redox_ref: f32,
    #[serde(rename = "ATP_REF")]
    pub atp_ref: f32,
    #[serde(rename = "MITO_REF")]
    pub mito_ref: f32,
    #[serde(rename = "DYN_REF")]
    pub dyn_ref: f32,
    #[serde(rename = "BIO_REF")]
    pub bio_ref: f32,
}

/// Thresholds for classification.
#[derive(Debug, Clone, Deserialize)]
pub struct Thresholds {
    #[serde(rename = "ROS_HIGH")]
    pub ros_high: f32,
    #[serde(rename = "BIOENERGETICS_LOW")]
    pub bioenergetics_low: f32,
    #[serde(rename = "BIOENERGETICS_HIGH")]
    pub bioenergetics_high: f32,
    #[serde(rename = "DYNAMICS_HIGH")]
    pub dynamics_high: f32,
    #[serde(rename = "REGULATION_LOW")]
    pub regulation_low: f32,
    #[serde(rename = "STRUCTURAL_BIO_MIN")]
    pub structural_bio_min: f32,
    #[serde(rename = "STRUCTURAL_DYN_MIN")]
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
