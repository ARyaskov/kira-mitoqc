//! Input abstraction layer.

use std::path::{Path, PathBuf};

use thiserror::Error;

pub mod expression;
pub mod gene_index;
pub mod spec;

pub use expression::{ExpressionContract, ExpressionUnit};
pub use gene_index::{
    GeneIndex, GeneResolution, GeneResolutionQC, ResolvedGeneSets, resolve_all_genesets,
};
pub use spec::{InputMode, InputSpec};

/// Input-layer errors.
#[derive(Debug, Error)]
pub enum InputError {
    #[error("invalid input path: {path:?}")]
    InvalidInputPath { path: PathBuf },

    #[error("feature list must be non-empty")]
    EmptyFeatureList,

    #[error("missing required file: {path:?}")]
    MissingFile { path: PathBuf },

    #[error("missing optional feature `{feature}` required to read {path:?}")]
    MissingFeature {
        feature: &'static str,
        path: PathBuf,
    },

    #[error("I/O error reading {path:?}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("matrix parse error in {path:?}: {message}")]
    MatrixParse { path: PathBuf, message: String },

    #[error(
        "dimension mismatch: matrix {rows}x{cols} vs features {features} and barcodes {barcodes}"
    )]
    DimensionMismatch {
        rows: usize,
        cols: usize,
        features: usize,
        barcodes: usize,
    },

    #[error("invalid cluster file {path:?}: {message}")]
    InvalidClusterFile { path: PathBuf, message: String },

    #[error("missing H5AD dataset: {path}")]
    MissingH5adDataset { path: String },

    #[error("unsupported H5AD matrix layout: {layout}")]
    UnsupportedH5adMatrix { layout: String },

    #[error("invalid gene symbol key: {key}")]
    InvalidGeneSymbolKey { key: String },

    #[error("h5ad feature not enabled")]
    H5adFeatureNotEnabled,

    #[error("no gene-symbol-like column detected in features.tsv")]
    GeneSymbolNotDetected,

    #[error("invalid gene symbol column: requested {requested}, available {available}")]
    InvalidGeneSymbolColumn { requested: usize, available: usize },

    #[error("invalid legacy genes.tsv: {reason}")]
    LegacyGenesTsvInvalid { reason: String },
}

/// Validate that an input path exists (no parsing performed).
pub fn validate_input_path(path: &Path) -> Result<(), InputError> {
    if path.exists() {
        Ok(())
    } else {
        Err(InputError::InvalidInputPath {
            path: path.to_path_buf(),
        })
    }
}
