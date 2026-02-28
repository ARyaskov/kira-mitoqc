//! Matrix Market (MTX) loader backed by kira-scio.

use std::path::{Path, PathBuf};

use kira_scio::api::Reader;
use kira_scio::error::ErrorCode;
use sprs::CsMat;
use tracing::info;

use crate::input::InputError;

/// Loaded MTX inputs.
#[derive(Debug, Clone)]
pub struct MtxInput {
    pub matrix: CsMat<f32>,
    pub features: Vec<String>,
    pub barcodes: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureFileKind {
    FeaturesTsv,
    GenesTsv,
}

#[derive(Debug, Clone)]
pub struct MtxDiscovery {
    pub input_dir: PathBuf,
    pub prefix: Option<String>,
    pub matrix_path: PathBuf,
    pub feature_path: PathBuf,
    pub feature_kind: FeatureFileKind,
    pub barcodes_path: PathBuf,
}

/// Load features.tsv/genes.tsv and barcodes.tsv without parsing the matrix.
pub fn load_mtx_metadata(
    path: &Path,
    gene_symbol_col: Option<usize>,
) -> Result<(Vec<String>, Vec<String>), InputError> {
    if let Some(col) = gene_symbol_col
        && col != 1
    {
        return Err(InputError::InvalidGeneSymbolColumn {
            requested: col + 1,
            available: 2,
        });
    }

    let metadata = Reader::new(path)
        .read_metadata()
        .map_err(|e| map_scio_error(path, e.code, e.message))?;
    Ok((metadata.gene_symbols, metadata.barcodes))
}

/// Load Matrix Market files from a directory.
pub fn load_mtx_dir(path: &Path, gene_symbol_col: Option<usize>) -> Result<MtxInput, InputError> {
    if let Some(col) = gene_symbol_col
        && col != 1
    {
        return Err(InputError::InvalidGeneSymbolColumn {
            requested: col + 1,
            available: 2,
        });
    }

    info!(path = ?path, "Loading MTX directory via kira-scio");
    let canonical = Reader::new(path)
        .read_all()
        .map_err(|e| map_scio_error(path, e.code, e.message))?;

    let matrix = CsMat::new_csc(
        (canonical.matrix.n_genes, canonical.matrix.n_cells),
        canonical.matrix.col_ptr,
        canonical.matrix.row_idx,
        canonical.matrix.values,
    );

    Ok(MtxInput {
        matrix,
        features: canonical.metadata.gene_symbols,
        barcodes: canonical.metadata.barcodes,
    })
}

/// Discover MTX/feature/barcode files in an input directory, including optional dataset prefix.
pub fn discover_dataset_files(path: &Path) -> Result<MtxDiscovery, InputError> {
    let ds = kira_scio::discover(path).map_err(|e| map_scio_error(path, e.code, e.message))?;

    let (feature_path, feature_kind) = if let Some(features) = ds.features {
        info!("Detected 10x v3+ format: using features.tsv");
        (features, FeatureFileKind::FeaturesTsv)
    } else if let Some(genes) = ds.genes {
        info!("Detected legacy 10x v2 format: using genes.tsv (column 2 as gene symbol)");
        (genes, FeatureFileKind::GenesTsv)
    } else {
        return Err(InputError::MissingFile {
            path: ds.input_dir.join("features.tsv"),
        });
    };

    let barcodes_path = ds.barcodes.ok_or_else(|| InputError::MissingFile {
        path: ds.input_dir.join("barcodes.tsv"),
    })?;

    Ok(MtxDiscovery {
        input_dir: ds.input_dir,
        prefix: ds.prefix,
        matrix_path: ds.matrix,
        feature_path,
        feature_kind,
        barcodes_path,
    })
}

/// Resolve shared pipeline cache filename from optional dataset prefix.
pub fn resolve_shared_cache_filename(prefix: Option<&str>) -> String {
    kira_scio::resolve_shared_cache_filename(prefix)
}

fn map_scio_error(path: &Path, code: ErrorCode, message: String) -> InputError {
    match code {
        ErrorCode::InvalidInputPath => InputError::InvalidInputPath {
            path: path.to_path_buf(),
        },
        ErrorCode::MissingFile => InputError::MissingFile {
            path: path.to_path_buf(),
        },
        ErrorCode::DimensionMismatch => InputError::DimensionMismatch {
            rows: 0,
            cols: 0,
            features: 0,
            barcodes: 0,
        },
        ErrorCode::FeatureDisabled => InputError::H5adFeatureNotEnabled,
        ErrorCode::Io => InputError::MatrixParse {
            path: path.to_path_buf(),
            message,
        },
        _ => InputError::MatrixParse {
            path: path.to_path_buf(),
            message,
        },
    }
}
