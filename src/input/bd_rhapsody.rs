//! BD Rhapsody dense expression reader backed by kira-scio.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;
use kira_scio::error::ErrorCode;
use sprs::CsMat;

use crate::data::aggregate::{AggregationMode, ClusterMap};
use crate::input::InputError;

/// Loaded BD Rhapsody dense input.
#[derive(Debug, Clone)]
pub struct BdRhapsodyInput {
    pub matrix: CsMat<f32>,
    pub features: Vec<String>,
    pub barcodes: Vec<String>,
}

/// Resolve BD/Rhapsody-like dense file inside input path.
pub fn resolve_bd_input_path(path: &Path) -> Result<PathBuf, InputError> {
    kira_scio::resolve_bd_input_path(path).map_err(|e| map_scio_error(path, e.code, e.message))
}

/// Load BD Rhapsody dense input into the canonical sparse matrix container.
pub fn load_bd_rhapsody(path: &Path) -> Result<BdRhapsodyInput, InputError> {
    let source = resolve_bd_input_path(path)?;
    let reader = Reader::with_options(
        &source,
        ReaderOptions {
            force_format: Some(DetectedFormat::BdRhapsodyWta),
            strict: true,
        },
    );
    let canonical = reader
        .read_all()
        .map_err(|e| map_scio_error(&source, e.code, e.message))?;

    let matrix = CsMat::new_csc(
        (canonical.matrix.n_genes, canonical.matrix.n_cells),
        canonical.matrix.col_ptr,
        canonical.matrix.row_idx,
        canonical.matrix.values,
    );

    Ok(BdRhapsodyInput {
        matrix,
        features: canonical.metadata.gene_symbols,
        barcodes: canonical.metadata.barcodes,
    })
}

/// Read only feature and barcode metadata from a BD Rhapsody dense input.
pub fn load_bd_rhapsody_metadata(path: &Path) -> Result<(Vec<String>, Vec<String>), InputError> {
    let source = resolve_bd_input_path(path)?;
    let reader = Reader::with_options(
        &source,
        ReaderOptions {
            force_format: Some(DetectedFormat::BdRhapsodyWta),
            strict: true,
        },
    );
    let md = reader
        .read_metadata()
        .map_err(|e| map_scio_error(&source, e.code, e.message))?;
    Ok((md.gene_symbols, md.barcodes))
}

/// Compute mitochondrial fraction for each output sample under an aggregation mode.
pub fn compute_mito_fraction_from_file(
    path: &Path,
    mito_symbols: &BTreeSet<String>,
    mode: AggregationMode,
    clusters: Option<&ClusterMap>,
) -> Result<Vec<f32>, InputError> {
    let bd = load_bd_rhapsody(path)?;
    let samples = match mode {
        AggregationMode::Sample => 1,
        AggregationMode::Cell => bd.barcodes.len(),
        AggregationMode::Cluster => {
            let map = clusters.ok_or_else(|| InputError::InvalidClusterFile {
                path: PathBuf::from("<in-memory>"),
                message: "cluster map required for cluster mode".to_string(),
            })?;
            if map.cell_to_cluster.len() != bd.barcodes.len() {
                return Err(InputError::InvalidClusterFile {
                    path: PathBuf::from("<in-memory>"),
                    message: "cluster map length does not match BD header cell count".to_string(),
                });
            }
            map.cluster_ids.len()
        }
    };

    let mut mito_totals = vec![0.0f32; samples];
    let mut all_totals = vec![0.0f32; samples];
    let mito_rows = bd
        .features
        .iter()
        .enumerate()
        .filter_map(|(idx, name)| mito_symbols.contains(name).then_some(idx))
        .collect::<BTreeSet<_>>();

    for (cell_idx, col) in bd.matrix.outer_iterator().enumerate() {
        let out_idx = match mode {
            AggregationMode::Sample => 0,
            AggregationMode::Cell => cell_idx,
            AggregationMode::Cluster => clusters.expect("checked above").cell_to_cluster[cell_idx],
        };
        for (row_idx, value) in col.iter() {
            all_totals[out_idx] += *value;
            if mito_rows.contains(&row_idx) {
                mito_totals[out_idx] += *value;
            }
        }
    }

    Ok(all_totals
        .iter()
        .zip(mito_totals.iter())
        .map(|(all, mito)| if *all <= 0.0 { 0.0 } else { *mito / *all })
        .collect())
}

fn map_scio_error(path: &Path, code: ErrorCode, message: String) -> InputError {
    match code {
        ErrorCode::InvalidInputPath => InputError::InvalidInputPath {
            path: path.to_path_buf(),
        },
        ErrorCode::MissingFile => InputError::MissingFile {
            path: path.to_path_buf(),
        },
        ErrorCode::FeatureDisabled => InputError::H5adFeatureNotEnabled,
        _ => InputError::MatrixParse {
            path: path.to_path_buf(),
            message,
        },
    }
}
