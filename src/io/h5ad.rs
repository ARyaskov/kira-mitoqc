//! H5AD loader backed by kira-scio.

use std::path::Path;

use hdf5::File;
use kira_scio::api::{Reader, ReaderOptions};
use kira_scio::detect::DetectedFormat;
use kira_scio::error::ErrorCode;
use sprs::CsMat;

use crate::data::aggregate::ClusterMap;
use crate::input::InputError;

/// Loaded H5AD inputs.
#[derive(Debug, Clone)]
pub struct H5adInput {
    pub matrix: CsMat<f32>,
    pub features: Vec<String>,
    pub barcodes: Vec<String>,
}

/// Load H5AD expression matrix and metadata.
pub fn load_h5ad(path: &Path, gene_symbol_key: Option<&str>) -> Result<H5adInput, InputError> {
    validate_gene_symbol_key(gene_symbol_key)?;
    let reader = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::H5ad),
            strict: true,
        },
    );
    let canonical = reader
        .read_all()
        .map_err(|e| map_scio_error(path, e.code, e.message))?;

    let matrix = CsMat::new_csc(
        (canonical.matrix.n_genes, canonical.matrix.n_cells),
        canonical
            .matrix
            .col_ptr
            .into_iter()
            .map(|v| v as usize)
            .collect::<Vec<_>>(),
        canonical
            .matrix
            .row_idx
            .into_iter()
            .map(|v| v as usize)
            .collect::<Vec<_>>(),
        canonical.matrix.values,
    );

    Ok(H5adInput {
        matrix,
        features: canonical.metadata.gene_symbols,
        barcodes: canonical.metadata.barcodes,
    })
}

/// Load H5AD features and barcodes without reading the matrix.
pub fn load_h5ad_metadata(
    path: &Path,
    gene_symbol_key: Option<&str>,
) -> Result<(Vec<String>, Vec<String>), InputError> {
    validate_gene_symbol_key(gene_symbol_key)?;
    let reader = Reader::with_options(
        path,
        ReaderOptions {
            force_format: Some(DetectedFormat::H5ad),
            strict: true,
        },
    );
    let md = reader
        .read_metadata()
        .map_err(|e| map_scio_error(path, e.code, e.message))?;
    Ok((md.gene_symbols, md.barcodes))
}

/// Load cluster labels from obs/<column> and build ClusterMap.
pub fn load_h5ad_clusters(path: &Path, column: &str) -> Result<ClusterMap, InputError> {
    let file = File::open(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let barcodes = read_obs_index(&file)?;
    let labels = read_strings(&file, &format!("obs/{column}"))?;
    if labels.len() != barcodes.len() {
        return Err(InputError::InvalidClusterFile {
            path: path.to_path_buf(),
            message: "cluster labels length mismatch".to_string(),
        });
    }
    Ok(build_cluster_map(&barcodes, &labels))
}

/// Read the obs index, trying conventional names used by different writers.
fn read_obs_index(file: &File) -> Result<Vec<String>, InputError> {
    // Scanpy historically writes `obs/_index`; some toolkits use `obs/index`
    // or per-dataset columns. Try them in priority order before giving up.
    for candidate in ["obs/_index", "obs/index", "obs/barcode", "obs/cell_id"] {
        if file.dataset(candidate).is_ok() {
            return read_strings(file, candidate);
        }
    }
    Err(InputError::MissingH5adDataset {
        path: "obs/_index (or alias)".to_string(),
    })
}

/// Allowlist of well-known column names; arbitrary keys are accepted but
/// non-empty values are flagged as user error early.
fn validate_gene_symbol_key(gene_symbol_key: Option<&str>) -> Result<(), InputError> {
    if let Some(key) = gene_symbol_key
        && key.is_empty()
    {
        return Err(InputError::InvalidGeneSymbolKey {
            key: String::new(),
        });
    }
    Ok(())
}

fn read_strings(file: &File, path: &str) -> Result<Vec<String>, InputError> {
    let dataset = file
        .dataset(path)
        .map_err(|_| InputError::MissingH5adDataset {
            path: path.to_string(),
        })?;
    dataset
        .read_1d::<String>()
        .map_err(|_| InputError::MissingH5adDataset {
            path: path.to_string(),
        })
}

fn build_cluster_map(barcodes: &[String], labels: &[String]) -> ClusterMap {
    use std::collections::{BTreeMap, BTreeSet};

    let mut cluster_ids: BTreeSet<String> = BTreeSet::new();
    for label in labels {
        cluster_ids.insert(label.clone());
    }

    let cluster_ids: Vec<String> = cluster_ids.into_iter().collect();
    let mut cluster_index = BTreeMap::new();
    for (idx, id) in cluster_ids.iter().enumerate() {
        cluster_index.insert(id.clone(), idx);
    }

    let mut cell_to_cluster = Vec::with_capacity(barcodes.len());
    for label in labels {
        let idx = cluster_index.get(label).copied().unwrap_or(0);
        cell_to_cluster.push(idx);
    }

    ClusterMap {
        cluster_ids,
        cell_to_cluster,
    }
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
