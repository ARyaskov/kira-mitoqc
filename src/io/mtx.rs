//! Matrix Market (MTX) loader.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use sprs::CsMat;
use tracing::{debug, info};

use crate::input::InputError;
use crate::io::feature_detect::detect_gene_symbol_column;

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

/// Load features.tsv and barcodes.tsv without parsing the matrix.
pub fn load_mtx_metadata(
    path: &Path,
    gene_symbol_col: Option<usize>,
) -> Result<(Vec<String>, Vec<String>), InputError> {
    let input_dir = resolve_input_dir(path)?.to_path_buf();
    let prefix = detect_prefix(&input_dir)?;
    let feature_path = candidate_path(&input_dir, prefix.as_deref(), "features.tsv");
    let genes_path = candidate_path(&input_dir, prefix.as_deref(), "genes.tsv");
    let barcodes_path = candidate_path(&input_dir, prefix.as_deref(), "barcodes.tsv");
    let (feature_path, feature_kind) = if exists_plain_or_gz(&feature_path) {
        (feature_path, FeatureFileKind::FeaturesTsv)
    } else if exists_plain_or_gz(&genes_path) {
        (genes_path, FeatureFileKind::GenesTsv)
    } else {
        return Err(InputError::MissingFile { path: feature_path });
    };

    let (features_reader, features_opened) = open_maybe_gz(&feature_path)?;
    debug!(path = ?features_opened, "opened features");
    let features = match feature_kind {
        FeatureFileKind::FeaturesTsv => {
            read_features_with_detection(&feature_path, features_reader, gene_symbol_col)?
        }
        FeatureFileKind::GenesTsv => {
            if gene_symbol_col.is_some() {
                return Err(InputError::LegacyGenesTsvInvalid {
                    reason: "gene symbol column override is not supported for genes.tsv"
                        .to_string(),
                });
            }
            read_genes_tsv_symbols(features_reader, &features_opened)?
        }
    };

    let (barcodes_reader, barcodes_opened) = open_maybe_gz(&barcodes_path)?;
    debug!(path = ?barcodes_opened, "opened barcodes");
    let barcodes = read_lines_from_reader(barcodes_reader, &barcodes_opened)?;

    Ok((features, barcodes))
}

/// Load Matrix Market files from a directory.
pub fn load_mtx_dir(path: &Path, gene_symbol_col: Option<usize>) -> Result<MtxInput, InputError> {
    info!(path = ?path, "Loading MTX directory");
    let discovery = discover_dataset_files(path)?;
    info!(
        matrix = ?discovery.matrix_path,
        features = ?discovery.feature_path,
        barcodes = ?discovery.barcodes_path,
        prefix = discovery.prefix.as_deref().unwrap_or(""),
        "Resolved MTX input files"
    );

    let (features_reader, features_opened) = open_maybe_gz(&discovery.feature_path)?;
    debug!(path = ?features_opened, "opened features");
    let features = match discovery.feature_kind {
        FeatureFileKind::FeaturesTsv => {
            read_features_with_detection(&discovery.feature_path, features_reader, gene_symbol_col)?
        }
        FeatureFileKind::GenesTsv => {
            if gene_symbol_col.is_some() {
                return Err(InputError::LegacyGenesTsvInvalid {
                    reason: "gene symbol column override is not supported for genes.tsv"
                        .to_string(),
                });
            }
            read_genes_tsv_symbols(features_reader, &features_opened)?
        }
    };

    let (barcodes_reader, barcodes_opened) = open_maybe_gz(&discovery.barcodes_path)?;
    debug!(path = ?barcodes_opened, "opened barcodes");
    let barcodes = read_lines_from_reader(barcodes_reader, &barcodes_opened)?;

    let (_, matrix_opened) = open_maybe_gz(&discovery.matrix_path)?;
    debug!(path = ?matrix_opened, "opened matrix");
    info!(path = ?discovery.matrix_path, "Parsing matrix market (streaming)");
    let matrix = read_matrix_market(&discovery.matrix_path)?;
    info!(
        rows = matrix.rows(),
        cols = matrix.cols(),
        nnz = matrix.nnz(),
        "Loaded matrix market into CSC"
    );

    if matrix.rows() != features.len() || matrix.cols() != barcodes.len() {
        return Err(InputError::DimensionMismatch {
            rows: matrix.rows(),
            cols: matrix.cols(),
            features: features.len(),
            barcodes: barcodes.len(),
        });
    }

    Ok(MtxInput {
        matrix,
        features,
        barcodes,
    })
}

/// Discover MTX/feature/barcode files in an input directory, including optional dataset prefix.
pub fn discover_dataset_files(path: &Path) -> Result<MtxDiscovery, InputError> {
    let input_dir = resolve_input_dir(path)?.to_path_buf();
    let prefix = detect_prefix(&input_dir)?;

    let matrix_path = candidate_path(&input_dir, prefix.as_deref(), "matrix.mtx");
    let features_path = candidate_path(&input_dir, prefix.as_deref(), "features.tsv");
    let genes_path = candidate_path(&input_dir, prefix.as_deref(), "genes.tsv");
    let barcodes_path = candidate_path(&input_dir, prefix.as_deref(), "barcodes.tsv");

    if !exists_plain_or_gz(&matrix_path) {
        return Err(InputError::MissingFile { path: matrix_path });
    }
    if !exists_plain_or_gz(&barcodes_path) {
        return Err(InputError::MissingFile {
            path: barcodes_path,
        });
    }

    let (feature_path, feature_kind) = if exists_plain_or_gz(&features_path) {
        info!("Detected 10x v3+ format: using features.tsv");
        (features_path, FeatureFileKind::FeaturesTsv)
    } else if exists_plain_or_gz(&genes_path) {
        info!("Detected legacy 10x v2 format: using genes.tsv (column 2 as gene symbol)");
        (genes_path, FeatureFileKind::GenesTsv)
    } else {
        return Err(InputError::MissingFile {
            path: features_path,
        });
    };

    Ok(MtxDiscovery {
        input_dir,
        prefix,
        matrix_path,
        feature_path,
        feature_kind,
        barcodes_path,
    })
}

fn candidate_path(input_dir: &Path, prefix: Option<&str>, name: &str) -> PathBuf {
    match prefix {
        Some(p) => input_dir.join(format!("{p}_{name}")),
        None => input_dir.join(name),
    }
}

/// Resolve shared pipeline cache filename from optional dataset prefix.
pub fn resolve_shared_cache_filename(prefix: Option<&str>) -> String {
    match prefix {
        Some(p) if !p.is_empty() => format!("{p}.kira-organelle.bin"),
        _ => "kira-organelle.bin".to_string(),
    }
}

fn open_maybe_gz(path: &Path) -> Result<(Box<dyn Read>, PathBuf), InputError> {
    if path.exists() {
        let file = File::open(path).map_err(|source| InputError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        return Ok((Box::new(file), path.to_path_buf()));
    }

    let gz_path = gz_path(path);
    if gz_path.exists() {
        let file = File::open(&gz_path).map_err(|source| InputError::Io {
            path: gz_path.clone(),
            source,
        })?;
        let decoder = flate2::read::GzDecoder::new(file);
        return Ok((Box::new(decoder), gz_path));
    }

    Err(InputError::MissingFile {
        path: path.to_path_buf(),
    })
}

fn resolve_input_dir(path: &Path) -> Result<&Path, InputError> {
    if path.is_dir() {
        if contains_dataset_files(path)? {
            return Ok(path);
        }
        return Err(InputError::InvalidInputPath {
            path: path.to_path_buf(),
        });
    }
    let parent = path.parent().ok_or_else(|| InputError::InvalidInputPath {
        path: path.to_path_buf(),
    })?;
    Ok(parent)
}

fn contains_dataset_files(path: &Path) -> Result<bool, InputError> {
    let entries = std::fs::read_dir(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    for entry in entries {
        let entry = entry.map_err(|source| InputError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        let name = entry.file_name();
        let Some(name) = name.to_str() else {
            continue;
        };
        if name == "matrix.mtx"
            || name == "matrix.mtx.gz"
            || name == "features.tsv"
            || name == "features.tsv.gz"
            || name == "genes.tsv"
            || name == "genes.tsv.gz"
            || name == "barcodes.tsv"
            || name == "barcodes.tsv.gz"
            || extract_prefix(name).is_some()
        {
            return Ok(true);
        }
    }
    Ok(false)
}

fn detect_prefix(input_dir: &Path) -> Result<Option<String>, InputError> {
    let mut prefixes = std::collections::BTreeSet::new();
    let entries = std::fs::read_dir(input_dir).map_err(|source| InputError::Io {
        path: input_dir.to_path_buf(),
        source,
    })?;
    for entry in entries {
        let entry = entry.map_err(|source| InputError::Io {
            path: input_dir.to_path_buf(),
            source,
        })?;
        let name_os = entry.file_name();
        let Some(name) = name_os.to_str() else {
            continue;
        };
        if let Some(prefix) = extract_prefix(name) {
            prefixes.insert(prefix.to_string());
        }
    }
    if prefixes.len() > 1 {
        return Err(InputError::InvalidInputPath {
            path: input_dir.to_path_buf(),
        });
    }
    Ok(prefixes.into_iter().next())
}

fn extract_prefix(name: &str) -> Option<&str> {
    let suffixes = [
        "_matrix.mtx",
        "_matrix.mtx.gz",
        "_features.tsv",
        "_features.tsv.gz",
        "_barcodes.tsv",
        "_barcodes.tsv.gz",
        "_genes.tsv",
        "_genes.tsv.gz",
    ];
    for suffix in suffixes {
        if let Some(prefix) = name.strip_suffix(suffix)
            && !prefix.is_empty()
        {
            return Some(prefix);
        }
    }
    None
}

fn exists_plain_or_gz(path: &Path) -> bool {
    path.exists() || gz_path(path).exists()
}

fn gz_path(path: &Path) -> PathBuf {
    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
        path.with_extension(format!("{ext}.gz"))
    } else {
        path.with_extension("gz")
    }
}

fn read_features_with_detection(
    path: &Path,
    reader: Box<dyn Read>,
    gene_symbol_col: Option<usize>,
) -> Result<Vec<String>, InputError> {
    let mut detect_reader = reader;
    let column = if let Some(col) = gene_symbol_col {
        col
    } else {
        detect_gene_symbol_column(&mut detect_reader, 100)?
    };

    let (reader, _) = open_maybe_gz(path)?;
    read_column_from_reader(reader, path, column)
}

fn read_genes_tsv_symbols(reader: Box<dyn Read>, path: &Path) -> Result<Vec<String>, InputError> {
    let reader = BufReader::new(reader);
    let mut values = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.map_err(|source| InputError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            return Err(InputError::LegacyGenesTsvInvalid {
                reason: format!("genes.tsv line {} has fewer than 2 columns", line_no + 1),
            });
        }
        values.push(parts[1].trim().to_string());
    }

    if values.is_empty() {
        return Err(InputError::LegacyGenesTsvInvalid {
            reason: "genes.tsv is empty".to_string(),
        });
    }

    Ok(values)
}

fn read_column_from_reader(
    reader: Box<dyn Read>,
    path: &Path,
    column: usize,
) -> Result<Vec<String>, InputError> {
    let reader = BufReader::new(reader);
    let mut values = Vec::new();
    let mut max_cols = 0usize;
    for line in reader.lines() {
        let line = line.map_err(|source| InputError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() > max_cols {
            max_cols = parts.len();
        }
        let value = parts.get(column).map(|v| v.trim()).unwrap_or("");
        values.push(value.to_string());
    }

    if max_cols > 0 && column >= max_cols {
        return Err(InputError::InvalidGeneSymbolColumn {
            requested: column + 1,
            available: max_cols,
        });
    }

    Ok(values)
}

fn read_lines_from_reader<R: Read>(reader: R, path: &Path) -> Result<Vec<String>, InputError> {
    let reader = BufReader::new(reader);
    let mut values = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|source| InputError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        if line.trim().is_empty() {
            continue;
        }
        values.push(line.trim().to_string());
    }
    Ok(values)
}

fn read_matrix_market(path: &Path) -> Result<CsMat<f32>, InputError> {
    // Fast path: parse as real-valued matrix directly.
    let (reader, opened) = open_maybe_gz(path)?;
    let mut reader = BufReader::new(reader);
    match sprs::io::read_matrix_market_from_bufread::<f32, usize, _>(&mut reader) {
        Ok(tri) => return Ok(tri.to_csc()),
        Err(err) => {
            let msg = err.to_string();
            if !msg.contains("integer file into real matrix") {
                return Err(InputError::MatrixParse {
                    path: opened,
                    message: msg,
                });
            }
        }
    }

    // Fallback for integer MM files without buffering the whole file.
    let (reader, opened) = open_maybe_gz(path)?;
    let mut reader = BufReader::new(reader);
    let tri = sprs::io::read_matrix_market_from_bufread::<i64, usize, _>(&mut reader).map_err(
        |source| InputError::MatrixParse {
            path: opened.clone(),
            message: source.to_string(),
        },
    )?;
    let csc = tri.to_csc();
    let data: Vec<f32> = csc.data().iter().map(|v| *v as f32).collect();
    let indptr = csc.indptr().raw_storage().to_vec();
    let indices = csc.indices().to_vec();
    Ok(CsMat::new_csc(csc.shape(), indptr, indices, data))
}
