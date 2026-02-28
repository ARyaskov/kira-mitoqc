use std::path::{Path, PathBuf};

use sprs::CsMat;
use thiserror::Error;

use crate::input::InputError;
use crate::io::mtx::{MtxInput, load_mtx_dir};

#[derive(Debug, Error)]
pub enum OrganelleCacheError {
    #[error("I/O error reading {path:?}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("cache format error in {path:?}: {message}")]
    Format { path: PathBuf, message: String },
    #[error("input parsing error: {0}")]
    Input(#[from] InputError),
}

#[derive(Debug)]
pub struct OrganelleCacheView {
    pub n_genes: usize,
    pub n_cells: usize,
    pub nnz: usize,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
    inner: kira_shared_sc_cache::SharedCacheMmap,
}

impl OrganelleCacheView {
    pub fn col_ptr(&self) -> &[u64] {
        self.inner.col_ptr()
    }

    pub fn row_idx(&self) -> &[u32] {
        self.inner.row_idx()
    }

    pub fn values_u32(&self) -> &[u32] {
        self.inner.values_u32()
    }

    pub fn to_mtx_input(&self) -> MtxInput {
        let indptr: Vec<usize> = self.col_ptr().iter().map(|v| *v as usize).collect();
        let indices: Vec<usize> = self.row_idx().iter().map(|v| *v as usize).collect();
        let data: Vec<f32> = self.values_u32().iter().map(|v| *v as f32).collect();
        let matrix = CsMat::new_csc((self.n_genes, self.n_cells), indptr, indices, data);
        MtxInput {
            matrix,
            features: self.genes.clone(),
            barcodes: self.barcodes.clone(),
        }
    }
}

pub fn write_organelle_bin_from_mtx(
    out_path: &Path,
    input_path: &Path,
    gene_symbol_col: Option<usize>,
) -> Result<(), OrganelleCacheError> {
    let mtx = load_mtx_dir(input_path, gene_symbol_col)?;
    write_organelle_bin(out_path, &mtx)
}

pub fn write_organelle_bin(path: &Path, input: &MtxInput) -> Result<(), OrganelleCacheError> {
    if input.features.len() != input.matrix.rows() {
        return Err(OrganelleCacheError::Format {
            path: path.to_path_buf(),
            message: "feature count does not match matrix rows".to_string(),
        });
    }
    if input.barcodes.len() != input.matrix.cols() {
        return Err(OrganelleCacheError::Format {
            path: path.to_path_buf(),
            message: "barcode count does not match matrix cols".to_string(),
        });
    }

    let mut row_idx = Vec::with_capacity(input.matrix.nnz());
    let mut values_u32 = Vec::with_capacity(input.matrix.nnz());
    for (col_idx, col) in input.matrix.outer_iterator().enumerate() {
        let mut prev_row: Option<usize> = None;
        for (row, value) in col.iter() {
            if let Some(prev) = prev_row
                && row <= prev
            {
                return Err(OrganelleCacheError::Format {
                    path: path.to_path_buf(),
                    message: format!("row_idx not strictly increasing in column {col_idx}"),
                });
            }
            prev_row = Some(row);
            if !value.is_finite()
                || *value < 0.0
                || value.fract() != 0.0
                || *value > u32::MAX as f32
            {
                return Err(OrganelleCacheError::Format {
                    path: path.to_path_buf(),
                    message: "values_u32 must be finite non-negative integer counts".to_string(),
                });
            }
            row_idx.push(row as u32);
            values_u32.push(*value as u32);
        }
    }

    let col_ptr: Vec<u64> = input
        .matrix
        .indptr()
        .raw_storage()
        .iter()
        .map(|v| *v as u64)
        .collect();

    let write_input = kira_shared_sc_cache::SharedCacheWriteInput {
        genes: &input.features,
        barcodes: &input.barcodes,
        col_ptr: &col_ptr,
        row_idx: &row_idx,
        values_u32: &values_u32,
    };
    kira_shared_sc_cache::write_shared_cache(path, &write_input).map_err(map_err)
}

pub fn mmap_organelle_bin(path: &Path) -> Result<OrganelleCacheView, OrganelleCacheError> {
    let inner = kira_shared_sc_cache::mmap_shared_cache(path).map_err(map_err)?;
    Ok(OrganelleCacheView {
        n_genes: inner.n_genes,
        n_cells: inner.n_cells,
        nnz: inner.nnz,
        genes: inner.genes.clone(),
        barcodes: inner.barcodes.clone(),
        inner,
    })
}

fn map_err(err: kira_shared_sc_cache::SharedCacheError) -> OrganelleCacheError {
    match err {
        kira_shared_sc_cache::SharedCacheError::Io { path, source } => {
            OrganelleCacheError::Io { path, source }
        }
        kira_shared_sc_cache::SharedCacheError::Format { path, message } => {
            OrganelleCacheError::Format { path, message }
        }
    }
}
