use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::data::ExpressionSoA;

pub use kira_shared_sc_cache::ExprCacheMode;

#[derive(Debug, Error)]
pub enum CacheError {
    #[error("I/O error reading {path:?}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("invalid magic in {path:?}")]
    InvalidMagic { path: PathBuf },
    #[error("unsupported version {version} in {path:?}")]
    UnsupportedVersion { path: PathBuf, version: u32 },
    #[error("truncated file {path:?}")]
    Truncated { path: PathBuf },
    #[error("size mismatch in {path:?}: expected {expected} bytes, got {actual}")]
    SizeMismatch {
        path: PathBuf,
        expected: usize,
        actual: usize,
    },
}

#[derive(Debug)]
pub struct ExpressionSoAView<'a> {
    pub values: &'a [f32],
    pub genes: usize,
    pub samples: usize,
    pub mode: ExprCacheMode,
    #[allow(dead_code)]
    inner: kira_shared_sc_cache::ExprBinMmap,
}

impl<'a> ExpressionSoAView<'a> {
    pub fn get(&self, gene: usize, sample: usize) -> f32 {
        let idx = gene * self.samples + sample;
        self.values[idx]
    }
}

pub fn write_expr_bin(path: &Path, soa: &ExpressionSoA) -> Result<(), CacheError> {
    write_expr_bin_with_mode(path, soa, ExprCacheMode::Unknown)
}

pub fn write_expr_bin_with_mode(
    path: &Path,
    soa: &ExpressionSoA,
    mode: ExprCacheMode,
) -> Result<(), CacheError> {
    kira_shared_sc_cache::write_expr_bin_with_mode(path, soa.genes, soa.samples, &soa.values, mode)
        .map_err(map_err)
}

pub fn mmap_expr_bin(path: &Path) -> Result<ExpressionSoAView<'static>, CacheError> {
    let inner = kira_shared_sc_cache::mmap_expr_bin(path).map_err(map_err)?;
    let genes = inner.genes;
    let samples = inner.samples;
    let mode = inner.mode;
    let values = inner.values();
    let values_static: &'static [f32] = unsafe { std::mem::transmute(values) };
    Ok(ExpressionSoAView {
        values: values_static,
        genes,
        samples,
        mode,
        inner,
    })
}

fn map_err(err: kira_shared_sc_cache::ExprBinError) -> CacheError {
    match err {
        kira_shared_sc_cache::ExprBinError::Io { path, source } => CacheError::Io { path, source },
        kira_shared_sc_cache::ExprBinError::InvalidMagic { path } => {
            CacheError::InvalidMagic { path }
        }
        kira_shared_sc_cache::ExprBinError::UnsupportedVersion { path, version } => {
            CacheError::UnsupportedVersion { path, version }
        }
        kira_shared_sc_cache::ExprBinError::Truncated { path } => CacheError::Truncated { path },
        kira_shared_sc_cache::ExprBinError::SizeMismatch {
            path,
            expected,
            actual,
        } => CacheError::SizeMismatch {
            path,
            expected,
            actual,
        },
    }
}
