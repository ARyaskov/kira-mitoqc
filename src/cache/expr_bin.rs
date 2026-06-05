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

/// Mmap-backed expression cache. `values()` borrows from `inner`; the
/// `'a` lifetime parameter is kept for backwards source compatibility.
#[derive(Debug)]
pub struct ExpressionSoAView<'a> {
    pub genes: usize,
    pub samples: usize,
    pub mode: ExprCacheMode,
    inner: kira_shared_sc_cache::ExprBinMmap,
    _marker: std::marker::PhantomData<&'a ()>,
}

impl<'a> ExpressionSoAView<'a> {
    #[inline]
    pub fn values(&self) -> &[f32] {
        self.inner.values()
    }

    #[inline]
    pub fn get(&self, gene: usize, sample: usize) -> f32 {
        let idx = gene * self.samples + sample;
        self.values()[idx]
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

    let expected = genes.saturating_mul(samples);
    if inner.values().len() != expected {
        return Err(CacheError::SizeMismatch {
            path: path.to_path_buf(),
            expected: expected.saturating_mul(std::mem::size_of::<f32>()),
            actual: inner.values().len() * std::mem::size_of::<f32>(),
        });
    }

    Ok(ExpressionSoAView {
        genes,
        samples,
        mode,
        inner,
        _marker: std::marker::PhantomData,
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
