//! Binary cache for prepared expression SoA.
//!
//! Format (little-endian):
//!
//! ```text
//! struct ExprBinHeader {
//!     magic: [u8; 8];        // "KIRAMTX\0"
//!     version: u32;          // 1
//!     genes: u32;            // number of genes
//!     samples: u32;          // number of samples
//!     flags: u32;            // cache metadata (mode in low 8 bits)
//! }
//! ```
//!
//! Immediately followed by `genes * samples` f32 values in row-major
//! `[gene][sample]` order with no padding.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use memmap2::Mmap;
use thiserror::Error;

use crate::data::ExpressionSoA;

const MAGIC: &[u8; 8] = b"KIRAMTX\0";
const VERSION: u32 = 1;
const HEADER_LEN: usize = 8 + 4 + 4 + 4 + 4;
const MODE_MASK: u32 = 0xFF;

/// Aggregation mode encoded in expression cache metadata.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ExprCacheMode {
    Unknown,
    Sample,
    Cluster,
    Cell,
}

impl ExprCacheMode {
    fn to_flags(self) -> u32 {
        match self {
            Self::Unknown => 0,
            Self::Sample => 1,
            Self::Cluster => 2,
            Self::Cell => 3,
        }
    }

    fn from_flags(flags: u32) -> Self {
        match flags & MODE_MASK {
            1 => Self::Sample,
            2 => Self::Cluster,
            3 => Self::Cell,
            _ => Self::Unknown,
        }
    }
}

/// Errors related to the expression cache format.
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

/// A zero-copy view into a memory-mapped expression cache.
#[derive(Debug)]
pub struct ExpressionSoAView<'a> {
    pub values: &'a [f32],
    pub genes: usize,
    pub samples: usize,
    pub mode: ExprCacheMode,
    #[allow(dead_code)]
    mmap: Mmap,
}

impl<'a> ExpressionSoAView<'a> {
    /// Get a value at (gene, sample).
    pub fn get(&self, gene: usize, sample: usize) -> f32 {
        let idx = gene * self.samples + sample;
        self.values[idx]
    }
}

/// Write an ExpressionSoA to the binary cache format.
pub fn write_expr_bin(path: &Path, soa: &ExpressionSoA) -> Result<(), CacheError> {
    write_expr_bin_with_mode(path, soa, ExprCacheMode::Unknown)
}

/// Write an ExpressionSoA with aggregation mode metadata.
pub fn write_expr_bin_with_mode(
    path: &Path,
    soa: &ExpressionSoA,
    mode: ExprCacheMode,
) -> Result<(), CacheError> {
    let mut file = File::create(path).map_err(|source| CacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    let genes = soa.genes as u32;
    let samples = soa.samples as u32;
    let values_len =
        (soa.genes)
            .checked_mul(soa.samples)
            .ok_or_else(|| CacheError::SizeMismatch {
                path: path.to_path_buf(),
                expected: usize::MAX,
                actual: soa.values.len(),
            })?;

    if values_len != soa.values.len() {
        return Err(CacheError::SizeMismatch {
            path: path.to_path_buf(),
            expected: values_len * 4 + HEADER_LEN,
            actual: soa.values.len() * 4 + HEADER_LEN,
        });
    }

    file.write_all(MAGIC).map_err(|source| CacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    file.write_all(&VERSION.to_le_bytes())
        .map_err(|source| CacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    file.write_all(&genes.to_le_bytes())
        .map_err(|source| CacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    file.write_all(&samples.to_le_bytes())
        .map_err(|source| CacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;
    file.write_all(&mode.to_flags().to_le_bytes())
        .map_err(|source| CacheError::Io {
            path: path.to_path_buf(),
            source,
        })?;

    let bytes: &[u8] = bytemuck_cast_slice(&soa.values);
    file.write_all(bytes).map_err(|source| CacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    file.sync_all().map_err(|source| CacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    Ok(())
}

/// Memory-map an expression cache file.
pub fn mmap_expr_bin(path: &Path) -> Result<ExpressionSoAView<'static>, CacheError> {
    let file = File::open(path).map_err(|source| CacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let metadata = file.metadata().map_err(|source| CacheError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.len() < HEADER_LEN as u64 {
        return Err(CacheError::Truncated {
            path: path.to_path_buf(),
        });
    }

    let mmap = unsafe {
        Mmap::map(&file).map_err(|source| CacheError::Io {
            path: path.to_path_buf(),
            source,
        })?
    };

    let header = &mmap[..HEADER_LEN];
    if &header[..8] != MAGIC {
        return Err(CacheError::InvalidMagic {
            path: path.to_path_buf(),
        });
    }

    let version = read_u32_le(header, 8);
    if version != VERSION {
        return Err(CacheError::UnsupportedVersion {
            path: path.to_path_buf(),
            version,
        });
    }

    let genes = read_u32_le(header, 12) as usize;
    let samples = read_u32_le(header, 16) as usize;
    let mode = ExprCacheMode::from_flags(read_u32_le(header, 20));

    let values_len = genes
        .checked_mul(samples)
        .ok_or_else(|| CacheError::SizeMismatch {
            path: path.to_path_buf(),
            expected: usize::MAX,
            actual: metadata.len() as usize,
        })?;

    let expected_bytes = HEADER_LEN + values_len * 4;
    let actual_bytes = metadata.len() as usize;
    if expected_bytes != actual_bytes {
        return Err(CacheError::SizeMismatch {
            path: path.to_path_buf(),
            expected: expected_bytes,
            actual: actual_bytes,
        });
    }

    let values_bytes = &mmap[HEADER_LEN..];
    let values =
        unsafe { std::slice::from_raw_parts(values_bytes.as_ptr() as *const f32, values_len) };
    let values_static: &'static [f32] = unsafe { std::mem::transmute(values) };

    Ok(ExpressionSoAView {
        values: values_static,
        genes,
        samples,
        mode,
        mmap,
    })
}

fn bytemuck_cast_slice(values: &[f32]) -> &[u8] {
    let byte_len = values.len() * std::mem::size_of::<f32>();
    unsafe { std::slice::from_raw_parts(values.as_ptr() as *const u8, byte_len) }
}

#[inline]
fn read_u32_le(buf: &[u8], offset: usize) -> u32 {
    let bytes = buf[offset..offset + 4]
        .as_array::<4>()
        .expect("validated header width");
    u32::from_le_bytes(*bytes)
}
