//! Input/output modules for concrete formats.

pub mod feature_detect;
pub mod mtx;

#[cfg(feature = "h5ad")]
pub mod h5ad;
