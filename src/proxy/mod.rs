//! Proxy computation APIs.

use thiserror::Error;

use crate::compute::PrimitiveSignals;
use crate::config::refs::RefsV1;
use crate::core::types::{ProxyKey, ProxyScores};
use crate::input::ResolvedGeneSets;

pub mod v1;
pub mod v2;

/// Proxy computation errors.
#[derive(Debug, Error)]
pub enum ProxyError {
    #[error("insufficient genes for {set}: found {found}, required {required}")]
    InsufficientGenes {
        set: &'static str,
        found: usize,
        required: usize,
    },

    #[error("invalid numeric result for {proxy} at sample {sample}")]
    InvalidNumeric { proxy: &'static str, sample: usize },
}

/// Compute v1 proxy scores from primitive signals.
pub fn compute_proxies_v1(
    primitives: &PrimitiveSignals,
    resolved: &ResolvedGeneSets,
    refs: &RefsV1,
) -> Result<ProxyScores, ProxyError> {
    v1::compute_proxies_v1(primitives, resolved, refs)
}

pub use v2::{OptionalOmicsInputs, ProxyKeyV2, ProxyScoresV2, compute_proxies_v2};

pub(crate) fn validate_no_nan(proxy: ProxyKey, values: &[f32]) -> Result<(), ProxyError> {
    for (idx, value) in values.iter().enumerate() {
        if value.is_nan() {
            return Err(ProxyError::InvalidNumeric {
                proxy: proxy.as_str(),
                sample: idx,
            });
        }
    }
    Ok(())
}
