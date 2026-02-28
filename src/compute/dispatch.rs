//! Runtime dispatch for compute kernels.

use crate::cache::ExpressionSoAView;
use crate::compute::{GeneOffsets, PrimitiveSignals};
use crate::input::ResolvedGeneSets;

pub fn compute_primitives(
    soa: &ExpressionSoAView,
    resolved: &ResolvedGeneSets,
) -> PrimitiveSignals {
    let offsets = GeneOffsets::from_resolved(resolved);

    if cfg!(all(target_arch = "x86_64", target_feature = "avx2")) {
        #[cfg(target_feature = "avx2")]
        {
            return crate::compute::avx2::compute_primitives_avx2(soa, &offsets);
        }
    }

    if cfg!(target_arch = "aarch64") {
        #[cfg(target_arch = "aarch64")]
        {
            return crate::compute::neon::compute_primitives_neon(soa, &offsets);
        }
    }

    crate::compute::scalar::compute_primitives_scalar(soa, &offsets)
}
