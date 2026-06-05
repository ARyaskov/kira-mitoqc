//! Runtime CPUID dispatch — detection cached via `OnceLock`.

use std::sync::OnceLock;

use crate::cache::ExpressionSoAView;
use crate::compute::{GeneOffsets, PrimitiveSignals};
use crate::input::ResolvedGeneSets;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimdKind {
    Avx2,
    Neon,
    Scalar,
}

fn detect_simd() -> SimdKind {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2")
            && std::arch::is_x86_feature_detected!("fma")
        {
            return SimdKind::Avx2;
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            return SimdKind::Neon;
        }
    }
    SimdKind::Scalar
}

pub fn active_kind() -> SimdKind {
    static KIND: OnceLock<SimdKind> = OnceLock::new();
    *KIND.get_or_init(detect_simd)
}

pub fn compute_primitives(
    soa: &ExpressionSoAView<'_>,
    resolved: &ResolvedGeneSets,
) -> PrimitiveSignals {
    let offsets = GeneOffsets::from_resolved(resolved);

    match active_kind() {
        #[cfg(target_arch = "x86_64")]
        SimdKind::Avx2 => unsafe { crate::compute::avx2::compute_primitives_avx2(soa, &offsets) },
        #[cfg(target_arch = "aarch64")]
        SimdKind::Neon => crate::compute::neon::compute_primitives_neon(soa, &offsets),
        _ => crate::compute::scalar::compute_primitives_scalar(soa, &offsets),
    }
}
