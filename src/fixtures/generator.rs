//! Synthetic primitive signal generator.

use crate::compute::PrimitiveSignals;
use crate::util::numeric::clamp01;

/// Synthetic factor model for fixtures.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SyntheticFactors {
    pub bioenergetic_damage: f32,
    pub ros_pressure: f32,
    pub dynamics_instability: f32,
    pub biogenesis_failure: f32,
    pub compensation: f32,
}

/// Generate a linear series of factors.
pub fn factors_linear(n: usize, f: impl Fn(f32) -> SyntheticFactors) -> Vec<SyntheticFactors> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![f(0.0)];
    }
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f32 / (n as f32 - 1.0);
        out.push(f(t));
    }
    out
}

/// Generate PrimitiveSignals (single sample) from synthetic factors.
pub fn generate_primitives_v1(factors: &SyntheticFactors) -> PrimitiveSignals {
    let b = factors.bioenergetic_damage;
    let r = factors.ros_pressure;
    let d = factors.dynamics_instability;
    let g = factors.biogenesis_failure;
    let c = factors.compensation;

    let base_c = 0.35_f32;

    let c_i = clamp01(base_c - 0.20 * b + 0.10 * c);
    let c_iii = clamp01(base_c + 0.05 * b);
    let c_iv = clamp01(base_c + 0.10 * b);
    let c_v = clamp01(base_c - 0.05 * b);

    let mtdna_mean = clamp01(0.40 - 0.25 * b + 0.05 * c);
    let nuclear_mean = clamp01(0.40 - 0.10 * b - 0.20 * g + 0.05 * c);

    let ros_mean = clamp01(0.20 + 0.90 * r - 0.20 * c);
    let mitophagy_mean = clamp01(0.15 + 0.90 * d + 0.30 * r - 0.20 * c);

    let fusion_mean = clamp01(0.35 - 0.25 * d + 0.10 * c);
    let fission_mean = clamp01(0.25 + 0.55 * d - 0.10 * c);

    let biogenesis_mean = clamp01(0.45 - 0.80 * g + 0.20 * c);

    let atp_mt = clamp01(0.35 - 0.30 * b + 0.10 * c);
    let atp_nu = clamp01(0.30 - 0.05 * b - 0.10 * g + 0.05 * c);

    let stoich_variance = variance4(c_i, c_iii, c_iv, c_v);

    PrimitiveSignals {
        mtdna_mean: vec![mtdna_mean],
        nuclear_mean: vec![nuclear_mean],
        c_i: vec![c_i],
        c_iii: vec![c_iii],
        c_iv: vec![c_iv],
        c_v: vec![c_v],
        ros_mean: vec![ros_mean],
        mitophagy_mean: vec![mitophagy_mean],
        fusion_mean: vec![fusion_mean],
        fission_mean: vec![fission_mean],
        biogenesis_mean: vec![biogenesis_mean],
        atp_mt: vec![atp_mt],
        atp_nu: vec![atp_nu],
        stoich_variance: vec![stoich_variance],
    }
}

fn variance4(a: f32, b: f32, c: f32, d: f32) -> f32 {
    let mean = (a + b + c + d) * 0.25;
    let da = a - mean;
    let db = b - mean;
    let dc = c - mean;
    let dd = d - mean;
    (da * da + db * db + dc * dc + dd * dd) * 0.25
}
