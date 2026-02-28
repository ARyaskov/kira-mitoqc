//! AVX2 implementation for primitive signals.

#![cfg(target_feature = "avx2")]

use wide::f32x8;

use crate::cache::ExpressionSoAView;
use crate::compute::{GeneOffsets, PrimitiveSignals};

pub fn compute_primitives_avx2(soa: &ExpressionSoAView, offsets: &GeneOffsets) -> PrimitiveSignals {
    let samples = soa.samples;

    let mut mtdna_mean = vec![0.0; samples];
    let mut nuclear_mean = vec![0.0; samples];
    let mut c_i = vec![0.0; samples];
    let mut c_iii = vec![0.0; samples];
    let mut c_iv = vec![0.0; samples];
    let mut c_v = vec![0.0; samples];
    let mut ros_mean = vec![0.0; samples];
    let mut mitophagy_mean = vec![0.0; samples];
    let mut fusion_mean = vec![0.0; samples];
    let mut fission_mean = vec![0.0; samples];
    let mut biogenesis_mean = vec![0.0; samples];
    let mut atp_mt = vec![0.0; samples];
    let mut atp_nu = vec![0.0; samples];
    let mut stoich_variance = vec![0.0; samples];

    mean_over_offsets_avx2(soa, &offsets.mtdna_all, &mut mtdna_mean);
    mean_over_offsets_avx2(soa, &offsets.nuclear_oxphos, &mut nuclear_mean);

    mean_over_offsets_avx2(soa, &offsets.complex_i, &mut c_i);
    mean_over_offsets_avx2(soa, &offsets.complex_iii, &mut c_iii);
    mean_over_offsets_avx2(soa, &offsets.complex_iv, &mut c_iv);
    mean_over_offsets_avx2(soa, &offsets.complex_v, &mut c_v);

    mean_over_offsets_avx2(soa, &offsets.ros, &mut ros_mean);
    mean_over_offsets_avx2(soa, &offsets.mitophagy, &mut mitophagy_mean);
    mean_over_offsets_avx2(soa, &offsets.fusion, &mut fusion_mean);
    mean_over_offsets_avx2(soa, &offsets.fission, &mut fission_mean);
    mean_over_offsets_avx2(soa, &offsets.biogenesis, &mut biogenesis_mean);

    mean_over_offsets_avx2(soa, &offsets.atp_mt, &mut atp_mt);
    copy_single_gene_avx2(soa, offsets.atp_nu, &mut atp_nu);

    variance4_avx2(&c_i, &c_iii, &c_iv, &c_v, &mut stoich_variance);

    PrimitiveSignals {
        mtdna_mean,
        nuclear_mean,
        c_i,
        c_iii,
        c_iv,
        c_v,
        ros_mean,
        mitophagy_mean,
        fusion_mean,
        fission_mean,
        biogenesis_mean,
        atp_mt,
        atp_nu,
        stoich_variance,
    }
}

fn mean_over_offsets_avx2(soa: &ExpressionSoAView, offsets: &[usize], out: &mut [f32]) {
    let samples = soa.samples;
    if offsets.is_empty() {
        for value in out.iter_mut() {
            *value = 0.0;
        }
        return;
    }

    let count = offsets.len() as f32;
    let denom = f32x8::new([count; 8]);
    let chunks = samples / 8 * 8;

    for s in (0..chunks).step_by(8) {
        let mut acc = f32x8::new([0.0; 8]);
        for &g in offsets {
            let base = g * samples + s;
            let v = load_f32x8(&soa.values, base);
            acc = acc + v;
        }
        let mean = acc / denom;
        store_f32x8(mean, out, s);
    }

    for s in chunks..samples {
        let mut sum = 0.0;
        for &g in offsets {
            let idx = g * samples + s;
            sum += soa.values[idx];
        }
        out[s] = sum / count;
    }
}

fn copy_single_gene_avx2(soa: &ExpressionSoAView, gene: usize, out: &mut [f32]) {
    let samples = soa.samples;
    if gene == usize::MAX {
        for value in out.iter_mut() {
            *value = 0.0;
        }
        return;
    }
    let start = gene * samples;
    let chunks = samples / 8 * 8;
    for s in (0..chunks).step_by(8) {
        let v = load_f32x8(&soa.values, start + s);
        store_f32x8(v, out, s);
    }
    for s in chunks..samples {
        out[s] = soa.values[start + s];
    }
}

fn variance4_avx2(a: &[f32], b: &[f32], c: &[f32], d: &[f32], out: &mut [f32]) {
    let samples = out.len();
    let chunks = samples / 8 * 8;
    let quarter = f32x8::new([0.25; 8]);

    for s in (0..chunks).step_by(8) {
        let av = load_f32x8(a, s);
        let bv = load_f32x8(b, s);
        let cv = load_f32x8(c, s);
        let dv = load_f32x8(d, s);
        let sum = (av + bv) + (cv + dv);
        let mean = sum * quarter;
        let da = av - mean;
        let db = bv - mean;
        let dc = cv - mean;
        let dd = dv - mean;
        let var = ((da * da) + (db * db) + (dc * dc) + (dd * dd)) * quarter;
        store_f32x8(var, out, s);
    }

    for s in chunks..samples {
        let mean = (a[s] + b[s] + c[s] + d[s]) * 0.25;
        let da = a[s] - mean;
        let db = b[s] - mean;
        let dc = c[s] - mean;
        let dd = d[s] - mean;
        out[s] = (da * da + db * db + dc * dc + dd * dd) * 0.25;
    }
}

#[inline]
fn load_f32x8(values: &[f32], start: usize) -> f32x8 {
    f32x8::new([
        values[start],
        values[start + 1],
        values[start + 2],
        values[start + 3],
        values[start + 4],
        values[start + 5],
        values[start + 6],
        values[start + 7],
    ])
}

#[inline]
fn store_f32x8(v: f32x8, out: &mut [f32], start: usize) {
    let lanes = v.to_array();
    out[start..start + 8].copy_from_slice(&lanes);
}
