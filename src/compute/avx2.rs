//! AVX2 implementation for primitive signals.

#![allow(unsafe_op_in_unsafe_fn)]
#![cfg(target_feature = "avx2")]

use std::arch::x86_64::*;

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

    unsafe {
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
    }

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

unsafe fn mean_over_offsets_avx2(soa: &ExpressionSoAView, offsets: &[usize], out: &mut [f32]) {
    let samples = soa.samples;
    if offsets.is_empty() {
        for value in out.iter_mut() {
            *value = 0.0;
        }
        return;
    }

    let count = offsets.len() as f32;
    let denom = _mm256_set1_ps(count);
    let chunks = samples / 8 * 8;

    for s in (0..chunks).step_by(8) {
        let mut acc = _mm256_setzero_ps();
        for &g in offsets {
            let base = g * samples + s;
            let ptr = soa.values.as_ptr().add(base);
            let v = _mm256_loadu_ps(ptr);
            acc = _mm256_add_ps(acc, v);
        }
        let mean = _mm256_div_ps(acc, denom);
        _mm256_storeu_ps(out.as_mut_ptr().add(s), mean);
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

unsafe fn copy_single_gene_avx2(soa: &ExpressionSoAView, gene: usize, out: &mut [f32]) {
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
        let ptr = soa.values.as_ptr().add(start + s);
        let v = _mm256_loadu_ps(ptr);
        _mm256_storeu_ps(out.as_mut_ptr().add(s), v);
    }
    for s in chunks..samples {
        out[s] = soa.values[start + s];
    }
}

unsafe fn variance4_avx2(a: &[f32], b: &[f32], c: &[f32], d: &[f32], out: &mut [f32]) {
    let samples = out.len();
    let chunks = samples / 8 * 8;
    let quarter = _mm256_set1_ps(0.25);

    for s in (0..chunks).step_by(8) {
        let av = _mm256_loadu_ps(a.as_ptr().add(s));
        let bv = _mm256_loadu_ps(b.as_ptr().add(s));
        let cv = _mm256_loadu_ps(c.as_ptr().add(s));
        let dv = _mm256_loadu_ps(d.as_ptr().add(s));
        let sum = _mm256_add_ps(_mm256_add_ps(av, bv), _mm256_add_ps(cv, dv));
        let mean = _mm256_mul_ps(sum, quarter);
        let da = _mm256_sub_ps(av, mean);
        let db = _mm256_sub_ps(bv, mean);
        let dc = _mm256_sub_ps(cv, mean);
        let dd = _mm256_sub_ps(dv, mean);
        let var = _mm256_mul_ps(
            _mm256_add_ps(
                _mm256_add_ps(_mm256_mul_ps(da, da), _mm256_mul_ps(db, db)),
                _mm256_add_ps(_mm256_mul_ps(dc, dc), _mm256_mul_ps(dd, dd)),
            ),
            quarter,
        );
        _mm256_storeu_ps(out.as_mut_ptr().add(s), var);
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
