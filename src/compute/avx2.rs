//! AVX2 + FMA primitives, tiled by 8-sample chunks. Entry points are
//! `unsafe fn` + `#[target_feature]`; runtime gated via `dispatch::active_kind`.

#![cfg(target_arch = "x86_64")]

use std::arch::x86_64::*;

use crate::cache::ExpressionSoAView;
use crate::compute::{GeneOffsets, PrimitiveSignals};

#[target_feature(enable = "avx2,fma")]
pub unsafe fn compute_primitives_avx2(
    soa: &ExpressionSoAView<'_>,
    offsets: &GeneOffsets,
) -> PrimitiveSignals {
    let samples = soa.samples;

    let mut mtdna_mean = vec![0.0_f32; samples];
    let mut nuclear_mean = vec![0.0_f32; samples];
    let mut c_i = vec![0.0_f32; samples];
    let mut c_iii = vec![0.0_f32; samples];
    let mut c_iv = vec![0.0_f32; samples];
    let mut c_v = vec![0.0_f32; samples];
    let mut ros_mean = vec![0.0_f32; samples];
    let mut mitophagy_mean = vec![0.0_f32; samples];
    let mut fusion_mean = vec![0.0_f32; samples];
    let mut fission_mean = vec![0.0_f32; samples];
    let mut biogenesis_mean = vec![0.0_f32; samples];
    let mut atp_mt = vec![0.0_f32; samples];
    let mut atp_nu = vec![0.0_f32; samples];
    let mut stoich_variance = vec![0.0_f32; samples];

    let values = soa.values();

    let mut panels: [(&[usize], &mut [f32]); 12] = [
        (offsets.mtdna_all.as_slice(), mtdna_mean.as_mut_slice()),
        (offsets.nuclear_oxphos.as_slice(), nuclear_mean.as_mut_slice()),
        (offsets.complex_i.as_slice(), c_i.as_mut_slice()),
        (offsets.complex_iii.as_slice(), c_iii.as_mut_slice()),
        (offsets.complex_iv.as_slice(), c_iv.as_mut_slice()),
        (offsets.complex_v.as_slice(), c_v.as_mut_slice()),
        (offsets.ros.as_slice(), ros_mean.as_mut_slice()),
        (offsets.mitophagy.as_slice(), mitophagy_mean.as_mut_slice()),
        (offsets.fusion.as_slice(), fusion_mean.as_mut_slice()),
        (offsets.fission.as_slice(), fission_mean.as_mut_slice()),
        (offsets.biogenesis.as_slice(), biogenesis_mean.as_mut_slice()),
        (offsets.atp_mt.as_slice(), atp_mt.as_mut_slice()),
    ];

    unsafe {
        tiled_means_avx2(values, samples, &mut panels);

        copy_single_avx2(values, samples, offsets.atp_nu, &mut atp_nu);
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

/// Tiled SoA pass: per 8-sample chunk, sum every panel before advancing.
#[target_feature(enable = "avx2,fma")]
unsafe fn tiled_means_avx2(
    values: &[f32],
    samples: usize,
    panels: &mut [(&[usize], &mut [f32])],
) {
    let chunks = samples & !7;
    let s_ptr_base = values.as_ptr();

    let mut invs: [__m256; 16] = [_mm256_setzero_ps(); 16];
    let mut counts: [f32; 16] = [0.0; 16];
    debug_assert!(panels.len() <= invs.len());
    for (i, (offs, _)) in panels.iter().enumerate() {
        if !offs.is_empty() {
            let n = offs.len() as f32;
            counts[i] = n;
            invs[i] = _mm256_set1_ps(1.0 / n);
        }
    }

    let mut s = 0;
    while s < chunks {
        for (i, (offs, out)) in panels.iter_mut().enumerate() {
            if offs.is_empty() {
                unsafe {
                    _mm256_storeu_ps(out.as_mut_ptr().add(s), _mm256_setzero_ps());
                }
                continue;
            }
            let mut acc = _mm256_setzero_ps();
            for &g in *offs {
                let v = unsafe { _mm256_loadu_ps(s_ptr_base.add(g * samples + s)) };
                acc = _mm256_add_ps(acc, v);
            }
            let mean = _mm256_mul_ps(acc, invs[i]);
            unsafe {
                _mm256_storeu_ps(out.as_mut_ptr().add(s), mean);
            }
        }
        s += 8;
    }

    // Scalar tail.
    for s in chunks..samples {
        for (i, (offs, out)) in panels.iter_mut().enumerate() {
            if offs.is_empty() {
                out[s] = 0.0;
                continue;
            }
            let mut sum = 0.0_f32;
            for &g in *offs {
                sum += values[g * samples + s];
            }
            out[s] = sum / counts[i];
        }
    }
}

#[target_feature(enable = "avx2,fma")]
unsafe fn copy_single_avx2(values: &[f32], samples: usize, gene: usize, out: &mut [f32]) {
    if gene == usize::MAX {
        out.fill(0.0);
        return;
    }
    let start = gene * samples;
    let chunks = samples & !7;
    unsafe {
        let mut s = 0;
        while s < chunks {
            let v = _mm256_loadu_ps(values.as_ptr().add(start + s));
            _mm256_storeu_ps(out.as_mut_ptr().add(s), v);
            s += 8;
        }
    }
    for s in chunks..samples {
        out[s] = values[start + s];
    }
}

#[target_feature(enable = "avx2,fma")]
unsafe fn variance4_avx2(a: &[f32], b: &[f32], c: &[f32], d: &[f32], out: &mut [f32]) {
    let samples = out.len();
    let chunks = samples & !7;
    let quarter = _mm256_set1_ps(0.25);

    unsafe {
        let mut s = 0;
        while s < chunks {
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
            // FMA chain: ss = da² + db² + dc² + dd²
            let ss = _mm256_fmadd_ps(
                da,
                da,
                _mm256_fmadd_ps(db, db, _mm256_fmadd_ps(dc, dc, _mm256_mul_ps(dd, dd))),
            );
            let var = _mm256_mul_ps(ss, quarter);
            _mm256_storeu_ps(out.as_mut_ptr().add(s), var);
            s += 8;
        }
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
