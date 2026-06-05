//! NEON primitives, tiled by 4-sample chunks. Mirror of `avx2.rs`.

#![cfg(target_arch = "aarch64")]

use std::arch::aarch64::*;

use crate::cache::ExpressionSoAView;
use crate::compute::{GeneOffsets, PrimitiveSignals};

pub fn compute_primitives_neon(
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
        tiled_means_neon(values, samples, &mut panels);
        copy_single_neon(values, samples, offsets.atp_nu, &mut atp_nu);
        variance4_neon(&c_i, &c_iii, &c_iv, &c_v, &mut stoich_variance);
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

#[target_feature(enable = "neon")]
unsafe fn tiled_means_neon(
    values: &[f32],
    samples: usize,
    panels: &mut [(&[usize], &mut [f32])],
) {
    let chunks = samples & !3;
    let s_ptr_base = values.as_ptr();

    let mut invs: [float32x4_t; 16] = [vdupq_n_f32(0.0); 16];
    let mut counts: [f32; 16] = [0.0; 16];
    debug_assert!(panels.len() <= invs.len());
    for (i, (offs, _)) in panels.iter().enumerate() {
        if !offs.is_empty() {
            let n = offs.len() as f32;
            counts[i] = n;
            invs[i] = vdupq_n_f32(1.0 / n);
        }
    }

    let mut s = 0;
    while s < chunks {
        for (i, (offs, out)) in panels.iter_mut().enumerate() {
            if offs.is_empty() {
                unsafe {
                    vst1q_f32(out.as_mut_ptr().add(s), vdupq_n_f32(0.0));
                }
                continue;
            }
            let mut acc = vdupq_n_f32(0.0);
            for &g in *offs {
                let v = unsafe { vld1q_f32(s_ptr_base.add(g * samples + s)) };
                acc = vaddq_f32(acc, v);
            }
            let mean = vmulq_f32(acc, invs[i]);
            unsafe {
                vst1q_f32(out.as_mut_ptr().add(s), mean);
            }
        }
        s += 4;
    }

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

#[target_feature(enable = "neon")]
unsafe fn copy_single_neon(values: &[f32], samples: usize, gene: usize, out: &mut [f32]) {
    if gene == usize::MAX {
        out.fill(0.0);
        return;
    }
    let start = gene * samples;
    let chunks = samples & !3;
    unsafe {
        let mut s = 0;
        while s < chunks {
            let v = vld1q_f32(values.as_ptr().add(start + s));
            vst1q_f32(out.as_mut_ptr().add(s), v);
            s += 4;
        }
    }
    for s in chunks..samples {
        out[s] = values[start + s];
    }
}

#[target_feature(enable = "neon")]
unsafe fn variance4_neon(a: &[f32], b: &[f32], c: &[f32], d: &[f32], out: &mut [f32]) {
    let samples = out.len();
    let chunks = samples & !3;
    let quarter = unsafe { vdupq_n_f32(0.25) };

    unsafe {
        let mut s = 0;
        while s < chunks {
            let av = vld1q_f32(a.as_ptr().add(s));
            let bv = vld1q_f32(b.as_ptr().add(s));
            let cv = vld1q_f32(c.as_ptr().add(s));
            let dv = vld1q_f32(d.as_ptr().add(s));
            let sum = vaddq_f32(vaddq_f32(av, bv), vaddq_f32(cv, dv));
            let mean = vmulq_f32(sum, quarter);
            let da = vsubq_f32(av, mean);
            let db = vsubq_f32(bv, mean);
            let dc = vsubq_f32(cv, mean);
            let dd = vsubq_f32(dv, mean);
            // FMA chain: ss = da*da + (db*db + (dc*dc + dd*dd))
            let ss = vfmaq_f32(
                vfmaq_f32(vfmaq_f32(vmulq_f32(dd, dd), dc, dc), db, db),
                da,
                da,
            );
            let var = vmulq_f32(ss, quarter);
            vst1q_f32(out.as_mut_ptr().add(s), var);
            s += 4;
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
