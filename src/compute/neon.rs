//! NEON implementation for primitive signals.

#![allow(unsafe_op_in_unsafe_fn)]
#![cfg(target_arch = "aarch64")]

use std::arch::aarch64::*;

use crate::cache::ExpressionSoAView;
use crate::compute::{GeneOffsets, PrimitiveSignals};

pub fn compute_primitives_neon(soa: &ExpressionSoAView, offsets: &GeneOffsets) -> PrimitiveSignals {
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
        mean_over_offsets_neon(soa, &offsets.mtdna_all, &mut mtdna_mean);
        mean_over_offsets_neon(soa, &offsets.nuclear_oxphos, &mut nuclear_mean);

        mean_over_offsets_neon(soa, &offsets.complex_i, &mut c_i);
        mean_over_offsets_neon(soa, &offsets.complex_iii, &mut c_iii);
        mean_over_offsets_neon(soa, &offsets.complex_iv, &mut c_iv);
        mean_over_offsets_neon(soa, &offsets.complex_v, &mut c_v);

        mean_over_offsets_neon(soa, &offsets.ros, &mut ros_mean);
        mean_over_offsets_neon(soa, &offsets.mitophagy, &mut mitophagy_mean);
        mean_over_offsets_neon(soa, &offsets.fusion, &mut fusion_mean);
        mean_over_offsets_neon(soa, &offsets.fission, &mut fission_mean);
        mean_over_offsets_neon(soa, &offsets.biogenesis, &mut biogenesis_mean);

        mean_over_offsets_neon(soa, &offsets.atp_mt, &mut atp_mt);
        copy_single_gene_neon(soa, offsets.atp_nu, &mut atp_nu);

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

unsafe fn mean_over_offsets_neon(soa: &ExpressionSoAView, offsets: &[usize], out: &mut [f32]) {
    let samples = soa.samples;
    if offsets.is_empty() {
        for value in out.iter_mut() {
            *value = 0.0;
        }
        return;
    }

    let count = offsets.len() as f32;
    let denom = vdupq_n_f32(count);
    let chunks = samples / 4 * 4;

    for s in (0..chunks).step_by(4) {
        let mut acc = vdupq_n_f32(0.0);
        for &g in offsets {
            let base = g * samples + s;
            let ptr = soa.values.as_ptr().add(base);
            let v = vld1q_f32(ptr);
            acc = vaddq_f32(acc, v);
        }
        let mean = vdivq_f32(acc, denom);
        vst1q_f32(out.as_mut_ptr().add(s), mean);
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

unsafe fn copy_single_gene_neon(soa: &ExpressionSoAView, gene: usize, out: &mut [f32]) {
    let samples = soa.samples;
    if gene == usize::MAX {
        for value in out.iter_mut() {
            *value = 0.0;
        }
        return;
    }
    let start = gene * samples;
    let chunks = samples / 4 * 4;
    for s in (0..chunks).step_by(4) {
        let ptr = soa.values.as_ptr().add(start + s);
        let v = vld1q_f32(ptr);
        vst1q_f32(out.as_mut_ptr().add(s), v);
    }
    for s in chunks..samples {
        out[s] = soa.values[start + s];
    }
}

unsafe fn variance4_neon(a: &[f32], b: &[f32], c: &[f32], d: &[f32], out: &mut [f32]) {
    let samples = out.len();
    let chunks = samples / 4 * 4;
    let quarter = vdupq_n_f32(0.25);

    for s in (0..chunks).step_by(4) {
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
        let var = vmulq_f32(
            vaddq_f32(
                vaddq_f32(vmulq_f32(da, da), vmulq_f32(db, db)),
                vaddq_f32(vmulq_f32(dc, dc), vmulq_f32(dd, dd)),
            ),
            quarter,
        );
        vst1q_f32(out.as_mut_ptr().add(s), var);
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
