//! Scalar reference implementation for primitive signals.

use crate::cache::ExpressionSoAView;
use crate::compute::{GeneOffsets, PrimitiveSignals};

pub fn compute_primitives_scalar(
    soa: &ExpressionSoAView,
    offsets: &GeneOffsets,
) -> PrimitiveSignals {
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

    mean_over_offsets(soa, &offsets.mtdna_all, &mut mtdna_mean);
    mean_over_offsets(soa, &offsets.nuclear_oxphos, &mut nuclear_mean);

    mean_over_offsets(soa, &offsets.complex_i, &mut c_i);
    mean_over_offsets(soa, &offsets.complex_iii, &mut c_iii);
    mean_over_offsets(soa, &offsets.complex_iv, &mut c_iv);
    mean_over_offsets(soa, &offsets.complex_v, &mut c_v);

    mean_over_offsets(soa, &offsets.ros, &mut ros_mean);
    mean_over_offsets(soa, &offsets.mitophagy, &mut mitophagy_mean);
    mean_over_offsets(soa, &offsets.fusion, &mut fusion_mean);
    mean_over_offsets(soa, &offsets.fission, &mut fission_mean);
    mean_over_offsets(soa, &offsets.biogenesis, &mut biogenesis_mean);

    mean_over_offsets(soa, &offsets.atp_mt, &mut atp_mt);
    copy_single_gene(soa, offsets.atp_nu, &mut atp_nu);

    compute_variance4(&c_i, &c_iii, &c_iv, &c_v, &mut stoich_variance);

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

fn mean_over_offsets(soa: &ExpressionSoAView, offsets: &[usize], out: &mut [f32]) {
    let samples = soa.samples;
    if offsets.is_empty() {
        for value in out.iter_mut() {
            *value = 0.0;
        }
        return;
    }

    let denom = offsets.len() as f32;
    for s in 0..samples {
        let mut sum = 0.0;
        for &g in offsets {
            let idx = g * samples + s;
            sum += soa.values[idx];
        }
        out[s] = sum / denom;
    }
}

fn copy_single_gene(soa: &ExpressionSoAView, gene: usize, out: &mut [f32]) {
    let samples = soa.samples;
    if gene == usize::MAX {
        for value in out.iter_mut() {
            *value = 0.0;
        }
        return;
    }
    let start = gene * samples;
    let slice = &soa.values[start..start + samples];
    out.copy_from_slice(slice);
}

fn compute_variance4(a: &[f32], b: &[f32], c: &[f32], d: &[f32], out: &mut [f32]) {
    for i in 0..out.len() {
        let mean = (a[i] + b[i] + c[i] + d[i]) * 0.25;
        let da = a[i] - mean;
        let db = b[i] - mean;
        let dc = c[i] - mean;
        let dd = d[i] - mean;
        out[i] = (da * da + db * db + dc * dc + dd * dd) * 0.25;
    }
}
