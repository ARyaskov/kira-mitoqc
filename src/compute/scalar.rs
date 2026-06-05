//! Scalar reference primitives. Same outer tiling as SIMD paths.

use crate::cache::ExpressionSoAView;
use crate::compute::{GeneOffsets, PrimitiveSignals};

pub fn compute_primitives_scalar(
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

    tiled_means_scalar(values, samples, &mut panels);
    copy_single_gene(values, samples, offsets.atp_nu, &mut atp_nu);
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

fn tiled_means_scalar(values: &[f32], samples: usize, panels: &mut [(&[usize], &mut [f32])]) {
    let mut counts: [f32; 16] = [0.0; 16];
    debug_assert!(panels.len() <= counts.len());
    for (i, (offs, _)) in panels.iter().enumerate() {
        counts[i] = offs.len() as f32;
    }

    for s in 0..samples {
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

fn copy_single_gene(values: &[f32], samples: usize, gene: usize, out: &mut [f32]) {
    if gene == usize::MAX {
        out.fill(0.0);
        return;
    }
    let start = gene * samples;
    out.copy_from_slice(&values[start..start + samples]);
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
