use std::collections::BTreeSet;

use serde::Serialize;

use crate::cache::ExpressionSoAView;
use crate::input::{GeneIndex, ResolvedGeneSets};
use crate::metrics::metabolic_extension::panels::{
    BIOGENESIS_PANEL, FAO_PANEL, GLYCOLYSIS_PANEL, OXPHOS_PANEL, ROS_PANEL, panel_alias,
    to_mouse_like,
};

const MIN_GENES: usize = 3;
const TRIM_FRAC: f32 = 0.1;
const EPS: f32 = 1e-6;

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct MetabolicThresholds {
    pub metabolic_rigid_high: f32,
    pub ros_high: f32,
    pub energetic_strain_high: f32,
    pub compensation_failure: f32,
    pub ogi_oxphos_dominant: f32,
    pub ogi_glycolysis_dominant: f32,
}

impl Default for MetabolicThresholds {
    fn default() -> Self {
        Self {
            metabolic_rigid_high: 2.0,
            ros_high: 2.0,
            energetic_strain_high: 1.5,
            compensation_failure: -1.5,
            ogi_oxphos_dominant: 1.5,
            ogi_glycolysis_dominant: -1.5,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MetabolicMissingness {
    pub oxphos_found: usize,
    pub oxphos_total: usize,
    pub glycolysis_found: usize,
    pub glycolysis_total: usize,
    pub fao_found: usize,
    pub fao_total: usize,
    pub ros_found: usize,
    pub ros_total: usize,
    pub biogenesis_found: usize,
    pub biogenesis_total: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetabolicMetrics {
    pub oxphos_core: Vec<f32>,
    pub gly_core: Vec<f32>,
    pub fao_core: Vec<f32>,
    pub ros_core: Vec<f32>,
    pub bio_core: Vec<f32>,
    pub mri: Vec<f32>,
    pub osl: Vec<f32>,
    pub ess: Vec<f32>,
    pub mcb: Vec<f32>,
    pub ogi: Vec<f32>,
    pub metabolic_rigid_high: Vec<bool>,
    pub ros_high: Vec<bool>,
    pub energetic_strain_high: Vec<bool>,
    pub compensation_failure: Vec<bool>,
    pub oxphos_dominant: Vec<bool>,
    pub glycolysis_dominant: Vec<bool>,
    pub thresholds: MetabolicThresholds,
    pub missingness: MetabolicMissingness,
}

pub fn compute_metabolic_metrics(
    soa: &ExpressionSoAView,
    gene_index: &GeneIndex,
    resolved: &ResolvedGeneSets,
    mito_stress_signal: &[f32],
) -> MetabolicMetrics {
    let base = existing_gene_count(resolved);
    let oxphos_rows = resolve_panel_rows(OXPHOS_PANEL, gene_index, base);
    let gly_rows = resolve_panel_rows(GLYCOLYSIS_PANEL, gene_index, base + OXPHOS_PANEL.len());
    let fao_rows = resolve_panel_rows(
        FAO_PANEL,
        gene_index,
        base + OXPHOS_PANEL.len() + GLYCOLYSIS_PANEL.len(),
    );
    let ros_rows = resolve_panel_rows(
        ROS_PANEL,
        gene_index,
        base + OXPHOS_PANEL.len() + GLYCOLYSIS_PANEL.len() + FAO_PANEL.len(),
    );
    let bio_rows = resolve_panel_rows(
        BIOGENESIS_PANEL,
        gene_index,
        base + OXPHOS_PANEL.len() + GLYCOLYSIS_PANEL.len() + FAO_PANEL.len() + ROS_PANEL.len(),
    );

    let oxphos_core = panel_trimmed_mean(soa, &oxphos_rows);
    let gly_core = panel_trimmed_mean(soa, &gly_rows);
    let fao_core = panel_trimmed_mean(soa, &fao_rows);
    let ros_core = panel_trimmed_mean(soa, &ros_rows);
    let bio_core = panel_trimmed_mean(soa, &bio_rows);

    let zox = robust_zscore(&oxphos_core);
    let zgly = robust_zscore(&gly_core);
    let zfao = robust_zscore(&fao_core);
    let zros = robust_zscore(&ros_core);
    let zbio = robust_zscore(&bio_core);
    let zmito = robust_zscore(mito_stress_signal);

    let n = soa.samples;
    let mut mri = vec![f32::NAN; n];
    let mut osl = vec![f32::NAN; n];
    let mut ess = vec![f32::NAN; n];
    let mut mcb = vec![f32::NAN; n];
    let mut ogi = vec![f32::NAN; n];

    for i in 0..n {
        let ox = zox[i];
        let gy = zgly[i];
        let fa = zfao[i];
        let rs = zros[i];
        let bg = zbio[i];

        if ox.is_finite() && gy.is_finite() && fa.is_finite() {
            let dom = ox.max(gy).max(fa);
            let mean = (ox + gy + fa) / 3.0;
            let var =
                ((ox - mean) * (ox - mean) + (gy - mean) * (gy - mean) + (fa - mean) * (fa - mean))
                    / 3.0;
            mri[i] = dom - var.max(0.0).sqrt();
        }

        if rs.is_finite() {
            let mito = zmito[i];
            osl[i] = if mito.is_finite() {
                0.6 * rs + 0.4 * mito
            } else {
                rs
            };
        }

        if rs.is_finite() && ox.is_finite() && gy.is_finite() {
            let supply = ox.max(gy);
            ess[i] = (rs - supply).max(0.0);
        }

        if bg.is_finite() && rs.is_finite() {
            mcb[i] = bg - rs;
        }

        if ox.is_finite() && gy.is_finite() {
            ogi[i] = ox - gy;
        }
    }

    let thresholds = MetabolicThresholds::default();
    let metabolic_rigid_high = mri
        .iter()
        .map(|v| v.is_finite() && *v >= thresholds.metabolic_rigid_high)
        .collect();
    let ros_high = osl
        .iter()
        .map(|v| v.is_finite() && *v >= thresholds.ros_high)
        .collect();
    let energetic_strain_high = ess
        .iter()
        .map(|v| v.is_finite() && *v >= thresholds.energetic_strain_high)
        .collect();
    let compensation_failure = mcb
        .iter()
        .map(|v| v.is_finite() && *v <= thresholds.compensation_failure)
        .collect();
    let oxphos_dominant = ogi
        .iter()
        .map(|v| v.is_finite() && *v >= thresholds.ogi_oxphos_dominant)
        .collect();
    let glycolysis_dominant = ogi
        .iter()
        .map(|v| v.is_finite() && *v <= thresholds.ogi_glycolysis_dominant)
        .collect();

    MetabolicMetrics {
        oxphos_core,
        gly_core,
        fao_core,
        ros_core,
        bio_core,
        mri,
        osl,
        ess,
        mcb,
        ogi,
        metabolic_rigid_high,
        ros_high,
        energetic_strain_high,
        compensation_failure,
        oxphos_dominant,
        glycolysis_dominant,
        thresholds,
        missingness: MetabolicMissingness {
            oxphos_found: oxphos_rows.len(),
            oxphos_total: OXPHOS_PANEL.len(),
            glycolysis_found: gly_rows.len(),
            glycolysis_total: GLYCOLYSIS_PANEL.len(),
            fao_found: fao_rows.len(),
            fao_total: FAO_PANEL.len(),
            ros_found: ros_rows.len(),
            ros_total: ROS_PANEL.len(),
            biogenesis_found: bio_rows.len(),
            biogenesis_total: BIOGENESIS_PANEL.len(),
        },
    }
}

fn existing_gene_count(resolved: &ResolvedGeneSets) -> usize {
    resolved.mtdna_complex_i.genes.len()
        + resolved.mtdna_complex_iii.genes.len()
        + resolved.mtdna_complex_iv.genes.len()
        + resolved.mtdna_complex_v.genes.len()
        + resolved.nuclear_oxphos_complex_i.genes.len()
        + resolved.nuclear_oxphos_complex_ii.genes.len()
        + resolved.nuclear_oxphos_complex_iii.genes.len()
        + resolved.nuclear_oxphos_complex_iv.genes.len()
        + resolved.nuclear_oxphos_complex_v.genes.len()
        + resolved.ros.genes.len()
        + resolved.mitophagy.genes.len()
        + resolved.fusion.genes.len()
        + resolved.fission.genes.len()
        + resolved.biogenesis.genes.len()
}

fn resolve_panel_rows(panel: &[&str], gene_index: &GeneIndex, base_offset: usize) -> Vec<usize> {
    let mut rows = Vec::with_capacity(panel.len());
    let mut seen = BTreeSet::new();

    for (pos, symbol) in panel.iter().enumerate() {
        let mut present = gene_index.get_index(symbol).is_some()
            || gene_index.get_index(&to_mouse_like(symbol)).is_some();

        if !present {
            if let Some(alias) = panel_alias(symbol) {
                present = gene_index.get_index(alias).is_some()
                    || gene_index.get_index(&to_mouse_like(alias)).is_some();
            }
        }

        if present {
            let row = base_offset + pos;
            if seen.insert(row) {
                rows.push(row);
            }
        }
    }
    rows
}

fn panel_trimmed_mean(soa: &ExpressionSoAView, rows: &[usize]) -> Vec<f32> {
    let n = soa.samples;
    let mut out = vec![f32::NAN; n];
    if rows.len() < MIN_GENES {
        return out;
    }

    let mut buf = Vec::<f32>::with_capacity(rows.len());
    for s in 0..n {
        buf.clear();
        for &g in rows {
            let idx = g * n + s;
            let v = soa.values[idx].max(0.0).ln_1p();
            if v.is_finite() {
                buf.push(v);
            }
        }
        if buf.len() < MIN_GENES {
            continue;
        }
        out[s] = trimmed_mean_in_place(&mut buf, TRIM_FRAC);
    }
    out
}

fn robust_zscore(values: &[f32]) -> Vec<f32> {
    let mut finite = Vec::new();
    for &v in values {
        if v.is_finite() {
            finite.push(v);
        }
    }

    let mut out = vec![f32::NAN; values.len()];
    if finite.is_empty() {
        return out;
    }

    let med = median(&mut finite);
    let mut deviations: Vec<f32> = values
        .iter()
        .filter_map(|v| {
            if v.is_finite() {
                Some((v - med).abs())
            } else {
                None
            }
        })
        .collect();
    let mad = median(&mut deviations);

    if mad <= 0.0 || !mad.is_finite() {
        for (i, &v) in values.iter().enumerate() {
            if v.is_finite() {
                out[i] = 0.0;
            }
        }
        return out;
    }

    let denom = 1.4826 * mad + EPS;
    for (i, &v) in values.iter().enumerate() {
        if v.is_finite() {
            out[i] = (v - med) / denom;
        }
    }
    out
}

fn trimmed_mean_in_place(values: &mut [f32], trim_frac: f32) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let n = values.len();
    let k = ((n as f32) * trim_frac).floor() as usize;
    let start = k.min(n);
    let end = n.saturating_sub(k);
    if start >= end {
        return f32::NAN;
    }
    let mut sum = 0.0;
    for v in &values[start..end] {
        sum += *v;
    }
    sum / ((end - start) as f32)
}

fn median(values: &mut [f32]) -> f32 {
    if values.is_empty() {
        return f32::NAN;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        0.5 * (values[n / 2 - 1] + values[n / 2])
    }
}

#[cfg(test)]
mod tests {
    use super::{median, trimmed_mean_in_place};

    #[test]
    fn trimmed_mean_drops_tails() {
        let mut v = vec![
            1.0, 2.0, 3.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0,
        ];
        let tm = trimmed_mean_in_place(&mut v, 0.1);
        assert!((tm - 263.125).abs() < 1e-6);
    }

    #[test]
    fn median_even_count() {
        let mut v = vec![1.0, 4.0, 2.0, 3.0];
        let m = median(&mut v);
        assert_eq!(m, 2.5);
    }
}
