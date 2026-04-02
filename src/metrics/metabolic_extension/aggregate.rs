use std::collections::BTreeMap;

use serde::Serialize;

use crate::metrics::metabolic_extension::panels::MITO_METABOLIC_PANEL_V1;
use crate::metrics::metabolic_extension::scores::{
    MetabolicMetrics, MetabolicMissingness, MetabolicThresholds,
};

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MetabolicGlobalStats {
    pub mri_median: Option<f64>,
    pub mri_mad: Option<f64>,
    pub osl_median: Option<f64>,
    pub osl_mad: Option<f64>,
    pub ess_median: Option<f64>,
    pub ess_mad: Option<f64>,
    pub mcb_median: Option<f64>,
    pub mcb_mad: Option<f64>,
    pub ogi_median: Option<f64>,
    pub ogi_mad: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MetabolicClusterStats {
    pub cluster_id: String,
    pub mri_median: Option<f64>,
    pub mri_p10: Option<f64>,
    pub mri_p90: Option<f64>,
    pub osl_median: Option<f64>,
    pub osl_p10: Option<f64>,
    pub osl_p90: Option<f64>,
    pub ess_median: Option<f64>,
    pub ess_p10: Option<f64>,
    pub ess_p90: Option<f64>,
    pub mcb_median: Option<f64>,
    pub mcb_p10: Option<f64>,
    pub mcb_p90: Option<f64>,
    pub metabolic_rigid_high_fraction: f64,
    pub ros_high_fraction: f64,
    pub energetic_strain_high_fraction: f64,
    pub compensation_failure_fraction: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TopClusterMetric {
    pub cluster_id: String,
    pub value: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MetabolicSummary {
    pub panel_version: &'static str,
    pub thresholds: MetabolicThresholds,
    pub global_stats: MetabolicGlobalStats,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub cluster_stats: Vec<MetabolicClusterStats>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub top_clusters_by_mri: Vec<TopClusterMetric>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub top_clusters_by_ess: Vec<TopClusterMetric>,
    pub missingness: MetabolicMissingness,
}

pub fn build_summary(
    metrics: &MetabolicMetrics,
    sample_ids: &[String],
    include_cluster_stats: bool,
) -> MetabolicSummary {
    let global_stats = MetabolicGlobalStats {
        mri_median: median_opt(&metrics.mri),
        mri_mad: mad_opt(&metrics.mri),
        osl_median: median_opt(&metrics.osl),
        osl_mad: mad_opt(&metrics.osl),
        ess_median: median_opt(&metrics.ess),
        ess_mad: mad_opt(&metrics.ess),
        mcb_median: median_opt(&metrics.mcb),
        mcb_mad: mad_opt(&metrics.mcb),
        ogi_median: median_opt(&metrics.ogi),
        ogi_mad: mad_opt(&metrics.ogi),
    };

    let mut cluster_stats = Vec::new();
    if include_cluster_stats {
        let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (idx, sample_id) in sample_ids.iter().enumerate() {
            groups.entry(sample_id.clone()).or_default().push(idx);
        }

        for (cluster_id, idxs) in groups {
            cluster_stats.push(cluster_stat_row(&cluster_id, &idxs, metrics));
        }
    }

    let top_clusters_by_mri = top_clusters(&cluster_stats, |r| r.mri_median, 5);
    let top_clusters_by_ess = top_clusters(&cluster_stats, |r| r.ess_median, 5);

    MetabolicSummary {
        panel_version: MITO_METABOLIC_PANEL_V1,
        thresholds: metrics.thresholds,
        global_stats,
        cluster_stats,
        top_clusters_by_mri,
        top_clusters_by_ess,
        missingness: metrics.missingness.clone(),
    }
}

fn cluster_stat_row(
    cluster_id: &str,
    idxs: &[usize],
    metrics: &MetabolicMetrics,
) -> MetabolicClusterStats {
    let mri = select(&metrics.mri, idxs);
    let osl = select(&metrics.osl, idxs);
    let ess = select(&metrics.ess, idxs);
    let mcb = select(&metrics.mcb, idxs);

    let denom = idxs.len() as f64;
    MetabolicClusterStats {
        cluster_id: cluster_id.to_string(),
        mri_median: median_opt(&mri),
        mri_p10: percentile_opt(&mri, 0.10),
        mri_p90: percentile_opt(&mri, 0.90),
        osl_median: median_opt(&osl),
        osl_p10: percentile_opt(&osl, 0.10),
        osl_p90: percentile_opt(&osl, 0.90),
        ess_median: median_opt(&ess),
        ess_p10: percentile_opt(&ess, 0.10),
        ess_p90: percentile_opt(&ess, 0.90),
        mcb_median: median_opt(&mcb),
        mcb_p10: percentile_opt(&mcb, 0.10),
        mcb_p90: percentile_opt(&mcb, 0.90),
        metabolic_rigid_high_fraction: fraction_true(&metrics.metabolic_rigid_high, idxs, denom),
        ros_high_fraction: fraction_true(&metrics.ros_high, idxs, denom),
        energetic_strain_high_fraction: fraction_true(&metrics.energetic_strain_high, idxs, denom),
        compensation_failure_fraction: fraction_true(&metrics.compensation_failure, idxs, denom),
    }
}

fn top_clusters(
    rows: &[MetabolicClusterStats],
    value_fn: impl Fn(&MetabolicClusterStats) -> Option<f64>,
    n: usize,
) -> Vec<TopClusterMetric> {
    let mut values: Vec<TopClusterMetric> = rows
        .iter()
        .filter_map(|row| {
            value_fn(row).map(|v| TopClusterMetric {
                cluster_id: row.cluster_id.clone(),
                value: round6(v),
            })
        })
        .collect();

    values.sort_by(|a, b| {
        b.value
            .total_cmp(&a.value)
            .then_with(|| a.cluster_id.cmp(&b.cluster_id))
    });
    values.truncate(n);
    values
}

fn select(values: &[f32], idxs: &[usize]) -> Vec<f32> {
    idxs.iter()
        .map(|&i| values.get(i).copied().unwrap_or(f32::NAN))
        .collect()
}

fn fraction_true(flags: &[bool], idxs: &[usize], denom: f64) -> f64 {
    if denom <= 0.0 {
        return 0.0;
    }
    let mut count = 0usize;
    for &idx in idxs {
        if flags.get(idx).copied().unwrap_or(false) {
            count += 1;
        }
    }
    round6((count as f64) / denom)
}

fn finite_values(values: &[f32]) -> Vec<f64> {
    values
        .iter()
        .filter_map(|v| if v.is_finite() { Some(*v as f64) } else { None })
        .collect()
}

fn median_opt(values: &[f32]) -> Option<f64> {
    let mut finite = finite_values(values);
    if finite.is_empty() {
        return None;
    }
    finite.sort_by(|a, b| a.total_cmp(b));
    let n = finite.len();
    let med = if n % 2 == 1 {
        finite[n / 2]
    } else {
        0.5 * (finite[n / 2 - 1] + finite[n / 2])
    };
    Some(round6(med))
}

fn mad_opt(values: &[f32]) -> Option<f64> {
    let mut finite = finite_values(values);
    if finite.is_empty() {
        return None;
    }
    finite.sort_by(|a, b| a.total_cmp(b));
    let n = finite.len();
    let med = if n % 2 == 1 {
        finite[n / 2]
    } else {
        0.5 * (finite[n / 2 - 1] + finite[n / 2])
    };

    let mut dev: Vec<f64> = finite.into_iter().map(|v| (v - med).abs()).collect();
    dev.sort_by(|a, b| a.total_cmp(b));
    let dn = dev.len();
    let mad = if dn % 2 == 1 {
        dev[dn / 2]
    } else {
        0.5 * (dev[dn / 2 - 1] + dev[dn / 2])
    };
    Some(round6(mad))
}

fn percentile_opt(values: &[f32], q: f64) -> Option<f64> {
    let mut finite = finite_values(values);
    if finite.is_empty() {
        return None;
    }
    finite.sort_by(|a, b| a.total_cmp(b));
    let n = finite.len();
    let idx = ((n - 1) as f64 * q).round() as usize;
    Some(round6(finite[idx]))
}

fn round6(v: f64) -> f64 {
    (v * 1_000_000.0).round() / 1_000_000.0
}
