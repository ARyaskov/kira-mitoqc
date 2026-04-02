use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::cache::ExpressionSoAView;
use crate::core::types::{ProxyKey, ProxyScores};
use crate::input::GeneIndex;
use crate::score::AxisScoresVec;
use crate::util::numeric::clamp01;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RedoxRegime {
    Baseline,
    CompensatedOxidativeStress,
    UnbufferedOxidativeStress,
    RedoxOverload,
}

impl RedoxRegime {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Baseline => "Baseline",
            Self::CompensatedOxidativeStress => "CompensatedOxidativeStress",
            Self::UnbufferedOxidativeStress => "UnbufferedOxidativeStress",
            Self::RedoxOverload => "RedoxOverload",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RedoxMetrics {
    pub mito_oxidative_stress_index: Vec<f32>,
    pub redox_buffering_capacity: Vec<f32>,
    pub mito_redox_mismatch: Vec<f32>,
    pub mitochondrial_stress_adaptation_score: Vec<f32>,
    pub redox_regime: Vec<RedoxRegime>,
    pub low_confidence: Vec<bool>,
}

#[derive(Debug, Error)]
pub enum RedoxError {
    #[error("resource read error at {path:?}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
}

const OXIDATIVE_PANEL_FALLBACK: &[&str] = &[
    "SOD2", "PRDX3", "TXN2", "GPX1", "GCLC", "GCLM", "HMOX1", "NQO1", "CAT", "SRXN1",
];
const BUFFERING_PANEL_FALLBACK: &[&str] = &[
    "TXN", "TXNRD1", "TXNRD2", "GSR", "SOD2", "PRDX3", "GPX4", "IDH2", "ME3", "NNT",
];

pub fn compute_redox_metrics(
    soa: &ExpressionSoAView,
    gene_index: &GeneIndex,
    axes: &AxisScoresVec,
    proxies: &ProxyScores,
) -> Result<RedoxMetrics, RedoxError> {
    let samples = soa.samples;

    let oxidative_genes = load_panel_genes("oxidative_stress_proxy.tsv", OXIDATIVE_PANEL_FALLBACK)?;
    let buffering_genes =
        load_panel_genes("redox_buffering_capacity.tsv", BUFFERING_PANEL_FALLBACK)?;

    let oxidative_offsets = resolve_panel_offsets(gene_index, &oxidative_genes, soa.genes);
    let buffering_offsets = resolve_panel_offsets(gene_index, &buffering_genes, soa.genes);

    let oxidative_cov = if oxidative_genes.is_empty() {
        0.0
    } else {
        (oxidative_offsets.len() as f32) / (oxidative_genes.len() as f32)
    };
    let buffering_cov = if buffering_genes.is_empty() {
        0.0
    } else {
        (buffering_offsets.len() as f32) / (buffering_genes.len() as f32)
    };

    let low_conf = oxidative_cov < 0.25 || buffering_cov < 0.25;

    let oxidative_mean = mean_over_offsets(soa, &oxidative_offsets);
    let buffering_mean = mean_over_offsets(soa, &buffering_offsets);

    let oxidative_norm = min_max_norm(&oxidative_mean);
    let buffering_norm = min_max_norm(&buffering_mean);

    let ros_proxy = get_proxy_or_zero(proxies, ProxyKey::RosResponseOverdrive, samples);
    let nadh_proxy = get_proxy_or_zero(proxies, ProxyKey::NadhImbalance, samples);
    let biogenesis_failure = get_proxy_or_zero(proxies, ProxyKey::BiogenesisFailure, samples);

    let mut mito_oxidative_stress_index = vec![0.0; samples];
    let mut redox_buffering_capacity = vec![0.0; samples];
    let mut mito_redox_mismatch = vec![0.0; samples];
    let mut mitochondrial_stress_adaptation_score = vec![0.0; samples];
    let mut redox_regime = vec![RedoxRegime::Baseline; samples];

    for i in 0..samples {
        let oxidative =
            clamp01(0.65 * oxidative_norm[i] + 0.20 * ros_proxy[i] + 0.15 * nadh_proxy[i]);
        let buffering = clamp01(
            0.70 * buffering_norm[i]
                + 0.20 * (1.0 - biogenesis_failure[i])
                + 0.10 * (1.0 - axes.ros[i]),
        );

        let mismatch = (oxidative - buffering).clamp(-1.0, 1.0);
        let mismatch01 = clamp01((mismatch + 1.0) * 0.5);
        let adaptation =
            clamp01(0.45 * axes.bioenergetics[i] + 0.25 * axes.regulation[i] + 0.30 * mismatch01);

        mito_oxidative_stress_index[i] = oxidative;
        redox_buffering_capacity[i] = buffering;
        mito_redox_mismatch[i] = mismatch;
        mitochondrial_stress_adaptation_score[i] = adaptation;

        redox_regime[i] = classify_redox_regime(oxidative, buffering, mismatch);
    }

    Ok(RedoxMetrics {
        mito_oxidative_stress_index,
        redox_buffering_capacity,
        mito_redox_mismatch,
        mitochondrial_stress_adaptation_score,
        redox_regime,
        low_confidence: vec![low_conf; samples],
    })
}

fn classify_redox_regime(oxidative: f32, buffering: f32, mismatch: f32) -> RedoxRegime {
    if mismatch >= 0.45 && oxidative >= 0.75 {
        return RedoxRegime::RedoxOverload;
    }
    if mismatch >= 0.20 && oxidative >= 0.55 {
        return RedoxRegime::UnbufferedOxidativeStress;
    }
    if oxidative >= 0.45 && buffering >= 0.45 {
        return RedoxRegime::CompensatedOxidativeStress;
    }
    RedoxRegime::Baseline
}

fn get_proxy_or_zero(proxies: &ProxyScores, key: ProxyKey, samples: usize) -> Vec<f32> {
    proxies
        .normalized
        .get(&key)
        .cloned()
        .unwrap_or_else(|| vec![0.0; samples])
}

fn min_max_norm(values: &[f32]) -> Vec<f32> {
    if values.is_empty() {
        return Vec::new();
    }
    let mut min_v = f32::INFINITY;
    let mut max_v = f32::NEG_INFINITY;
    for &v in values {
        if v < min_v {
            min_v = v;
        }
        if v > max_v {
            max_v = v;
        }
    }
    let span = max_v - min_v;
    if !span.is_finite() || span.abs() < 1e-12 {
        return vec![0.0; values.len()];
    }
    values.iter().map(|v| clamp01((v - min_v) / span)).collect()
}

fn mean_over_offsets(soa: &ExpressionSoAView, offsets: &[usize]) -> Vec<f32> {
    let samples = soa.samples;
    let mut out = vec![0.0; samples];
    if offsets.is_empty() {
        return out;
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
    out
}

fn resolve_panel_offsets(
    gene_index: &GeneIndex,
    genes: &[String],
    max_gene_rows: usize,
) -> Vec<usize> {
    let mut out = Vec::new();
    for gene in genes {
        if let Some(idx) = gene_index.get_index(gene) {
            if idx < max_gene_rows {
                out.push(idx);
            }
            continue;
        }
        let title = to_mouse_like(gene);
        if let Some(idx) = gene_index.get_index(&title) {
            if idx < max_gene_rows {
                out.push(idx);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn to_mouse_like(gene: &str) -> String {
    if gene.is_empty() {
        return String::new();
    }
    let mut chars = gene.chars();
    let first = chars.next().map(|c| c.to_ascii_uppercase()).unwrap_or(' ');
    let rest = chars.as_str().to_ascii_lowercase();
    format!("{first}{rest}")
}

fn load_panel_genes(file_name: &str, fallback: &[&str]) -> Result<Vec<String>, RedoxError> {
    if let Some(path) = resource_candidates(file_name)
        .into_iter()
        .find(|p| p.is_file())
    {
        let raw = std::fs::read_to_string(&path).map_err(|source| RedoxError::Io {
            path: path.clone(),
            source,
        })?;
        let mut out = Vec::new();
        for line in raw.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            out.push(trimmed.to_string());
        }
        if !out.is_empty() {
            return Ok(out);
        }
    }
    Ok(fallback.iter().map(|s| s.to_string()).collect())
}

fn resource_candidates(file_name: &str) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let base_manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let rel = Path::new("resources")
        .join("mitochondria")
        .join("redox")
        .join(file_name);

    out.push(base_manifest.join(&rel));
    if let Ok(cwd) = std::env::current_dir() {
        out.push(cwd.join(&rel));
    }
    out
}

#[cfg(test)]
mod tests {
    use crate::cache::{ExprCacheMode, mmap_expr_bin, write_expr_bin_with_mode};
    use crate::core::types::{ProxyKey, ProxyScores};
    use crate::data::ExpressionSoA;
    use crate::input::GeneIndex;
    use crate::score::AxisScoresVec;

    use super::{RedoxRegime, compute_redox_metrics};

    #[test]
    fn redox_stage_deterministic() {
        let genes = vec![
            "SOD2".to_string(),
            "PRDX3".to_string(),
            "TXN2".to_string(),
            "TXN".to_string(),
            "GSR".to_string(),
        ];
        let samples = 3;
        let mut values = Vec::new();
        for g in 0..genes.len() {
            for s in 0..samples {
                values.push((g as f32) * 0.1 + (s as f32) * 0.01);
            }
        }
        let soa = ExpressionSoA {
            values,
            genes: genes.len(),
            samples,
        };
        let path = std::env::temp_dir().join("kira_mitoqc_redox_stage.bin");
        write_expr_bin_with_mode(&path, &soa, ExprCacheMode::Cell).expect("write");
        let view = mmap_expr_bin(&path).expect("mmap");

        let index = GeneIndex::from_feature_list(&genes);
        let axes = AxisScoresVec {
            bioenergetics: vec![0.3, 0.6, 0.9],
            ros: vec![0.2, 0.6, 0.8],
            dynamics: vec![0.1, 0.1, 0.1],
            regulation: vec![0.3, 0.6, 0.9],
        };
        let mut proxies = ProxyScores::default();
        proxies.set(ProxyKey::RosResponseOverdrive, vec![0.2, 0.5, 0.9]);
        proxies.set(ProxyKey::NadhImbalance, vec![0.2, 0.4, 0.8]);
        proxies.set(ProxyKey::BiogenesisFailure, vec![0.2, 0.3, 0.6]);

        let a = compute_redox_metrics(&view, &index, &axes, &proxies).expect("redox");
        let b = compute_redox_metrics(&view, &index, &axes, &proxies).expect("redox");
        assert_eq!(a, b);
    }

    #[test]
    fn redox_regime_boundaries() {
        let genes = vec![
            "SOD2".to_string(),
            "PRDX3".to_string(),
            "TXN2".to_string(),
            "TXN".to_string(),
            "GSR".to_string(),
        ];
        let samples = 3;
        let mut values = Vec::new();
        for g in 0..genes.len() {
            for s in 0..samples {
                values.push(
                    if s == 0 {
                        0.8
                    } else if s == 1 {
                        0.6
                    } else {
                        0.2
                    } + (g as f32) * 0.01,
                );
            }
        }
        let soa = ExpressionSoA {
            values,
            genes: genes.len(),
            samples,
        };
        let path = std::env::temp_dir().join("kira_mitoqc_redox_regime.bin");
        write_expr_bin_with_mode(&path, &soa, ExprCacheMode::Cell).expect("write");
        let view = mmap_expr_bin(&path).expect("mmap");

        let index = GeneIndex::from_feature_list(&genes);
        let axes = AxisScoresVec {
            bioenergetics: vec![0.9, 0.7, 0.2],
            ros: vec![0.9, 0.7, 0.2],
            dynamics: vec![0.1, 0.1, 0.1],
            regulation: vec![0.9, 0.7, 0.2],
        };
        let mut proxies = ProxyScores::default();
        proxies.set(ProxyKey::RosResponseOverdrive, vec![0.95, 0.7, 0.2]);
        proxies.set(ProxyKey::NadhImbalance, vec![0.95, 0.6, 0.2]);
        proxies.set(ProxyKey::BiogenesisFailure, vec![0.8, 0.5, 0.1]);

        let metrics = compute_redox_metrics(&view, &index, &axes, &proxies).expect("redox");
        assert!(matches!(
            metrics.redox_regime[0],
            RedoxRegime::RedoxOverload | RedoxRegime::UnbufferedOxidativeStress
        ));
        assert!(matches!(
            metrics.redox_regime[2],
            RedoxRegime::Baseline | RedoxRegime::CompensatedOxidativeStress
        ));
    }
}
