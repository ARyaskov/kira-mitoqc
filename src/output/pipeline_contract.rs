//! Minimal aggregator contract outputs for pipeline mode.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use serde::Serialize;

use crate::input::ExpressionSource;
use crate::metrics::metabolic_extension::aggregate::MetabolicSummary;
use crate::metrics::metabolic_extension::scores::MetabolicMetrics;
use crate::output::OutputError;
use crate::output::profile::MitoProfileV1;
use crate::redox::{RedoxMetrics, RedoxRegime};

#[derive(Debug, Serialize)]
struct SummaryInput<'a> {
    mode: &'a str,
    n_samples: usize,
    input_format: &'a str,
    expression_type: ExpressionSource,
}

#[derive(Debug, Serialize)]
struct SummaryDecay {
    decay_score_median: f64,
    robustness_margin_median: f64,
}

#[derive(Debug, Serialize)]
struct SummaryAxesMedian {
    bioenergetics: f64,
    ros: f64,
    dynamics: f64,
    regulation: f64,
}

#[derive(Debug, Serialize)]
struct SummaryRedox {
    regime_fractions: BTreeMap<String, f64>,
    mean_mito_redox_mismatch: f64,
    high_redox_overload_fraction: f64,
}

#[derive(Debug, Serialize)]
struct SummaryJson<'a> {
    tool: &'a str,
    input: SummaryInput<'a>,
    mitochondrial_state_distribution: BTreeMap<String, f64>,
    decay: SummaryDecay,
    axes_median: SummaryAxesMedian,
    #[serde(skip_serializing_if = "Option::is_none")]
    redox: Option<SummaryRedox>,
    #[serde(skip_serializing_if = "Option::is_none")]
    mitochondrial_metabolic: Option<&'a MetabolicSummary>,
}

#[derive(Debug, Serialize)]
struct SampleMetricsSpec<'a> {
    file: &'a str,
    id_column: &'a str,
    state_column: &'a str,
}

#[derive(Debug, Serialize)]
struct PipelineArtifacts<'a> {
    summary: &'a str,
    primary_metrics: &'a str,
    shared_cache: &'a str,
}

#[derive(Debug, Serialize)]
struct PipelineStepJson<'a> {
    tool: &'a str,
    mode: &'a str,
    artifacts: PipelineArtifacts<'a>,
    sample_metrics: SampleMetricsSpec<'a>,
    axes: [&'a str; 4],
}

/// Write `summary.json` for pipeline aggregator ingestion.
pub fn write_summary_json(
    out_dir: &Path,
    profiles: &[MitoProfileV1],
    input_format: &str,
    expression_type: ExpressionSource,
) -> Result<(), OutputError> {
    write_summary_json_with_redox(out_dir, profiles, input_format, expression_type, None, None)
}

/// Write `summary.json` for pipeline aggregator ingestion with optional redox extension.
pub fn write_summary_json_with_redox(
    out_dir: &Path,
    profiles: &[MitoProfileV1],
    input_format: &str,
    expression_type: ExpressionSource,
    redox: Option<&RedoxMetrics>,
    metabolic_summary: Option<&MetabolicSummary>,
) -> Result<(), OutputError> {
    fs::create_dir_all(out_dir).map_err(|source| OutputError::CreateDir {
        path: out_dir.to_path_buf(),
        source,
    })?;

    let path = out_dir.join("summary.json");
    let file = File::create(&path).map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;

    let summary = SummaryJson {
        tool: "kira-mitoqc",
        input: SummaryInput {
            mode: "pipeline",
            n_samples: profiles.len(),
            input_format,
            expression_type,
        },
        mitochondrial_state_distribution: state_distribution(profiles),
        decay: SummaryDecay {
            decay_score_median: median(profiles.iter().map(|p| p.decay_score)),
            robustness_margin_median: median(profiles.iter().map(|p| p.robustness_margin)),
        },
        axes_median: SummaryAxesMedian {
            bioenergetics: median(profiles.iter().map(|p| p.axes.bioenergetics)),
            ros: median(profiles.iter().map(|p| p.axes.ros)),
            dynamics: median(profiles.iter().map(|p| p.axes.dynamics)),
            regulation: median(profiles.iter().map(|p| p.axes.regulation)),
        },
        redox: redox.map(build_redox_summary),
        mitochondrial_metabolic: metabolic_summary,
    };

    serde_json::to_writer_pretty(file, &summary)
        .map_err(|source| OutputError::Serialize { path, source })?;
    Ok(())
}

/// Write `mito_metrics.tsv` for downstream aggregator consumption.
pub fn write_mito_metrics_tsv(
    out_dir: &Path,
    barcodes: &[String],
    profiles: &[MitoProfileV1],
    metabolic: &MetabolicMetrics,
) -> Result<(), OutputError> {
    fs::create_dir_all(out_dir).map_err(|source| OutputError::CreateDir {
        path: out_dir.to_path_buf(),
        source,
    })?;

    let path = out_dir.join("mito_metrics.tsv");
    let mut file = File::create(&path).map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;

    writeln!(
        file,
        "cell_id\tmitochondrial_state\tdecay_score\trobustness_margin\tbioenergetics\tros\tdynamics\tregulation\toxphos_core\tgly_core\tfao_core\tros_core\tbio_core\tMRI\tOSL\tESS\tMCB\tOGI\tmetabolic_rigid_high\tros_high\tenergetic_strain_high\tcompensation_failure"
    )
    .map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;

    for (sample, profile) in profiles.iter().enumerate() {
        writeln!(
            file,
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}",
            barcodes[sample],
            profile.mitochondrial_state,
            profile.decay_score,
            profile.robustness_margin,
            profile.axes.bioenergetics,
            profile.axes.ros,
            profile.axes.dynamics,
            profile.axes.regulation,
            metabolic.oxphos_core.get(sample).copied().unwrap_or(f32::NAN),
            metabolic.gly_core.get(sample).copied().unwrap_or(f32::NAN),
            metabolic.fao_core.get(sample).copied().unwrap_or(f32::NAN),
            metabolic.ros_core.get(sample).copied().unwrap_or(f32::NAN),
            metabolic.bio_core.get(sample).copied().unwrap_or(f32::NAN),
            metabolic.mri.get(sample).copied().unwrap_or(f32::NAN),
            metabolic.osl.get(sample).copied().unwrap_or(f32::NAN),
            metabolic.ess.get(sample).copied().unwrap_or(f32::NAN),
            metabolic.mcb.get(sample).copied().unwrap_or(f32::NAN),
            metabolic.ogi.get(sample).copied().unwrap_or(f32::NAN),
            metabolic
                .metabolic_rigid_high
                .get(sample)
                .copied()
                .unwrap_or(false),
            metabolic.ros_high.get(sample).copied().unwrap_or(false),
            metabolic
                .energetic_strain_high
                .get(sample)
                .copied()
                .unwrap_or(false),
            metabolic
                .compensation_failure
                .get(sample)
                .copied()
                .unwrap_or(false)
        )
        .map_err(|source| OutputError::WriteFile {
            path: path.clone(),
            source,
        })?;
    }

    Ok(())
}

/// Write `pipeline_step.json` manifest for pipeline mode.
pub fn write_pipeline_step_json(out_dir: &Path, shared_cache: &str) -> Result<(), OutputError> {
    fs::create_dir_all(out_dir).map_err(|source| OutputError::CreateDir {
        path: out_dir.to_path_buf(),
        source,
    })?;
    let path = out_dir.join("pipeline_step.json");
    let file = File::create(&path).map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;

    let payload = PipelineStepJson {
        tool: "kira-mitoqc",
        mode: "pipeline",
        artifacts: PipelineArtifacts {
            summary: "summary.json",
            primary_metrics: "mito_metrics.tsv",
            shared_cache,
        },
        sample_metrics: SampleMetricsSpec {
            file: "mito_metrics.tsv",
            id_column: "cell_id",
            state_column: "mitochondrial_state",
        },
        axes: ["bioenergetics", "ros", "dynamics", "regulation"],
    };
    serde_json::to_writer_pretty(file, &payload)
        .map_err(|source| OutputError::Serialize { path, source })?;
    Ok(())
}

fn state_distribution(profiles: &[MitoProfileV1]) -> BTreeMap<String, f64> {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for profile in profiles {
        *counts
            .entry(profile.mitochondrial_state.clone())
            .or_insert(0usize) += 1;
    }

    let n_samples = profiles.len() as f64;
    let mut out = BTreeMap::new();
    for (state, count) in counts {
        let frac = if n_samples == 0.0 {
            0.0
        } else {
            (count as f64) / n_samples
        };
        out.insert(state, round6(frac));
    }
    out
}

fn build_redox_summary(redox: &RedoxMetrics) -> SummaryRedox {
    let n = redox.redox_regime.len() as f64;
    let mut counts = BTreeMap::<String, usize>::new();
    counts.insert(
        RedoxRegime::CompensatedOxidativeStress.as_str().to_string(),
        0,
    );
    counts.insert(
        RedoxRegime::UnbufferedOxidativeStress.as_str().to_string(),
        0,
    );
    counts.insert(RedoxRegime::RedoxOverload.as_str().to_string(), 0);

    for regime in &redox.redox_regime {
        match regime {
            RedoxRegime::CompensatedOxidativeStress
            | RedoxRegime::UnbufferedOxidativeStress
            | RedoxRegime::RedoxOverload => {
                *counts.entry(regime.as_str().to_string()).or_insert(0) += 1;
            }
            RedoxRegime::Baseline => {}
        }
    }

    let mut regime_fractions = BTreeMap::new();
    for (k, v) in counts {
        let frac = if n == 0.0 { 0.0 } else { v as f64 / n };
        regime_fractions.insert(k, round6(frac));
    }

    let mean_mismatch = if redox.mito_redox_mismatch.is_empty() {
        0.0
    } else {
        let sum: f64 = redox.mito_redox_mismatch.iter().map(|v| *v as f64).sum();
        round6(sum / redox.mito_redox_mismatch.len() as f64)
    };

    let overload_count = redox
        .redox_regime
        .iter()
        .filter(|r| matches!(r, RedoxRegime::RedoxOverload))
        .count();
    let overload_frac = if n == 0.0 {
        0.0
    } else {
        round6(overload_count as f64 / n)
    };

    SummaryRedox {
        regime_fractions,
        mean_mito_redox_mismatch: mean_mismatch,
        high_redox_overload_fraction: overload_frac,
    }
}

fn median(values: impl Iterator<Item = f32>) -> f64 {
    let mut sorted = values.map(|v| v as f64).collect::<Vec<_>>();
    if sorted.is_empty() {
        return 0.0;
    }
    sorted.sort_by(f64::total_cmp);

    let mid = sorted.len() / 2;
    let raw = if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    };
    round6(raw)
}

fn round6(v: f64) -> f64 {
    (v * 1_000_000.0).round() / 1_000_000.0
}
