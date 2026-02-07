//! Minimal aggregator contract outputs for pipeline mode.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use serde::Serialize;

use crate::output::OutputError;
use crate::output::profile::MitoProfileV1;

#[derive(Debug, Serialize)]
struct SummaryInput<'a> {
    mode: &'a str,
    n_samples: usize,
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
struct SummaryJson<'a> {
    tool: &'a str,
    input: SummaryInput<'a>,
    mitochondrial_state_distribution: BTreeMap<String, f64>,
    decay: SummaryDecay,
    axes_median: SummaryAxesMedian,
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
pub fn write_summary_json(out_dir: &Path, profiles: &[MitoProfileV1]) -> Result<(), OutputError> {
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
        "cell_id\tmitochondrial_state\tdecay_score\trobustness_margin\tbioenergetics\tros\tdynamics\tregulation"
    )
    .map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;

    for (sample, profile) in profiles.iter().enumerate() {
        writeln!(
            file,
            "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            barcodes[sample],
            profile.mitochondrial_state,
            profile.decay_score,
            profile.robustness_margin,
            profile.axes.bioenergetics,
            profile.axes.ros,
            profile.axes.dynamics,
            profile.axes.regulation
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
