//! Output writers for kira-mitoqc.

use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::core::types::{ProxyKey, ProxyScores};
use crate::output::profile::MitoProfileV1;
use crate::score::{AxisScoresVec, DecayScoreVec};

pub mod profile;
pub mod v2;

/// Output errors.
#[derive(Debug, Error)]
pub enum OutputError {
    #[error("cannot create output directory {path:?}: {source}")]
    CreateDir {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("cannot write file {path:?}: {source}")]
    WriteFile {
        path: PathBuf,
        source: std::io::Error,
    },

    #[error("serialization error for {path:?}: {source}")]
    Serialize {
        path: PathBuf,
        source: serde_json::Error,
    },
}

/// Write JSON output to `mitochondrial_profile.json`.
pub fn write_json(out_dir: &Path, profiles: &[MitoProfileV1]) -> Result<(), OutputError> {
    fs::create_dir_all(out_dir).map_err(|source| OutputError::CreateDir {
        path: out_dir.to_path_buf(),
        source,
    })?;
    let path = out_dir.join("mitochondrial_profile.json");
    let file = File::create(&path).map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;
    serde_json::to_writer_pretty(file, profiles)
        .map_err(|source| OutputError::Serialize { path, source })?;
    Ok(())
}

/// Write axis scores to `axes.tsv`.
pub fn write_axes_tsv(out_dir: &Path, axes: &AxisScoresVec) -> Result<(), OutputError> {
    fs::create_dir_all(out_dir).map_err(|source| OutputError::CreateDir {
        path: out_dir.to_path_buf(),
        source,
    })?;
    let path = out_dir.join("axes.tsv");
    let mut file = File::create(&path).map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;

    writeln!(file, "sample\tbioenergetics\tros\tdynamics\tregulation").map_err(|source| {
        OutputError::WriteFile {
            path: path.clone(),
            source,
        }
    })?;

    for i in 0..axes.bioenergetics.len() {
        writeln!(
            file,
            "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
            i, axes.bioenergetics[i], axes.ros[i], axes.dynamics[i], axes.regulation[i]
        )
        .map_err(|source| OutputError::WriteFile {
            path: path.clone(),
            source,
        })?;
    }

    Ok(())
}

/// Write decay scores to `decay.tsv`.
pub fn write_decay_tsv(out_dir: &Path, decay: &DecayScoreVec) -> Result<(), OutputError> {
    fs::create_dir_all(out_dir).map_err(|source| OutputError::CreateDir {
        path: out_dir.to_path_buf(),
        source,
    })?;
    let path = out_dir.join("decay.tsv");
    let mut file = File::create(&path).map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;

    writeln!(file, "sample\tdecay_score\trobustness_margin").map_err(|source| {
        OutputError::WriteFile {
            path: path.clone(),
            source,
        }
    })?;

    for i in 0..decay.decay.len() {
        writeln!(
            file,
            "{}\t{:.6}\t{:.6}",
            i, decay.decay[i], decay.robustness_margin[i]
        )
        .map_err(|source| OutputError::WriteFile {
            path: path.clone(),
            source,
        })?;
    }

    Ok(())
}

/// Write normalized proxy scores to `proxies.tsv`.
pub fn write_proxies_tsv(out_dir: &Path, proxies: &ProxyScores) -> Result<(), OutputError> {
    fs::create_dir_all(out_dir).map_err(|source| OutputError::CreateDir {
        path: out_dir.to_path_buf(),
        source,
    })?;
    let path = out_dir.join("proxies.tsv");
    let mut file = File::create(&path).map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;

    let proxy_order = [
        ProxyKey::ETCStoichiometryLoss,
        ProxyKey::MtdnaExpressionUncoupling,
        ProxyKey::AtpCouplingLoss,
        ProxyKey::RosResponseOverdrive,
        ProxyKey::NadhImbalance,
        ProxyKey::DynamicsImbalance,
        ProxyKey::MitophagyExcess,
        ProxyKey::BiogenesisFailure,
    ];

    let header = proxy_order
        .iter()
        .map(|k| k.as_str())
        .collect::<Vec<_>>()
        .join("\t");
    writeln!(file, "sample\t{header}").map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;

    let values = proxy_order
        .iter()
        .map(|k| proxies.normalized.get(k).expect("missing proxy values"))
        .collect::<Vec<_>>();

    let samples = values.first().map(|v| v.len()).unwrap_or(0);
    for i in 0..samples {
        write!(file, "{}", i).map_err(|source| OutputError::WriteFile {
            path: path.clone(),
            source,
        })?;
        for vec in &values {
            write!(file, "\t{:.6}", vec[i]).map_err(|source| OutputError::WriteFile {
                path: path.clone(),
                source,
            })?;
        }
        writeln!(file).map_err(|source| OutputError::WriteFile {
            path: path.clone(),
            source,
        })?;
    }

    Ok(())
}

/// Write v2 JSON output to `mitochondrial_profile.v2.json`.
pub fn write_json_v2(
    out_dir: &Path,
    profiles: &[v2::MitoProfileBundleV2],
) -> Result<(), OutputError> {
    fs::create_dir_all(out_dir).map_err(|source| OutputError::CreateDir {
        path: out_dir.to_path_buf(),
        source,
    })?;
    let path = out_dir.join("mitochondrial_profile.v2.json");
    let file = File::create(&path).map_err(|source| OutputError::WriteFile {
        path: path.clone(),
        source,
    })?;
    serde_json::to_writer_pretty(file, profiles)
        .map_err(|source| OutputError::Serialize { path, source })?;
    Ok(())
}
