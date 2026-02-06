//! Aggregation helpers for expression matrices.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use sprs::CsMat;

use crate::input::InputError;

/// Aggregation mode for input data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AggregationMode {
    Sample,
    Cell,
    Cluster,
}

/// Mapping from cell index to cluster index with stable cluster ordering.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClusterMap {
    pub cluster_ids: Vec<String>,
    pub cell_to_cluster: Vec<usize>,
}

/// Aggregated dense matrix in gene-major order.
#[derive(Debug, Clone, PartialEq)]
pub struct AggregatedMatrix {
    pub values: Vec<f32>,
    pub genes: usize,
    pub samples: usize,
}

/// Load a cluster assignment file and align to barcodes.
pub fn load_cluster_map(path: &Path, barcodes: &[String]) -> Result<ClusterMap, InputError> {
    let file = File::open(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let reader = BufReader::new(file);

    let mut mapping: BTreeMap<String, String> = BTreeMap::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line.map_err(|source| InputError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 2 {
            return Err(InputError::InvalidClusterFile {
                path: path.to_path_buf(),
                message: format!("line {}: expected 2 columns", line_no + 1),
            });
        }
        let barcode = parts[0].trim();
        let cluster = parts[1].trim();
        if barcode.is_empty() || cluster.is_empty() {
            return Err(InputError::InvalidClusterFile {
                path: path.to_path_buf(),
                message: format!("line {}: empty barcode or cluster", line_no + 1),
            });
        }
        if mapping.contains_key(barcode) {
            return Err(InputError::InvalidClusterFile {
                path: path.to_path_buf(),
                message: format!("line {}: duplicate barcode {barcode}", line_no + 1),
            });
        }
        mapping.insert(barcode.to_string(), cluster.to_string());
    }

    if mapping.is_empty() {
        return Err(InputError::InvalidClusterFile {
            path: path.to_path_buf(),
            message: "no rows found".to_string(),
        });
    }

    let mut cluster_ids: BTreeSet<String> = BTreeSet::new();
    for barcode in barcodes {
        let cluster = mapping
            .get(barcode)
            .ok_or_else(|| InputError::InvalidClusterFile {
                path: path.to_path_buf(),
                message: format!("missing cluster for barcode {barcode}"),
            })?;
        cluster_ids.insert(cluster.clone());
    }

    let cluster_ids: Vec<String> = cluster_ids.into_iter().collect();
    let mut cluster_index: BTreeMap<String, usize> = BTreeMap::new();
    for (idx, id) in cluster_ids.iter().enumerate() {
        cluster_index.insert(id.clone(), idx);
    }

    let mut cell_to_cluster = Vec::with_capacity(barcodes.len());
    for barcode in barcodes {
        let cluster = mapping
            .get(barcode)
            .ok_or_else(|| InputError::InvalidClusterFile {
                path: path.to_path_buf(),
                message: format!("missing cluster for barcode {barcode}"),
            })?;
        let idx = cluster_index
            .get(cluster)
            .ok_or_else(|| InputError::InvalidClusterFile {
                path: path.to_path_buf(),
                message: format!("unknown cluster id {cluster}"),
            })?;
        cell_to_cluster.push(*idx);
    }

    Ok(ClusterMap {
        cluster_ids,
        cell_to_cluster,
    })
}

/// Aggregate a sparse matrix into dense gene-major order.
pub fn aggregate(
    matrix: &CsMat<f32>,
    mode: AggregationMode,
    clusters: Option<&ClusterMap>,
) -> AggregatedMatrix {
    let csc = matrix.to_csc();
    let genes = csc.rows();
    let cols = csc.cols();

    match mode {
        AggregationMode::Sample => {
            let samples = 1;
            let mut values = vec![0.0; genes];
            for col in csc.outer_iterator() {
                for (row, val) in col.iter() {
                    values[row] += val;
                }
            }
            if cols > 0 {
                let denom = cols as f32;
                for value in &mut values {
                    *value /= denom;
                }
            }
            AggregatedMatrix {
                values,
                genes,
                samples,
            }
        }
        AggregationMode::Cell => {
            let samples = cols;
            let mut values = vec![0.0; genes * samples];
            for (col_idx, col) in csc.outer_iterator().enumerate() {
                for (row, val) in col.iter() {
                    let idx = row * samples + col_idx;
                    values[idx] = *val;
                }
            }
            AggregatedMatrix {
                values,
                genes,
                samples,
            }
        }
        AggregationMode::Cluster => {
            let clusters = clusters.expect("cluster map required for cluster mode");
            let samples = clusters.cluster_ids.len();
            let mut values = vec![0.0; genes * samples];
            let mut counts = vec![0usize; samples];

            for (col_idx, col) in csc.outer_iterator().enumerate() {
                let cluster_idx = clusters.cell_to_cluster[col_idx];
                counts[cluster_idx] += 1;
                for (row, val) in col.iter() {
                    let idx = row * samples + cluster_idx;
                    values[idx] += val;
                }
            }

            for cluster_idx in 0..samples {
                let count = counts[cluster_idx];
                if count == 0 {
                    continue;
                }
                let denom = count as f32;
                for gene in 0..genes {
                    let idx = gene * samples + cluster_idx;
                    values[idx] /= denom;
                }
            }

            AggregatedMatrix {
                values,
                genes,
                samples,
            }
        }
    }
}
