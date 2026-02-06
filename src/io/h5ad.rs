//! H5AD loader (feature-gated).

use std::path::Path;

use hdf5::File;
use sprs::CsMat;

use crate::data::aggregate::ClusterMap;
use crate::input::InputError;

/// Loaded H5AD inputs.
#[derive(Debug, Clone)]
pub struct H5adInput {
    pub matrix: CsMat<f32>,
    pub features: Vec<String>,
    pub barcodes: Vec<String>,
}

/// Load H5AD expression matrix and metadata.
pub fn load_h5ad(path: &Path, gene_symbol_key: Option<&str>) -> Result<H5adInput, InputError> {
    let file = File::open(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    let (features, barcodes) = load_h5ad_metadata_from_file(&file, gene_symbol_key)?;
    let matrix = load_matrix(&file)?;

    if matrix.rows() != features.len() || matrix.cols() != barcodes.len() {
        return Err(InputError::DimensionMismatch {
            rows: matrix.rows(),
            cols: matrix.cols(),
            features: features.len(),
            barcodes: barcodes.len(),
        });
    }

    Ok(H5adInput {
        matrix,
        features,
        barcodes,
    })
}

/// Load H5AD features and barcodes without reading the matrix.
pub fn load_h5ad_metadata(
    path: &Path,
    gene_symbol_key: Option<&str>,
) -> Result<(Vec<String>, Vec<String>), InputError> {
    let file = File::open(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    load_h5ad_metadata_from_file(&file, gene_symbol_key)
}

fn load_h5ad_metadata_from_file(
    file: &File,
    gene_symbol_key: Option<&str>,
) -> Result<(Vec<String>, Vec<String>), InputError> {
    let barcodes = read_strings(file, "obs/_index")?;

    let features = if let Some(key) = gene_symbol_key {
        read_strings(file, &format!("var/{key}")).map_err(|_| InputError::InvalidGeneSymbolKey {
            key: key.to_string(),
        })?
    } else if let Ok(values) = read_strings(file, "var/gene_symbols") {
        values
    } else {
        read_strings(file, "var/_index")?
    };

    Ok((features, barcodes))
}

/// Load cluster labels from obs/<column> and build ClusterMap.
pub fn load_h5ad_clusters(path: &Path, column: &str) -> Result<ClusterMap, InputError> {
    let file = File::open(path).map_err(|source| InputError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let barcodes = read_strings(&file, "obs/_index")?;
    let labels = read_strings(&file, &format!("obs/{column}"))?;
    if labels.len() != barcodes.len() {
        return Err(InputError::InvalidClusterFile {
            path: path.to_path_buf(),
            message: "cluster labels length mismatch".to_string(),
        });
    }
    Ok(build_cluster_map(&barcodes, &labels))
}

fn read_strings(file: &File, path: &str) -> Result<Vec<String>, InputError> {
    let dataset = file
        .dataset(path)
        .map_err(|_| InputError::MissingH5adDataset {
            path: path.to_string(),
        })?;
    dataset
        .read_1d::<String>()
        .map_err(|_| InputError::MissingH5adDataset {
            path: path.to_string(),
        })
}

fn load_matrix(file: &File) -> Result<CsMat<f32>, InputError> {
    if let Ok(dataset) = file.dataset("X") {
        let shape = dataset.shape();
        if shape.len() != 2 {
            return Err(InputError::UnsupportedH5adMatrix {
                layout: "dense: non-2d".to_string(),
            });
        }
        let rows = shape[0] as usize;
        let cols = shape[1] as usize;
        let data = read_f32_vec(&dataset)?;
        return Ok(dense_to_csr(rows, cols, &data));
    }

    let group = file
        .group("X")
        .map_err(|_| InputError::MissingH5adDataset {
            path: "X".to_string(),
        })?;

    let encoding = read_attr_string(&group, "encoding-type")
        .or_else(|_| read_attr_string(&group, "encoding_type"))
        .unwrap_or_else(|_| "".to_string());

    let data =
        read_f32_vec(
            &group
                .dataset("data")
                .map_err(|_| InputError::MissingH5adDataset {
                    path: "X/data".to_string(),
                })?,
        )?;
    let indices =
        read_usize_vec(
            &group
                .dataset("indices")
                .map_err(|_| InputError::MissingH5adDataset {
                    path: "X/indices".to_string(),
                })?,
        )?;
    let indptr =
        read_usize_vec(
            &group
                .dataset("indptr")
                .map_err(|_| InputError::MissingH5adDataset {
                    path: "X/indptr".to_string(),
                })?,
        )?;

    let shape = read_shape_attr(&group)?;
    let rows = shape.0;
    let cols = shape.1;

    match encoding.as_str() {
        "csr_matrix" => Ok(CsMat::new_csr((rows, cols), indptr, indices, data)),
        "csc_matrix" => Ok(CsMat::new_csc((rows, cols), indptr, indices, data)),
        _ => Err(InputError::UnsupportedH5adMatrix { layout: encoding }),
    }
}

fn read_attr_string(group: &hdf5::Group, name: &str) -> Result<String, InputError> {
    let attr = group
        .attr(name)
        .map_err(|_| InputError::MissingH5adDataset {
            path: format!("X/@{name}"),
        })?;
    attr.read_scalar::<String>()
        .map_err(|_| InputError::MissingH5adDataset {
            path: format!("X/@{name}"),
        })
}

fn read_shape_attr(group: &hdf5::Group) -> Result<(usize, usize), InputError> {
    let attr = group
        .attr("shape")
        .map_err(|_| InputError::MissingH5adDataset {
            path: "X/@shape".to_string(),
        })?;
    let shape: Vec<i64> = attr
        .read_raw()
        .map_err(|_| InputError::MissingH5adDataset {
            path: "X/@shape".to_string(),
        })?;
    if shape.len() != 2 {
        return Err(InputError::UnsupportedH5adMatrix {
            layout: "shape attr length".to_string(),
        });
    }
    Ok((shape[0] as usize, shape[1] as usize))
}

fn read_f32_vec(dataset: &hdf5::Dataset) -> Result<Vec<f32>, InputError> {
    if let Ok(values) = dataset.read_raw::<f32>() {
        return Ok(values);
    }
    let values = dataset
        .read_raw::<f64>()
        .map_err(|_| InputError::MissingH5adDataset {
            path: dataset.name().unwrap_or("dataset").to_string(),
        })?;
    Ok(values.into_iter().map(|v| v as f32).collect())
}

fn read_usize_vec(dataset: &hdf5::Dataset) -> Result<Vec<usize>, InputError> {
    if let Ok(values) = dataset.read_raw::<i64>() {
        return Ok(values.into_iter().map(|v| v as usize).collect());
    }
    let values = dataset
        .read_raw::<i32>()
        .map_err(|_| InputError::MissingH5adDataset {
            path: dataset.name().unwrap_or("dataset").to_string(),
        })?;
    Ok(values.into_iter().map(|v| v as usize).collect())
}

fn dense_to_csr(rows: usize, cols: usize, data: &[f32]) -> CsMat<f32> {
    let mut indptr = Vec::with_capacity(rows + 1);
    let mut indices = Vec::with_capacity(rows * cols);
    let mut values = Vec::with_capacity(rows * cols);

    indptr.push(0);
    for r in 0..rows {
        let row_start = r * cols;
        for c in 0..cols {
            indices.push(c);
            values.push(data[row_start + c]);
        }
        indptr.push(values.len());
    }

    CsMat::new_csr((rows, cols), indptr, indices, values)
}

fn build_cluster_map(barcodes: &[String], labels: &[String]) -> ClusterMap {
    use std::collections::{BTreeMap, BTreeSet};

    let mut cluster_ids: BTreeSet<String> = BTreeSet::new();
    for label in labels {
        cluster_ids.insert(label.clone());
    }

    let cluster_ids: Vec<String> = cluster_ids.into_iter().collect();
    let mut cluster_index = BTreeMap::new();
    for (idx, id) in cluster_ids.iter().enumerate() {
        cluster_index.insert(id.clone(), idx);
    }

    let mut cell_to_cluster = Vec::with_capacity(barcodes.len());
    for label in labels {
        let idx = cluster_index.get(label).copied().unwrap_or(0);
        cell_to_cluster.push(idx);
    }

    ClusterMap {
        cluster_ids,
        cell_to_cluster,
    }
}
