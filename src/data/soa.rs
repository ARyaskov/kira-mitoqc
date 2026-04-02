//! Structure-of-arrays expression representation.

use crate::data::aggregate::{AggregatedMatrix, AggregationMode, ClusterMap, aggregate};
use crate::input::{GeneIndex, GeneResolution, InputError, ResolvedGeneSets};
use crate::io::mtx::MtxInput;
use crate::metrics::metabolic_extension::panels::{
    BIOGENESIS_PANEL, FAO_PANEL, GLYCOLYSIS_PANEL, OXPHOS_PANEL, ROS_PANEL,
};

/// Dense expression values laid out as [gene][sample].
#[derive(Debug, Clone, PartialEq)]
pub struct ExpressionSoA {
    /// values laid out as [gene][sample]
    pub values: Vec<f32>,
    pub genes: usize,
    pub samples: usize,
}

impl ExpressionSoA {
    /// Get a value at (gene, sample).
    pub fn get(&self, gene: usize, sample: usize) -> f32 {
        let idx = gene * self.samples + sample;
        self.values[idx]
    }
}

/// Prepared expression values ready for downstream computation.
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedExpression {
    pub soa: ExpressionSoA,
    pub resolved: ResolvedGeneSets,
}

/// Prepare expression SoA for the resolved gene sets and aggregation mode.
pub fn prepare_expression(
    mtx: &MtxInput,
    resolved: &ResolvedGeneSets,
    mode: AggregationMode,
) -> Result<PreparedExpression, InputError> {
    if matches!(mode, AggregationMode::Cluster) {
        return Err(InputError::InvalidClusterFile {
            path: std::path::PathBuf::from("<in-memory>"),
            message: "cluster map required for cluster mode".to_string(),
        });
    }
    prepare_expression_with_clusters(mtx, resolved, mode, None)
}

/// Prepare expression SoA with an optional cluster map.
pub fn prepare_expression_with_clusters(
    mtx: &MtxInput,
    resolved: &ResolvedGeneSets,
    mode: AggregationMode,
    clusters: Option<&ClusterMap>,
) -> Result<PreparedExpression, InputError> {
    if matches!(mode, AggregationMode::Cluster) && clusters.is_none() {
        return Err(InputError::InvalidClusterFile {
            path: std::path::PathBuf::from("<in-memory>"),
            message: "cluster map required for cluster mode".to_string(),
        });
    }
    let aggregated = aggregate(&mtx.matrix, mode, clusters);
    let gene_index = GeneIndex::from_feature_list(&mtx.features);

    let ordered_genes = ordered_gene_symbols(resolved);
    let samples = aggregated.samples;
    let total_genes = ordered_genes.len();
    let mut values = vec![0.0; total_genes * samples];

    for (out_gene, symbol) in ordered_genes.iter().enumerate() {
        if let Some(idx) = gene_index.get_index(symbol) {
            copy_gene_row(&aggregated, idx, out_gene, &mut values);
        }
    }

    Ok(PreparedExpression {
        soa: ExpressionSoA {
            values,
            genes: total_genes,
            samples,
        },
        resolved: resolved.clone(),
    })
}

fn copy_gene_row(
    aggregated: &AggregatedMatrix,
    input_gene: usize,
    out_gene: usize,
    out: &mut [f32],
) {
    let samples = aggregated.samples;
    let src_start = input_gene * samples;
    let dst_start = out_gene * samples;
    let src = &aggregated.values[src_start..src_start + samples];
    let dst = &mut out[dst_start..dst_start + samples];
    dst.copy_from_slice(src);
}

fn ordered_gene_symbols(resolved: &ResolvedGeneSets) -> Vec<String> {
    let mut ordered = Vec::new();
    append_genes(&mut ordered, &resolved.mtdna_complex_i);
    append_genes(&mut ordered, &resolved.mtdna_complex_iii);
    append_genes(&mut ordered, &resolved.mtdna_complex_iv);
    append_genes(&mut ordered, &resolved.mtdna_complex_v);
    append_genes(&mut ordered, &resolved.nuclear_oxphos_complex_i);
    append_genes(&mut ordered, &resolved.nuclear_oxphos_complex_ii);
    append_genes(&mut ordered, &resolved.nuclear_oxphos_complex_iii);
    append_genes(&mut ordered, &resolved.nuclear_oxphos_complex_iv);
    append_genes(&mut ordered, &resolved.nuclear_oxphos_complex_v);
    append_genes(&mut ordered, &resolved.ros);
    append_genes(&mut ordered, &resolved.mitophagy);
    append_genes(&mut ordered, &resolved.fusion);
    append_genes(&mut ordered, &resolved.fission);
    append_genes(&mut ordered, &resolved.biogenesis);
    ordered.extend(OXPHOS_PANEL.iter().map(|g| (*g).to_string()));
    ordered.extend(GLYCOLYSIS_PANEL.iter().map(|g| (*g).to_string()));
    ordered.extend(FAO_PANEL.iter().map(|g| (*g).to_string()));
    ordered.extend(ROS_PANEL.iter().map(|g| (*g).to_string()));
    ordered.extend(BIOGENESIS_PANEL.iter().map(|g| (*g).to_string()));
    ordered
}

fn append_genes(target: &mut Vec<String>, resolution: &GeneResolution) {
    target.extend(resolution.genes.iter().cloned());
}
