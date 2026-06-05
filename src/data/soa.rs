//! Structure-of-arrays expression representation.

use rustc_hash::FxHashMap;

use crate::data::aggregate::{AggregationMode, ClusterMap};
use crate::input::{GeneIndex, GeneResolution, InputError, ResolvedGeneSets};
use crate::io::mtx::MtxInput;
use crate::metrics::metabolic_extension::panels::{
    BIOGENESIS_PANEL, FAO_PANEL, GLYCOLYSIS_PANEL, OXPHOS_PANEL, ROS_PANEL, panel_alias,
    to_mouse_like,
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

/// Maps gene symbols to row indices in the prepared SoA.
///
/// **Must** be used instead of `GeneIndex` (which indexes input features).
#[derive(Debug, Clone, Default)]
pub struct SoaIndex {
    map: rustc_hash::FxHashMap<String, usize>,
    rows: usize,
}

impl SoaIndex {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn get(&self, symbol: &str) -> Option<usize> {
        self.map.get(symbol).copied()
    }

    pub fn from_resolved(resolved: &ResolvedGeneSets) -> Self {
        Self::from_ordered(&ordered_gene_symbols(resolved))
    }

    /// Build from an arbitrary ordered list. Position = SoA row; first
    /// occurrence wins.
    pub fn from_ordered<S: AsRef<str>>(ordered: &[S]) -> Self {
        let mut map: rustc_hash::FxHashMap<String, usize> =
            rustc_hash::FxHashMap::with_capacity_and_hasher(ordered.len(), Default::default());
        for (row, name) in ordered.iter().enumerate() {
            map.entry(name.as_ref().to_string()).or_insert(row);
        }
        Self {
            rows: ordered.len(),
            map,
        }
    }
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

/// Streams the input CSC once and writes only the ~200 resolved rows
/// (not the full ~30k feature set).
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
    let gene_index = GeneIndex::from_feature_list(&mtx.features);

    let ordered_genes = ordered_gene_symbols(resolved);
    let total_genes = ordered_genes.len();

    // Determine output sample dimension.
    let cols = mtx.matrix.cols();
    let samples = match mode {
        AggregationMode::Sample => 1,
        AggregationMode::Cell => cols,
        AggregationMode::Cluster => clusters
            .map(|c| c.cluster_ids.len())
            .expect("checked above for Cluster mode"),
    };

    // input_row -> SoA output row. First out-row wins on collision.
    let mut input_to_out: FxHashMap<usize, usize> =
        FxHashMap::with_capacity_and_hasher(total_genes, Default::default());
    for (out_row, symbol) in ordered_genes.iter().enumerate() {
        if let Some(input_row) = resolve_symbol(&gene_index, symbol) {
            input_to_out.entry(input_row).or_insert(out_row);
        }
    }

    let mut values = vec![0.0_f32; total_genes * samples];

    let csc = mtx.matrix.to_csc();
    let cluster_assign = clusters.map(|c| c.cell_to_cluster.as_slice());

    for (cell_idx, col) in csc.outer_iterator().enumerate() {
        let sample_idx = match mode {
            AggregationMode::Sample => 0,
            AggregationMode::Cell => cell_idx,
            AggregationMode::Cluster => cluster_assign
                .expect("checked above for Cluster mode")[cell_idx],
        };

        for (row, value) in col.iter() {
            if let Some(&out_row) = input_to_out.get(&row) {
                let idx = out_row * samples + sample_idx;
                match mode {
                    AggregationMode::Cell => values[idx] = *value,
                    _ => values[idx] += *value,
                }
            }
        }
    }

    match mode {
        AggregationMode::Sample => {
            if cols > 0 {
                let denom = cols as f32;
                for v in &mut values {
                    *v /= denom;
                }
            }
        }
        AggregationMode::Cluster => {
            let assign =
                cluster_assign.expect("checked above for Cluster mode");
            let mut counts = vec![0usize; samples];
            for &cl in assign {
                counts[cl] += 1;
            }
            for (cluster_idx, count) in counts.iter().enumerate() {
                if *count == 0 {
                    continue;
                }
                let denom = *count as f32;
                for gene in 0..total_genes {
                    values[gene * samples + cluster_idx] /= denom;
                }
            }
        }
        AggregationMode::Cell => {}
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

/// Lookup chain: exact → alias → mouse-case → mouse-case(alias).
fn resolve_symbol(gene_index: &GeneIndex, symbol: &str) -> Option<usize> {
    if let Some(idx) = gene_index.get_index(symbol) {
        return Some(idx);
    }
    if let Some(alias) = panel_alias(symbol) {
        if let Some(idx) = gene_index.get_index(alias) {
            return Some(idx);
        }
        let mouse_alias = to_mouse_like(alias);
        if let Some(idx) = gene_index.get_index(&mouse_alias) {
            return Some(idx);
        }
    }
    let mouse = to_mouse_like(symbol);
    gene_index.get_index(&mouse)
}

pub(crate) fn ordered_gene_symbols(resolved: &ResolvedGeneSets) -> Vec<String> {
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
