//! Gene symbol to column index resolver.

use std::collections::BTreeMap;

use tracing::warn;

use crate::core::types::GeneSet;
use crate::input::InputError;

/// Gene symbol to column index mapping.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneIndex {
    /// Gene symbol -> column index in input matrix.
    map: BTreeMap<String, usize>,
}

impl GeneIndex {
    /// Build from a list of feature symbols. First occurrence wins.
    pub fn from_feature_list(features: &[String]) -> Self {
        let mut map = BTreeMap::new();
        for (idx, symbol) in features.iter().enumerate() {
            if map.contains_key(symbol) {
                warn!(
                    symbol = symbol.as_str(),
                    "Duplicate gene symbol in feature list"
                );
                continue;
            }
            map.insert(symbol.clone(), idx);
        }
        Self { map }
    }

    /// Build from a list of feature symbols with validation.
    pub fn try_from_feature_list(features: &[String]) -> Result<Self, InputError> {
        if features.is_empty() {
            return Err(InputError::EmptyFeatureList);
        }
        Ok(Self::from_feature_list(features))
    }

    /// Resolve a list of gene symbols into column indices.
    pub fn resolve(&self, genes: &[String]) -> GeneResolution {
        let mut found = Vec::with_capacity(genes.len());
        let mut missing = Vec::new();

        for gene in genes {
            if let Some(&idx) = self.map.get(gene) {
                found.push(idx);
            } else {
                missing.push(gene.clone());
            }
        }

        GeneResolution {
            genes: genes.to_vec(),
            found,
            missing,
        }
    }

    /// Get a column index for a gene symbol.
    pub fn get_index(&self, gene: &str) -> Option<usize> {
        self.map.get(gene).copied()
    }
}

/// Resolution result for a gene set.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneResolution {
    pub genes: Vec<String>,
    pub found: Vec<usize>,
    pub missing: Vec<String>,
}

/// QC counters derived from a resolution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeneResolutionQC {
    pub found: usize,
    pub missing: usize,
}

impl GeneResolutionQC {
    /// Compute QC counters from a gene resolution.
    pub fn from_resolution(resolution: &GeneResolution) -> Self {
        Self {
            found: resolution.found.len(),
            missing: resolution.missing.len(),
        }
    }
}

/// Resolved gene sets for all v1 sections.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedGeneSets {
    pub mtdna_complex_i: GeneResolution,
    pub mtdna_complex_iii: GeneResolution,
    pub mtdna_complex_iv: GeneResolution,
    pub mtdna_complex_v: GeneResolution,
    pub mtdna_all: GeneResolution,
    pub nuclear_oxphos_complex_i: GeneResolution,
    pub nuclear_oxphos_complex_ii: GeneResolution,
    pub nuclear_oxphos_complex_iii: GeneResolution,
    pub nuclear_oxphos_complex_iv: GeneResolution,
    pub nuclear_oxphos_complex_v: GeneResolution,
    pub nuclear_oxphos_all: GeneResolution,
    pub ros: GeneResolution,
    pub mitophagy: GeneResolution,
    pub fusion: GeneResolution,
    pub fission: GeneResolution,
    pub biogenesis: GeneResolution,
}

/// Resolve all v1 gene sets at once.
pub fn resolve_all_genesets(gene_index: &GeneIndex, geneset: &GeneSet) -> ResolvedGeneSets {
    let mtdna_all: Vec<String> = geneset
        .all_mtdna()
        .into_iter()
        .map(|s| s.to_string())
        .collect();
    let nuclear_oxphos_all: Vec<String> = geneset
        .all_nuclear_oxphos()
        .into_iter()
        .map(|s| s.to_string())
        .collect();

    ResolvedGeneSets {
        mtdna_complex_i: gene_index.resolve(&geneset.mtdna_complex_i),
        mtdna_complex_iii: gene_index.resolve(&geneset.mtdna_complex_iii),
        mtdna_complex_iv: gene_index.resolve(&geneset.mtdna_complex_iv),
        mtdna_complex_v: gene_index.resolve(&geneset.mtdna_complex_v),
        mtdna_all: gene_index.resolve(&mtdna_all),
        nuclear_oxphos_complex_i: gene_index.resolve(&geneset.nuclear_oxphos_complex_i),
        nuclear_oxphos_complex_ii: gene_index.resolve(&geneset.nuclear_oxphos_complex_ii),
        nuclear_oxphos_complex_iii: gene_index.resolve(&geneset.nuclear_oxphos_complex_iii),
        nuclear_oxphos_complex_iv: gene_index.resolve(&geneset.nuclear_oxphos_complex_iv),
        nuclear_oxphos_complex_v: gene_index.resolve(&geneset.nuclear_oxphos_complex_v),
        nuclear_oxphos_all: gene_index.resolve(&nuclear_oxphos_all),
        ros: gene_index.resolve(&geneset.ros_detox_genes),
        mitophagy: gene_index.resolve(&geneset.mitophagy_genes),
        fusion: gene_index.resolve(&geneset.dynamics_fusion),
        fission: gene_index.resolve(&geneset.dynamics_fission),
        biogenesis: gene_index.resolve(&geneset.biogenesis_genes),
    }
}
