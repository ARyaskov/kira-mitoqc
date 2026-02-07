//! Numeric compute kernels for primitive signals.

use crate::cache::ExpressionSoAView;
use crate::input::ResolvedGeneSets;

pub mod avx2;
pub mod dispatch;
pub mod neon;
pub mod scalar;

/// Tolerance for SIMD vs scalar comparisons.
pub const SIMD_EPS: f32 = 1e-5;

/// Primitive per-sample signals derived from expression values.
#[derive(Debug, Clone, PartialEq)]
pub struct PrimitiveSignals {
    pub mtdna_mean: Vec<f32>,
    pub nuclear_mean: Vec<f32>,

    pub c_i: Vec<f32>,
    pub c_iii: Vec<f32>,
    pub c_iv: Vec<f32>,
    pub c_v: Vec<f32>,

    pub ros_mean: Vec<f32>,
    pub mitophagy_mean: Vec<f32>,
    pub fusion_mean: Vec<f32>,
    pub fission_mean: Vec<f32>,
    pub biogenesis_mean: Vec<f32>,

    pub atp_mt: Vec<f32>,
    pub atp_nu: Vec<f32>,

    pub stoich_variance: Vec<f32>,
}

/// Precomputed offsets into the SoA layout.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneOffsets {
    pub mtdna_all: Vec<usize>,
    pub nuclear_oxphos: Vec<usize>,

    pub complex_i: Vec<usize>,
    pub complex_iii: Vec<usize>,
    pub complex_iv: Vec<usize>,
    pub complex_v: Vec<usize>,

    pub ros: Vec<usize>,
    pub mitophagy: Vec<usize>,
    pub fusion: Vec<usize>,
    pub fission: Vec<usize>,
    pub biogenesis: Vec<usize>,

    pub atp_mt: Vec<usize>,
    pub atp_nu: usize,
}

impl GeneOffsets {
    /// Build offsets based on the canonical SoA ordering.
    pub fn from_resolved(resolved: &ResolvedGeneSets) -> Self {
        let mut index_map = std::collections::BTreeMap::new();
        let mut cursor = 0usize;

        let complex_i = push_list(&mut cursor, &mut index_map, &resolved.mtdna_complex_i.genes);
        let complex_iii = push_list(
            &mut cursor,
            &mut index_map,
            &resolved.mtdna_complex_iii.genes,
        );
        let complex_iv = push_list(
            &mut cursor,
            &mut index_map,
            &resolved.mtdna_complex_iv.genes,
        );
        let complex_v = push_list(&mut cursor, &mut index_map, &resolved.mtdna_complex_v.genes);

        let nuclear_oxphos = concat_lists(&[
            &resolved.nuclear_oxphos_complex_i.genes,
            &resolved.nuclear_oxphos_complex_ii.genes,
            &resolved.nuclear_oxphos_complex_iii.genes,
            &resolved.nuclear_oxphos_complex_iv.genes,
            &resolved.nuclear_oxphos_complex_v.genes,
        ]);
        let nuclear_offsets = push_list(&mut cursor, &mut index_map, &nuclear_oxphos);

        let ros = push_list(&mut cursor, &mut index_map, &resolved.ros.genes);
        let mitophagy = push_list(&mut cursor, &mut index_map, &resolved.mitophagy.genes);
        let fusion = push_list(&mut cursor, &mut index_map, &resolved.fusion.genes);
        let fission = push_list(&mut cursor, &mut index_map, &resolved.fission.genes);
        let biogenesis = push_list(&mut cursor, &mut index_map, &resolved.biogenesis.genes);

        let mut mtdna_all = Vec::new();
        mtdna_all.extend(complex_i.iter().copied());
        mtdna_all.extend(complex_iii.iter().copied());
        mtdna_all.extend(complex_iv.iter().copied());
        mtdna_all.extend(complex_v.iter().copied());

        let atp_mt: Vec<usize> = complex_v.iter().copied().take(2).collect();
        let atp_nu = nuclear_offsets.first().copied().unwrap_or(usize::MAX);

        Self {
            mtdna_all,
            nuclear_oxphos: nuclear_offsets,
            complex_i,
            complex_iii,
            complex_iv,
            complex_v,
            ros,
            mitophagy,
            fusion,
            fission,
            biogenesis,
            atp_mt,
            atp_nu,
        }
    }
}

fn concat_lists(lists: &[&[String]]) -> Vec<String> {
    let mut merged = Vec::new();
    for list in lists {
        merged.extend(list.iter().cloned());
    }
    merged
}

fn push_list(
    cursor: &mut usize,
    index_map: &mut std::collections::BTreeMap<String, usize>,
    list: &[String],
) -> Vec<usize> {
    let mut offsets = Vec::with_capacity(list.len());
    for gene in list {
        offsets.push(*cursor);
        index_map.insert(gene.clone(), *cursor);
        *cursor += 1;
    }
    offsets
}

/// Compute primitive signals with SIMD dispatch.
pub fn compute_primitives(
    soa: &ExpressionSoAView,
    resolved: &ResolvedGeneSets,
) -> PrimitiveSignals {
    dispatch::compute_primitives(soa, resolved)
}
