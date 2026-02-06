//! v2 proxy computation with optional multi-omics inputs.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::compute::PrimitiveSignals;
use crate::config::refs_v2::RefsV2;
use crate::core::types::ProxyScores;
use crate::util::numeric::clamp01;

/// v2 proxy keys.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ProxyKeyV2 {
    MtDnaCopyNumberInstability,
    MtDnaHeteroplasmyBurden,
    MtDnaDeletionBurden,
    ProteomicsEtcStoichiometryLoss,
    ProteomicsAtpCouplingLoss,
}

impl ProxyKeyV2 {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::MtDnaCopyNumberInstability => "mtDNA_copy_number_instability",
            Self::MtDnaHeteroplasmyBurden => "mtDNA_heteroplasmy_burden",
            Self::MtDnaDeletionBurden => "mtDNA_deletion_burden",
            Self::ProteomicsEtcStoichiometryLoss => "proteomics_ETC_stoichiometry_loss",
            Self::ProteomicsAtpCouplingLoss => "proteomics_ATP_coupling_loss",
        }
    }
}

impl Serialize for ProxyKeyV2 {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(self.as_str())
    }
}

impl<'de> Deserialize<'de> for ProxyKeyV2 {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        match s.as_str() {
            "mtDNA_copy_number_instability" => Ok(Self::MtDnaCopyNumberInstability),
            "mtDNA_heteroplasmy_burden" => Ok(Self::MtDnaHeteroplasmyBurden),
            "mtDNA_deletion_burden" => Ok(Self::MtDnaDeletionBurden),
            "proteomics_ETC_stoichiometry_loss" => Ok(Self::ProteomicsEtcStoichiometryLoss),
            "proteomics_ATP_coupling_loss" => Ok(Self::ProteomicsAtpCouplingLoss),
            _ => Err(serde::de::Error::custom("unknown v2 proxy key")),
        }
    }
}

/// v2 proxy score bundle.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ProxyScoresV2 {
    pub v1: ProxyScores,
    pub v2_raw: BTreeMap<ProxyKeyV2, Vec<f32>>,
    pub v2_normalized: BTreeMap<ProxyKeyV2, Vec<f32>>,
}

/// Optional multi-omics inputs.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct OptionalOmicsInputs {
    pub mt_dna_copy_number: Option<Vec<f32>>,
    pub heteroplasmy: Option<Vec<f32>>,
    pub mt_dna_deletions: Option<Vec<f32>>,
    pub proteomics_etc: Option<Vec<f32>>,
    pub proteomics_atp: Option<Vec<f32>>,
}

/// Compute v2 proxies, extending v1 proxies with optional omics data.
pub fn compute_proxies_v2(
    primitives: &PrimitiveSignals,
    v1: &ProxyScores,
    refs_v2: &RefsV2,
    extra: &OptionalOmicsInputs,
) -> ProxyScoresV2 {
    let samples = primitives.mtdna_mean.len();

    let mut v2_raw = BTreeMap::new();
    let mut v2_normalized = BTreeMap::new();

    if let Some(values) = extra.mt_dna_copy_number.as_ref() {
        assert_eq!(values.len(), samples);
        let mean = values.iter().sum::<f32>() / samples as f32;
        let raw: Vec<f32> = values.iter().map(|v| (v - mean).abs()).collect();
        let norm: Vec<f32> = raw
            .iter()
            .map(|v| clamp01(*v / refs_v2.refs.cn_ref))
            .collect();
        v2_raw.insert(ProxyKeyV2::MtDnaCopyNumberInstability, raw);
        v2_normalized.insert(ProxyKeyV2::MtDnaCopyNumberInstability, norm);
    }

    if let Some(values) = extra.heteroplasmy.as_ref() {
        assert_eq!(values.len(), samples);
        let raw = values.clone();
        let norm: Vec<f32> = raw
            .iter()
            .map(|v| clamp01(*v / refs_v2.refs.het_ref))
            .collect();
        v2_raw.insert(ProxyKeyV2::MtDnaHeteroplasmyBurden, raw);
        v2_normalized.insert(ProxyKeyV2::MtDnaHeteroplasmyBurden, norm);
    }

    if let Some(values) = extra.mt_dna_deletions.as_ref() {
        assert_eq!(values.len(), samples);
        let raw = values.clone();
        let norm: Vec<f32> = raw
            .iter()
            .map(|v| clamp01(*v / refs_v2.refs.del_ref))
            .collect();
        v2_raw.insert(ProxyKeyV2::MtDnaDeletionBurden, raw);
        v2_normalized.insert(ProxyKeyV2::MtDnaDeletionBurden, norm);
    }

    if let Some(values) = extra.proteomics_etc.as_ref() {
        assert_eq!(values.len(), samples);
        let raw = values.clone();
        let norm: Vec<f32> = raw
            .iter()
            .map(|v| clamp01(*v / refs_v2.refs.prot_stoich_ref))
            .collect();
        v2_raw.insert(ProxyKeyV2::ProteomicsEtcStoichiometryLoss, raw);
        v2_normalized.insert(ProxyKeyV2::ProteomicsEtcStoichiometryLoss, norm);
    }

    if let Some(values) = extra.proteomics_atp.as_ref() {
        assert_eq!(values.len(), samples);
        let raw = values.clone();
        let norm: Vec<f32> = raw
            .iter()
            .map(|v| clamp01(*v / refs_v2.refs.prot_atp_ref))
            .collect();
        v2_raw.insert(ProxyKeyV2::ProteomicsAtpCouplingLoss, raw);
        v2_normalized.insert(ProxyKeyV2::ProteomicsAtpCouplingLoss, norm);
    }

    ProxyScoresV2 {
        v1: v1.clone(),
        v2_raw,
        v2_normalized,
    }
}
