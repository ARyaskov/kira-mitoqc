//! v1 proxy computation formulas.

use tracing::warn;

use crate::compute::PrimitiveSignals;
use crate::config::refs::RefsV1;
use crate::core::types::{ProxyKey, ProxyScores};
use crate::input::ResolvedGeneSets;
use crate::proxy::{ProxyError, validate_no_nan};
use crate::util::numeric::{clamp01, safe_div};

/// Compute v1 proxy scores from primitive signals.
pub fn compute_proxies_v1(
    primitives: &PrimitiveSignals,
    resolved: &ResolvedGeneSets,
    refs: &RefsV1,
) -> Result<ProxyScores, ProxyError> {
    enforce_qc(resolved, refs)?;

    let samples = primitives.mtdna_mean.len();

    let mut scores = ProxyScores::default();

    let stoich_raw = primitives.stoich_variance.clone();
    let stoich_norm = stoich_raw
        .iter()
        .map(|v| clamp01(*v / refs.refs.stoich_ref))
        .collect::<Vec<_>>();
    validate_no_nan(ProxyKey::ETCStoichiometryLoss, &stoich_norm)?;
    scores.set_raw(ProxyKey::ETCStoichiometryLoss, stoich_raw);
    scores.set(ProxyKey::ETCStoichiometryLoss, stoich_norm);

    let mut uncoupling_raw = vec![0.0; samples];
    let mut uncoupling_norm = vec![0.0; samples];
    for i in 0..samples {
        let raw = (primitives.mtdna_mean[i] - primitives.nuclear_mean[i]).abs();
        let norm = clamp01(raw / refs.refs.uncoupling_ref);
        uncoupling_raw[i] = raw;
        uncoupling_norm[i] = norm;
    }
    validate_no_nan(ProxyKey::MtdnaExpressionUncoupling, &uncoupling_norm)?;
    scores.set_raw(ProxyKey::MtdnaExpressionUncoupling, uncoupling_raw);
    scores.set(ProxyKey::MtdnaExpressionUncoupling, uncoupling_norm);

    let ros_raw = primitives.ros_mean.clone();
    let ros_norm = ros_raw
        .iter()
        .map(|v| clamp01(*v / refs.refs.ros_ref))
        .collect::<Vec<_>>();
    validate_no_nan(ProxyKey::RosResponseOverdrive, &ros_norm)?;
    scores.set_raw(ProxyKey::RosResponseOverdrive, ros_raw);
    scores.set(ProxyKey::RosResponseOverdrive, ros_norm);

    let mut redox_raw = vec![0.0; samples];
    let mut redox_norm = vec![0.0; samples];
    for i in 0..samples {
        let redox_proxy = safe_div(primitives.c_i[i], primitives.ros_mean[i], refs.eps.value);
        let raw = 1.0 - redox_proxy;
        let norm = clamp01(raw / refs.refs.redox_ref);
        redox_raw[i] = raw;
        redox_norm[i] = norm;
    }
    validate_no_nan(ProxyKey::NadhImbalance, &redox_norm)?;
    scores.set_raw(ProxyKey::NadhImbalance, redox_raw);
    scores.set(ProxyKey::NadhImbalance, redox_norm);

    let mut atp_raw = vec![0.0; samples];
    let mut atp_norm = vec![0.0; samples];
    for i in 0..samples {
        let ratio = safe_div(primitives.atp_mt[i], primitives.atp_nu[i], refs.eps.value);
        let raw = 1.0 - ratio;
        let norm = clamp01(raw / refs.refs.atp_ref);
        atp_raw[i] = raw;
        atp_norm[i] = norm;
    }
    validate_no_nan(ProxyKey::AtpCouplingLoss, &atp_norm)?;
    scores.set_raw(ProxyKey::AtpCouplingLoss, atp_raw);
    scores.set(ProxyKey::AtpCouplingLoss, atp_norm);

    let mitophagy_raw = primitives.mitophagy_mean.clone();
    let mitophagy_norm = mitophagy_raw
        .iter()
        .map(|v| clamp01(*v / refs.refs.mito_ref))
        .collect::<Vec<_>>();
    validate_no_nan(ProxyKey::MitophagyExcess, &mitophagy_norm)?;
    scores.set_raw(ProxyKey::MitophagyExcess, mitophagy_raw);
    scores.set(ProxyKey::MitophagyExcess, mitophagy_norm);

    let mut dynamics_raw = vec![0.0; samples];
    let mut dynamics_norm = vec![0.0; samples];
    for i in 0..samples {
        let ratio = safe_div(
            primitives.fusion_mean[i],
            primitives.fission_mean[i],
            refs.eps.value,
        );
        let raw = ratio.log2().abs();
        let norm = clamp01(raw / refs.refs.dyn_ref);
        dynamics_raw[i] = raw;
        dynamics_norm[i] = norm;
    }
    validate_no_nan(ProxyKey::DynamicsImbalance, &dynamics_norm)?;
    scores.set_raw(ProxyKey::DynamicsImbalance, dynamics_raw);
    scores.set(ProxyKey::DynamicsImbalance, dynamics_norm);

    let mut bio_raw = vec![0.0; samples];
    let mut bio_norm = vec![0.0; samples];
    for i in 0..samples {
        let raw = 1.0 - primitives.biogenesis_mean[i];
        let norm = clamp01(raw / refs.refs.bio_ref);
        bio_raw[i] = raw;
        bio_norm[i] = norm;
    }
    validate_no_nan(ProxyKey::BiogenesisFailure, &bio_norm)?;
    scores.set_raw(ProxyKey::BiogenesisFailure, bio_raw);
    scores.set(ProxyKey::BiogenesisFailure, bio_norm);

    Ok(scores)
}

fn enforce_qc(resolved: &ResolvedGeneSets, refs: &RefsV1) -> Result<(), ProxyError> {
    let mtdna_found = resolved.mtdna_all.found.len();
    if mtdna_found < refs.qc.min_mtdna_genes_found {
        return Err(ProxyError::InsufficientGenes {
            set: "mtDNA",
            found: mtdna_found,
            required: refs.qc.min_mtdna_genes_found,
        });
    }

    let nuc_found = resolved.nuclear_oxphos_all.found.len();
    if nuc_found < refs.qc.min_nuclear_oxphos_found {
        return Err(ProxyError::InsufficientGenes {
            set: "nuclear_oxphos",
            found: nuc_found,
            required: refs.qc.min_nuclear_oxphos_found,
        });
    }

    let ros_found = resolved.ros.found.len();
    if ros_found < refs.qc.min_ros_genes_found {
        warn!(
            ros_found,
            required = refs.qc.min_ros_genes_found,
            "ROS gene count below QC minimum"
        );
    }

    Ok(())
}
