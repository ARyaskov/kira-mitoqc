//! Weights configuration parsing and validation.

use serde::Deserialize;

use crate::util::numeric::approx_eq;

#[derive(Debug, Clone, Deserialize)]
pub struct Metadata {
    #[allow(dead_code)]
    pub version: String,
    #[allow(dead_code)]
    pub description: String,
    #[allow(dead_code)]
    pub notes: String,
}

/// Weights for a single axis (bioenergetics, ros, dynamics, regulation).
///
/// Field names use the TOML conventions from `weights_v1.toml` — preserve
/// the explicit renames since the casing mixes acronyms (ETC, mtDNA, ATP).
#[derive(Debug, Clone, Deserialize)]
pub struct AxisWeights {
    #[serde(rename = "ETC_stoichiometry_loss")]
    pub etc_stoichiometry_loss: f32,
    #[serde(rename = "mtDNA_expression_uncoupling")]
    pub mtdna_expression_uncoupling: f32,
    #[serde(rename = "ATP_coupling_loss")]
    pub atp_coupling_loss: f32,
}

/// Weights for the ROS axis.
#[derive(Debug, Clone, Deserialize)]
pub struct RosWeights {
    #[serde(rename = "ROS_response_overdrive")]
    pub ros_response_overdrive: f32,
    #[serde(rename = "NADH_imbalance")]
    pub nadh_imbalance: f32,
}

/// Weights for the dynamics axis.
#[derive(Debug, Clone, Deserialize)]
pub struct DynamicsWeights {
    #[serde(rename = "dynamics_imbalance")]
    pub dynamics_imbalance: f32,
    #[serde(rename = "mitophagy_excess")]
    pub mitophagy_excess: f32,
}

/// Weights for the regulation axis.
#[derive(Debug, Clone, Deserialize)]
pub struct RegulationWeights {
    #[serde(rename = "biogenesis_failure")]
    pub biogenesis_failure: f32,
}

/// Axis-specific proxy weights.
#[derive(Debug, Clone, Deserialize)]
pub struct AxisGroupWeights {
    #[serde(rename = "bioenergetics")]
    pub bioenergetics: AxisWeights,
    #[serde(rename = "ros")]
    pub ros: RosWeights,
    #[serde(rename = "dynamics")]
    pub dynamics: DynamicsWeights,
    #[serde(rename = "regulation")]
    pub regulation: RegulationWeights,
}

/// Global axis weights.
#[derive(Debug, Clone, Deserialize)]
pub struct GlobalWeights {
    pub bioenergetics: f32,
    pub ros: f32,
    pub dynamics: f32,
    pub regulation: f32,
}

/// Explainability configuration.
#[derive(Debug, Clone, Deserialize)]
pub struct Explainability {
    pub max_drivers: usize,
    pub min_abs_contribution: f32,
}

/// Weights configuration schema for v1.
#[derive(Debug, Clone, Deserialize)]
pub struct WeightsV1 {
    #[allow(dead_code)]
    pub metadata: Metadata,
    pub axis: AxisGroupWeights,
    pub global: GlobalWeights,
    pub explainability: Explainability,
}

impl WeightsV1 {
    /// Validate axis and global weights sums.
    pub fn validate(&self, tol: f32) -> Result<(), String> {
        let bio_sum = self.axis.bioenergetics.etc_stoichiometry_loss
            + self.axis.bioenergetics.mtdna_expression_uncoupling
            + self.axis.bioenergetics.atp_coupling_loss;
        if !approx_eq(bio_sum, 1.0, tol) {
            return Err(format!(
                "bioenergetics axis weights must sum to 1.0 (got {bio_sum})"
            ));
        }

        let ros_sum = self.axis.ros.ros_response_overdrive + self.axis.ros.nadh_imbalance;
        if !approx_eq(ros_sum, 1.0, tol) {
            return Err(format!("ros axis weights must sum to 1.0 (got {ros_sum})"));
        }

        let dyn_sum = self.axis.dynamics.dynamics_imbalance + self.axis.dynamics.mitophagy_excess;
        if !approx_eq(dyn_sum, 1.0, tol) {
            return Err(format!(
                "dynamics axis weights must sum to 1.0 (got {dyn_sum})"
            ));
        }

        let reg_sum = self.axis.regulation.biogenesis_failure;
        if !approx_eq(reg_sum, 1.0, tol) {
            return Err(format!(
                "regulation axis weights must sum to 1.0 (got {reg_sum})"
            ));
        }

        let global_sum = self.global.bioenergetics
            + self.global.ros
            + self.global.dynamics
            + self.global.regulation;
        if !approx_eq(global_sum, 1.0, tol) {
            return Err(format!(
                "global axis weights must sum to 1.0 (got {global_sum})"
            ));
        }

        Ok(())
    }
}
