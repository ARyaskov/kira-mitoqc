//! v2 weights configuration.

use std::fs;
use std::path::Path;

use serde::Deserialize;

use crate::config::ConfigError;
use crate::util::numeric::approx_eq;

const EMBEDDED_WEIGHTS_V2: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/assets/weights_v2.toml"
));

#[derive(Debug, Clone, Deserialize)]
pub struct AxisBioenergeticsV2 {
    #[serde(rename = "ETC_stoichiometry_loss")]
    pub etc_stoichiometry_loss: f32,
    #[serde(rename = "mtDNA_expression_uncoupling")]
    pub mtdna_expression_uncoupling: f32,
    #[serde(rename = "ATP_coupling_loss")]
    pub atp_coupling_loss: f32,
    #[serde(rename = "mtDNA_copy_number_instability")]
    pub mtdna_copy_number_instability: f32,
    #[serde(rename = "mtDNA_heteroplasmy_burden")]
    pub mtdna_heteroplasmy_burden: f32,
}

#[derive(Debug, Clone, Deserialize)]
pub struct AxisRegulationV2 {
    #[serde(rename = "biogenesis_failure")]
    pub biogenesis_failure: f32,
    #[serde(rename = "mtDNA_copy_number_instability")]
    pub mtdna_copy_number_instability: f32,
    #[serde(rename = "mtDNA_heteroplasmy_burden")]
    pub mtdna_heteroplasmy_burden: f32,
}

#[derive(Debug, Clone, Deserialize)]
pub struct AxisGroupWeightsV2 {
    pub bioenergetics: AxisBioenergeticsV2,
    pub ros: super::weights::RosWeights,
    pub dynamics: super::weights::DynamicsWeights,
    pub regulation: AxisRegulationV2,
}

#[derive(Debug, Clone, Deserialize)]
pub struct WeightsV2 {
    pub axis: AxisGroupWeightsV2,
    pub global: super::weights::GlobalWeights,
    pub explainability: super::weights::Explainability,
}

impl WeightsV2 {
    pub fn validate(&self, tol: f32) -> Result<(), String> {
        let bio_sum = self.axis.bioenergetics.etc_stoichiometry_loss
            + self.axis.bioenergetics.mtdna_expression_uncoupling
            + self.axis.bioenergetics.atp_coupling_loss
            + self.axis.bioenergetics.mtdna_copy_number_instability
            + self.axis.bioenergetics.mtdna_heteroplasmy_burden;
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

        let reg_sum = self.axis.regulation.biogenesis_failure
            + self.axis.regulation.mtdna_copy_number_instability
            + self.axis.regulation.mtdna_heteroplasmy_burden;
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

pub fn load_weights_v2(path: &Path) -> Result<WeightsV2, ConfigError> {
    let contents = fs::read_to_string(path).map_err(|source| ConfigError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let weights: WeightsV2 = toml::from_str(&contents).map_err(|source| ConfigError::Toml {
        path: path.to_path_buf(),
        source,
    })?;
    weights.validate(1e-6).map_err(ConfigError::Validation)?;
    Ok(weights)
}

pub fn load_weights_v2_embedded() -> Result<WeightsV2, ConfigError> {
    let weights: WeightsV2 =
        toml::from_str(EMBEDDED_WEIGHTS_V2).map_err(|source| ConfigError::Toml {
            path: std::path::PathBuf::from("embedded://weights_v2.toml"),
            source,
        })?;
    weights.validate(1e-6).map_err(ConfigError::Validation)?;
    Ok(weights)
}
