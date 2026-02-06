use std::path::Path;

use kira_mitoqc::config::ConfigV1;

#[test]
fn load_assets_and_validate() {
    let config = ConfigV1::load_from_assets_dir(Path::new("assets"))
        .expect("assets should load and validate");

    let geneset = &config.geneset;
    let refs = &config.refs;

    assert!(
        geneset.all_mtdna().len() >= refs.qc.min_mtdna_genes_found,
        "mtDNA gene count below QC minimum"
    );
    assert!(
        geneset.all_nuclear_oxphos().len() >= refs.qc.min_nuclear_oxphos_found,
        "nuclear OXPHOS gene count below QC minimum"
    );
    assert!(
        geneset.ros_detox_genes.len() >= refs.qc.min_ros_genes_found,
        "ROS detox gene count below QC minimum"
    );

    let weights = &config.weights;
    let bio_sum = weights.axis.bioenergetics.etc_stoichiometry_loss
        + weights.axis.bioenergetics.mtdna_expression_uncoupling
        + weights.axis.bioenergetics.atp_coupling_loss;
    let ros_sum = weights.axis.ros.ros_response_overdrive + weights.axis.ros.nadh_imbalance;
    let dyn_sum = weights.axis.dynamics.dynamics_imbalance + weights.axis.dynamics.mitophagy_excess;
    let reg_sum = weights.axis.regulation.biogenesis_failure;
    let global_sum = weights.global.bioenergetics
        + weights.global.ros
        + weights.global.dynamics
        + weights.global.regulation;

    let tol = 1e-6;
    assert!((bio_sum - 1.0).abs() <= tol, "bioenergetics sum != 1");
    assert!((ros_sum - 1.0).abs() <= tol, "ros sum != 1");
    assert!((dyn_sum - 1.0).abs() <= tol, "dynamics sum != 1");
    assert!((reg_sum - 1.0).abs() <= tol, "regulation sum != 1");
    assert!((global_sum - 1.0).abs() <= tol, "global sum != 1");
}
