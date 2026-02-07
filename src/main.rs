use std::path::PathBuf;

use clap::Parser;
use kira_mitoqc::cache::{
    ExprCacheMode, mmap_expr_bin, mmap_organelle_bin, write_expr_bin_with_mode,
};
use kira_mitoqc::classify::classify_v1;
use kira_mitoqc::compute::compute_primitives;
use kira_mitoqc::config::ConfigV1;
use kira_mitoqc::config::refs_v2::load_refs_v2;
use kira_mitoqc::config::weights_v2::load_weights_v2;
use kira_mitoqc::data::{AggregationMode, load_cluster_map, prepare_expression_with_clusters};
use kira_mitoqc::explain::explain_v1;
use kira_mitoqc::input::{
    GeneIndex, GeneResolutionQC, InputMode, InputSpec, resolve_all_genesets, validate_input_path,
};
use kira_mitoqc::io::mtx::{
    discover_dataset_files, load_mtx_dir, load_mtx_metadata, resolve_shared_cache_filename,
};
use kira_mitoqc::output::v2::assemble_profiles_v2;
use kira_mitoqc::output::{
    pipeline_contract::{write_mito_metrics_tsv, write_pipeline_step_json, write_summary_json},
    profile::assemble_profiles_v1,
    write_axes_tsv, write_decay_tsv, write_json, write_json_v2, write_proxies_tsv,
};
use kira_mitoqc::proxy::{OptionalOmicsInputs, compute_proxies_v1, compute_proxies_v2};
use kira_mitoqc::score::{compute_axes_v2, compute_decay_v2, score_profile_v1};
use tracing::info;
use tracing_subscriber::EnvFilter;

#[derive(Parser, Debug)]
#[command(
    name = "kira-mitoqc",
    version,
    about = "Deterministic mitochondrial QC scoring"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(clap::Subcommand, Debug)]
enum Commands {
    /// Run the pipeline on an input dataset (parsing not yet implemented).
    Run(RunArgs),
}

#[derive(Parser, Debug)]
struct RunArgs {
    /// Input path (unused in stage 01).
    #[arg(long)]
    input: PathBuf,

    /// Output directory (unused in stage 01).
    #[arg(long)]
    out: PathBuf,

    /// Aggregation mode (parsed but unused in stage 01).
    #[arg(long, value_enum, default_value = "sample")]
    mode: InputMode,

    /// Assets directory containing geneset/weights/refs TOML.
    #[arg(long, default_value = "assets")]
    assets: PathBuf,

    /// Cluster assignment file (TSV: barcode<TAB>cluster).
    #[arg(long)]
    clusters: Option<PathBuf>,

    /// Cache directory for expression binary.
    #[arg(long)]
    cache: Option<PathBuf>,

    /// Pipeline version.
    #[arg(long, value_enum, default_value = "v1")]
    version: RunVersion,

    /// Optional mtDNA copy number vector (one value per line).
    #[arg(long)]
    mtcopy: Option<PathBuf>,

    /// Optional heteroplasmy vector (one value per line).
    #[arg(long)]
    heteroplasmy: Option<PathBuf>,

    /// Optional mtDNA deletions vector (one value per line).
    #[arg(long)]
    mtdeletions: Option<PathBuf>,

    /// Optional proteomics ETC stoichiometry vector (one value per line).
    #[arg(long)]
    proteomics_etc: Option<PathBuf>,

    /// Optional proteomics ATP coupling vector (one value per line).
    #[arg(long)]
    proteomics_atp: Option<PathBuf>,

    /// H5AD gene symbol column in var (default: auto).
    #[arg(long)]
    gene_symbol_key: Option<String>,

    /// Gene symbol column for MTX features.tsv (1-based).
    #[arg(long)]
    gene_symbol_col: Option<usize>,

    /// Execution mode.
    #[arg(long, value_enum, default_value = "standalone")]
    run_mode: RunMode,
}

#[derive(clap::ValueEnum, Debug, Copy, Clone, PartialEq, Eq)]
enum RunVersion {
    V1,
    V2,
}

#[derive(clap::ValueEnum, Debug, Copy, Clone, PartialEq, Eq)]
enum RunMode {
    Standalone,
    Pipeline,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let filter = EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info"));
    tracing_subscriber::fmt().with_env_filter(filter).init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Run(args) => {
            validate_input_path(&args.input)?;
            let input_spec = InputSpec {
                input_path: args.input.clone(),
                mode: args.mode,
            };

            let cache_dir = args.cache.as_ref().unwrap_or(&args.out);
            let cache_path = cache_dir.join("expr.bin");
            let requested_cache_mode = cache_mode_from_input_mode(args.mode);

            let gene_symbol_col = args.gene_symbol_col.map(|v| v.saturating_sub(1));

            let is_h5ad = args
                .input
                .extension()
                .and_then(|ext| ext.to_str())
                .map(|ext| ext.eq_ignore_ascii_case("h5ad"))
                .unwrap_or(false);

            let mut cached_result: Option<(
                Vec<String>,
                Vec<String>,
                kira_mitoqc::cache::ExpressionSoAView<'static>,
                ConfigV1,
            )> = None;
            if cache_path.exists() {
                let (features, barcodes) = if is_h5ad {
                    load_h5ad_metadata_guarded(&args.input, args.gene_symbol_key.as_deref())?
                } else {
                    load_mtx_metadata(&args.input, gene_symbol_col)?
                };
                let config = load_config_autodetect(&args.assets, &features)?;
                let view = mmap_expr_bin(&cache_path)?;
                if view.mode == requested_cache_mode {
                    info!(path = ?cache_path, mode = ?view.mode, "using cached expression");
                    info!(
                        genes = view.genes,
                        samples = view.samples,
                        "Mapped expression cache"
                    );
                    cached_result = Some((features, barcodes, view, config));
                } else {
                    info!(
                        path = ?cache_path,
                        cached_mode = ?view.mode,
                        requested_mode = ?requested_cache_mode,
                        "cache mode mismatch, rebuilding expression cache"
                    );
                }
            }

            let (features, barcodes, view, config) = if let Some(cached) = cached_result {
                cached
            } else {
                info!(path = ?cache_path, "building expression cache");
                let (matrix, features, barcodes) = if is_h5ad {
                    if args.run_mode == RunMode::Pipeline {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidInput,
                            "pipeline run-mode is only supported for MTX inputs",
                        )
                        .into());
                    }
                    let (matrix, features, barcodes) =
                        load_h5ad_guarded(&args.input, args.gene_symbol_key.as_deref())?;
                    let gene_source = args
                        .gene_symbol_key
                        .as_deref()
                        .map(|k| format!("var/{k}"))
                        .unwrap_or_else(|| "var/gene_symbols or var/_index".to_string());
                    let sparse_format = if matrix.is_csc() {
                        "csc"
                    } else if matrix.is_csr() {
                        "csr"
                    } else {
                        "unknown"
                    };
                    info!(
                        features = features.len(),
                        barcodes = barcodes.len(),
                        format = sparse_format,
                        gene_symbols = gene_source.as_str(),
                        "Loaded H5AD input"
                    );
                    (matrix, features, barcodes)
                } else {
                    if args.run_mode == RunMode::Pipeline {
                        let discovery = discover_dataset_files(&args.input)?;
                        let shared_name =
                            resolve_shared_cache_filename(discovery.prefix.as_deref());
                        let shared_path = discovery.input_dir.join(shared_name);
                        info!(path = ?args.input, "Pipeline mode: loading MTX input for shared cache");
                        let mtx = load_mtx_dir(&args.input, gene_symbol_col)?;
                        info!(
                            rows = mtx.matrix.rows(),
                            cols = mtx.matrix.cols(),
                            nnz = mtx.matrix.nnz(),
                            "Pipeline mode: MTX loaded"
                        );
                        info!(path = ?shared_path, "Pipeline mode: writing shared organelle cache");
                        kira_mitoqc::cache::write_organelle_bin(&shared_path, &mtx)?;
                        info!(path = ?shared_path, "Wrote shared pipeline cache");
                        info!(path = ?shared_path, "Pipeline mode: opening shared cache via mmap");
                        let shared = mmap_organelle_bin(&shared_path)?;
                        info!(
                            path = ?shared_path,
                            genes = shared.n_genes,
                            cells = shared.n_cells,
                            nnz = shared.nnz,
                            "Using mmap shared cache for MTX input"
                        );
                        let mtx = shared.to_mtx_input();
                        (mtx.matrix, mtx.features, mtx.barcodes)
                    } else {
                        let mtx = load_mtx_dir(&args.input, gene_symbol_col)?;
                        info!(
                            features = mtx.features.len(),
                            barcodes = mtx.barcodes.len(),
                            "Loaded MTX input"
                        );
                        (mtx.matrix, mtx.features, mtx.barcodes)
                    }
                };

                let config = load_config_autodetect(&args.assets, &features)?;
                let geneset = &config.geneset;

                let gene_index = GeneIndex::try_from_feature_list(&features)?;
                let resolved = resolve_all_genesets(&gene_index, geneset);

                let mtdna_qc = GeneResolutionQC::from_resolution(&resolved.mtdna_all);
                let nuc_qc = GeneResolutionQC::from_resolution(&resolved.nuclear_oxphos_all);
                let ros_qc = GeneResolutionQC::from_resolution(&resolved.ros);
                let mito_qc = GeneResolutionQC::from_resolution(&resolved.mitophagy);
                let fusion_qc = GeneResolutionQC::from_resolution(&resolved.fusion);
                let fission_qc = GeneResolutionQC::from_resolution(&resolved.fission);
                let bio_qc = GeneResolutionQC::from_resolution(&resolved.biogenesis);

                info!(
                    mtdna_genes = geneset.all_mtdna().len(),
                    nuclear_oxphos_genes = geneset.all_nuclear_oxphos().len(),
                    ros_genes = geneset.ros_detox_genes.len(),
                    mitophagy_genes = geneset.mitophagy_genes.len(),
                    fusion_genes = geneset.dynamics_fusion.len(),
                    fission_genes = geneset.dynamics_fission.len(),
                    biogenesis_genes = geneset.biogenesis_genes.len(),
                    "Loaded config and gene sets"
                );

                info!(mode = ?input_spec.mode, "Selected input mode");

                info!(
                    mtdna_found = mtdna_qc.found,
                    mtdna_missing = mtdna_qc.missing,
                    nuclear_found = nuc_qc.found,
                    nuclear_missing = nuc_qc.missing,
                    ros_found = ros_qc.found,
                    ros_missing = ros_qc.missing,
                    mitophagy_found = mito_qc.found,
                    mitophagy_missing = mito_qc.missing,
                    fusion_found = fusion_qc.found,
                    fusion_missing = fusion_qc.missing,
                    fission_found = fission_qc.found,
                    fission_missing = fission_qc.missing,
                    biogenesis_found = bio_qc.found,
                    biogenesis_missing = bio_qc.missing,
                    "Resolved gene sets"
                );

                let agg_mode = match input_spec.mode {
                    InputMode::Sample => AggregationMode::Sample,
                    InputMode::Cluster => AggregationMode::Cluster,
                    InputMode::Cell => AggregationMode::Cell,
                };

                let cluster_map = if agg_mode == AggregationMode::Cluster {
                    let path = args.clusters.as_ref().ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidInput,
                            "--clusters is required for cluster mode",
                        )
                    })?;
                    if is_h5ad {
                        let column = parse_h5ad_cluster_arg(path)?;
                        Some(load_h5ad_clusters_guarded(&args.input, &column)?)
                    } else {
                        Some(load_cluster_map(path, &barcodes)?)
                    }
                } else {
                    None
                };

                let prepared = prepare_expression_with_clusters(
                    &kira_mitoqc::io::mtx::MtxInput {
                        matrix,
                        features: features.clone(),
                        barcodes: barcodes.clone(),
                    },
                    &resolved,
                    agg_mode,
                    cluster_map.as_ref(),
                )?;

                info!(
                    genes = prepared.soa.genes,
                    samples = prepared.soa.samples,
                    "Prepared expression SoA"
                );

                std::fs::create_dir_all(cache_dir)?;
                write_expr_bin_with_mode(&cache_path, &prepared.soa, requested_cache_mode)?;
                info!(
                    path = ?cache_path,
                    mode = ?requested_cache_mode,
                    "Wrote expression cache"
                );

                let view = mmap_expr_bin(&cache_path)?;
                (features, barcodes, view, config)
            };
            let geneset = &config.geneset;

            let gene_index = GeneIndex::try_from_feature_list(&features)?;
            let resolved = resolve_all_genesets(&gene_index, geneset);

            let mtdna_qc = GeneResolutionQC::from_resolution(&resolved.mtdna_all);
            let nuc_qc = GeneResolutionQC::from_resolution(&resolved.nuclear_oxphos_all);
            let ros_qc = GeneResolutionQC::from_resolution(&resolved.ros);
            let mito_qc = GeneResolutionQC::from_resolution(&resolved.mitophagy);
            let fusion_qc = GeneResolutionQC::from_resolution(&resolved.fusion);
            let fission_qc = GeneResolutionQC::from_resolution(&resolved.fission);
            let bio_qc = GeneResolutionQC::from_resolution(&resolved.biogenesis);

            info!(
                mtdna_found = mtdna_qc.found,
                mtdna_missing = mtdna_qc.missing,
                nuclear_found = nuc_qc.found,
                nuclear_missing = nuc_qc.missing,
                ros_found = ros_qc.found,
                ros_missing = ros_qc.missing,
                mitophagy_found = mito_qc.found,
                mitophagy_missing = mito_qc.missing,
                fusion_found = fusion_qc.found,
                fusion_missing = fusion_qc.missing,
                fission_found = fission_qc.found,
                fission_missing = fission_qc.missing,
                biogenesis_found = bio_qc.found,
                biogenesis_missing = bio_qc.missing,
                "Resolved gene sets"
            );

            info!("Computing primitives");
            let primitives = compute_primitives(&view, &resolved);

            info!("Computing proxies");
            let proxies = compute_proxies_v1(&primitives, &resolved, &config.refs)?;

            info!("Computing axes and decay");
            let scored = score_profile_v1(&proxies, &config.weights);

            info!("Classifying failure modes");
            let states = classify_v1(&scored.axes, &scored.decay, &config.refs);

            info!("Computing explainability");
            let explain = explain_v1(
                &proxies,
                &scored.axes,
                &scored.decay,
                &states,
                &config.weights,
            );

            info!("Assembling profiles");
            let profiles =
                assemble_profiles_v1(&states, &scored.decay, &scored.axes, &proxies, &explain);

            info!(path = ?args.out, "Writing outputs");
            write_json(&args.out, &profiles)?;
            write_axes_tsv(&args.out, &barcodes, &scored.axes)?;
            write_decay_tsv(&args.out, &scored.decay)?;
            write_proxies_tsv(&args.out, &proxies)?;

            if args.run_mode == RunMode::Pipeline {
                let shared_cache_name = if is_h5ad {
                    "kira-organelle.bin".to_string()
                } else {
                    discover_dataset_files(&args.input)
                        .map(|d| resolve_shared_cache_filename(d.prefix.as_deref()).to_string())
                        .unwrap_or_else(|_| "kira-organelle.bin".to_string())
                };

                write_summary_json(&args.out, &profiles)?;
                write_mito_metrics_tsv(&args.out, &barcodes, &profiles)?;
                write_pipeline_step_json(&args.out, &shared_cache_name)?;
            }

            if args.version == RunVersion::V2 {
                info!("Computing v2 proxies and axes");
                let refs_v2 = load_refs_v2(&args.assets.join("refs_v2.toml"))?;
                let weights_v2 = load_weights_v2(&args.assets.join("weights_v2.toml"))?;

                let omics = OptionalOmicsInputs {
                    mt_dna_copy_number: read_vec_file(args.mtcopy.as_ref())?,
                    heteroplasmy: read_vec_file(args.heteroplasmy.as_ref())?,
                    mt_dna_deletions: read_vec_file(args.mtdeletions.as_ref())?,
                    proteomics_etc: read_vec_file(args.proteomics_etc.as_ref())?,
                    proteomics_atp: read_vec_file(args.proteomics_atp.as_ref())?,
                };

                let proxies_v2 = compute_proxies_v2(&primitives, &proxies, &refs_v2, &omics);
                let axes_v2 = compute_axes_v2(&proxies_v2, &weights_v2, &refs_v2);
                let decay_v2 = compute_decay_v2(&axes_v2, &weights_v2);
                let v2_profiles = assemble_profiles_v2(&profiles, &axes_v2, &decay_v2, &proxies_v2);
                write_json_v2(&args.out, &v2_profiles)?;
            }

            let _ = (input_spec, barcodes);
        }
    }

    Ok(())
}

fn load_config_autodetect(
    assets_dir: &std::path::Path,
    features: &[String],
) -> Result<ConfigV1, Box<dyn std::error::Error>> {
    let human = ConfigV1::load_from_assets_dir_with_geneset(assets_dir, "geneset_v1.toml")?;
    let mouse = ConfigV1::load_from_assets_dir_with_geneset(assets_dir, "geneset_mouse_v1.toml")?;

    let human_hits = geneset_overlap_hits(features, &human.geneset);
    let mouse_hits = geneset_overlap_hits(features, &mouse.geneset);

    if mouse_hits > human_hits {
        info!(
            geneset = "geneset_mouse_v1.toml",
            human_hits, mouse_hits, "Selected geneset"
        );
        Ok(mouse)
    } else {
        info!(
            geneset = "geneset_v1.toml",
            human_hits, mouse_hits, "Selected geneset"
        );
        Ok(human)
    }
}

fn cache_mode_from_input_mode(mode: InputMode) -> ExprCacheMode {
    match mode {
        InputMode::Sample => ExprCacheMode::Sample,
        InputMode::Cluster => ExprCacheMode::Cluster,
        InputMode::Cell => ExprCacheMode::Cell,
    }
}

fn geneset_overlap_hits(features: &[String], geneset: &kira_mitoqc::core::types::GeneSet) -> usize {
    let feature_set: std::collections::BTreeSet<&str> =
        features.iter().map(std::string::String::as_str).collect();

    let mut hit_count = 0usize;
    for gene in geneset.all_mtdna() {
        if feature_set.contains(gene) {
            hit_count += 1;
        }
    }
    for gene in geneset.all_nuclear_oxphos() {
        if feature_set.contains(gene) {
            hit_count += 1;
        }
    }
    for gene in &geneset.ros_detox_genes {
        if feature_set.contains(gene.as_str()) {
            hit_count += 1;
        }
    }
    for gene in &geneset.mitophagy_genes {
        if feature_set.contains(gene.as_str()) {
            hit_count += 1;
        }
    }
    for gene in &geneset.dynamics_fusion {
        if feature_set.contains(gene.as_str()) {
            hit_count += 1;
        }
    }
    for gene in &geneset.dynamics_fission {
        if feature_set.contains(gene.as_str()) {
            hit_count += 1;
        }
    }
    for gene in &geneset.biogenesis_genes {
        if feature_set.contains(gene.as_str()) {
            hit_count += 1;
        }
    }
    hit_count
}

fn read_vec_file(path: Option<&PathBuf>) -> Result<Option<Vec<f32>>, Box<dyn std::error::Error>> {
    let Some(path) = path else {
        return Ok(None);
    };
    let contents = std::fs::read_to_string(path)?;
    let mut values = Vec::new();
    for (line_no, line) in contents.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let value: f32 = trimmed.parse().map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("invalid float at line {} in {:?}", line_no + 1, path),
            )
        })?;
        values.push(value);
    }
    Ok(Some(values))
}

fn parse_h5ad_cluster_arg(path: &PathBuf) -> Result<String, Box<dyn std::error::Error>> {
    let s = path.to_string_lossy();
    if let Some(rest) = s.strip_prefix("obs:") {
        if rest.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "missing obs column name for clusters",
            )
            .into());
        }
        Ok(rest.to_string())
    } else {
        Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "for h5ad, use --clusters obs:<column>",
        )
        .into())
    }
}

#[cfg(feature = "h5ad")]
fn load_h5ad_guarded(
    path: &PathBuf,
    gene_symbol_key: Option<&str>,
) -> Result<(sprs::CsMat<f32>, Vec<String>, Vec<String>), Box<dyn std::error::Error>> {
    let h5ad = kira_mitoqc::io::h5ad::load_h5ad(path, gene_symbol_key)?;
    Ok((h5ad.matrix, h5ad.features, h5ad.barcodes))
}

#[cfg(not(feature = "h5ad"))]
fn load_h5ad_guarded(
    path: &PathBuf,
    _gene_symbol_key: Option<&str>,
) -> Result<(sprs::CsMat<f32>, Vec<String>, Vec<String>), Box<dyn std::error::Error>> {
    Err(kira_mitoqc::input::InputError::MissingFeature {
        feature: "h5ad",
        path: path.clone(),
    }
    .into())
}

#[cfg(feature = "h5ad")]
fn load_h5ad_metadata_guarded(
    path: &PathBuf,
    gene_symbol_key: Option<&str>,
) -> Result<(Vec<String>, Vec<String>), Box<dyn std::error::Error>> {
    Ok(kira_mitoqc::io::h5ad::load_h5ad_metadata(
        path,
        gene_symbol_key,
    )?)
}

#[cfg(not(feature = "h5ad"))]
fn load_h5ad_metadata_guarded(
    path: &PathBuf,
    _gene_symbol_key: Option<&str>,
) -> Result<(Vec<String>, Vec<String>), Box<dyn std::error::Error>> {
    Err(kira_mitoqc::input::InputError::MissingFeature {
        feature: "h5ad",
        path: path.clone(),
    }
    .into())
}

#[cfg(feature = "h5ad")]
fn load_h5ad_clusters_guarded(
    path: &PathBuf,
    column: &str,
) -> Result<kira_mitoqc::data::ClusterMap, Box<dyn std::error::Error>> {
    Ok(kira_mitoqc::io::h5ad::load_h5ad_clusters(path, column)?)
}

#[cfg(not(feature = "h5ad"))]
fn load_h5ad_clusters_guarded(
    path: &PathBuf,
    _column: &str,
) -> Result<kira_mitoqc::data::ClusterMap, Box<dyn std::error::Error>> {
    Err(kira_mitoqc::input::InputError::MissingFeature {
        feature: "h5ad",
        path: path.clone(),
    }
    .into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn run_mode_default_is_standalone() {
        let cli = Cli::parse_from(["kira-mitoqc", "run", "--input", "in", "--out", "out"]);
        let Commands::Run(args) = cli.command;
        assert_eq!(args.run_mode, RunMode::Standalone);
    }

    #[test]
    fn run_mode_pipeline_is_accepted() {
        let cli = Cli::parse_from([
            "kira-mitoqc",
            "run",
            "--input",
            "in",
            "--out",
            "out",
            "--run-mode",
            "pipeline",
        ]);
        let Commands::Run(args) = cli.command;
        assert_eq!(args.run_mode, RunMode::Pipeline);
    }

    #[test]
    fn clap_schema_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn geneset_overlap_prefers_mouse_symbols() {
        let human = ConfigV1::load_from_assets_dir_with_geneset(
            std::path::Path::new("assets"),
            "geneset_v1.toml",
        )
        .expect("human config");
        let mouse = ConfigV1::load_from_assets_dir_with_geneset(
            std::path::Path::new("assets"),
            "geneset_mouse_v1.toml",
        )
        .expect("mouse config");
        let features = vec![
            "mt-Nd1".to_string(),
            "Atp5f1a".to_string(),
            "Sod2".to_string(),
        ];
        let h = geneset_overlap_hits(&features, &human.geneset);
        let m = geneset_overlap_hits(&features, &mouse.geneset);
        assert!(m > h);
    }
}
