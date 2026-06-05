use std::path::PathBuf;

use clap::Parser;
use kira_mitoqc::cache::{
    ExprCacheMode, mmap_expr_bin, mmap_organelle_bin, write_expr_bin_with_mode,
};
use kira_mitoqc::classify::classify_v1_with_redox;
use kira_mitoqc::compute::compute_primitives;
use kira_mitoqc::config::ConfigV1;
use kira_mitoqc::config::refs_v2::{load_refs_v2, load_refs_v2_embedded};
use kira_mitoqc::config::weights_v2::{load_weights_v2, load_weights_v2_embedded};
use kira_mitoqc::data::{
    AggregationMode, SoaIndex, load_cluster_map, prepare_expression_with_clusters,
};
use kira_mitoqc::explain::explain_v1;
use kira_mitoqc::input::bd_rhapsody::{
    compute_mito_fraction_from_file, load_bd_rhapsody, load_bd_rhapsody_metadata,
    resolve_bd_input_path,
};
use kira_mitoqc::input::{
    DetectedInputFormat, ExpressionSource, GeneIndex, GeneResolutionQC, InputFormat, InputMode,
    InputSpec, detect_input_format, resolve_all_genesets, validate_input_path,
};
use kira_mitoqc::io::mtx::{
    discover_dataset_files, load_mtx_dir, load_mtx_metadata, resolve_shared_cache_filename,
};
use kira_mitoqc::metrics::metabolic_extension::aggregate::build_summary as build_metabolic_summary;
use kira_mitoqc::metrics::metabolic_extension::panels::{
    BIOGENESIS_PANEL, FAO_PANEL, GLYCOLYSIS_PANEL, OXPHOS_PANEL, ROS_PANEL,
};
use kira_mitoqc::metrics::metabolic_extension::scores::compute_metabolic_metrics;
use kira_mitoqc::output::v2::assemble_profiles_v2;
use kira_mitoqc::output::{
    pipeline_contract::{
        write_mito_metrics_tsv, write_pipeline_step_json, write_summary_json_with_redox,
    },
    profile::assemble_profiles_v1,
    write_axes_tsv, write_decay_tsv, write_json, write_json_v2, write_proxies_tsv,
    write_redox_metrics_tsv,
};
use kira_mitoqc::proxy::{OptionalOmicsInputs, compute_proxies_v1, compute_proxies_v2};
use kira_mitoqc::redox::compute_redox_metrics;
use kira_mitoqc::score::{compute_axes_v2, compute_decay_v2, score_profile_v1};
use tracing::info;
use tracing_subscriber::EnvFilter;

#[derive(Parser, Debug)]
#[command(
    name = "kira-mitoqc",
    version,
    about = "Deterministic mitochondrial QC scoring"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(clap::Subcommand, Debug)]
pub enum Commands {
    /// Run the pipeline on an input dataset (parsing not yet implemented).
    Run(RunArgs),
}

#[derive(Parser, Debug)]
pub struct RunArgs {
    /// Input path (unused in stage 01).
    #[arg(long)]
    pub input: PathBuf,

    /// Output directory (unused in stage 01).
    #[arg(long)]
    pub out: PathBuf,

    /// Aggregation mode (parsed but unused in stage 01).
    #[arg(long, value_enum, default_value = "sample")]
    pub mode: InputMode,

    /// Assets directory containing geneset/weights/refs TOML.
    #[arg(long, default_value = "assets")]
    pub assets: PathBuf,

    /// Cluster assignment file (TSV: barcode<TAB>cluster).
    #[arg(long)]
    pub clusters: Option<PathBuf>,

    /// Cache directory for expression binary.
    #[arg(long)]
    pub cache: Option<PathBuf>,

    /// Pipeline version.
    #[arg(long, value_enum, default_value = "v1")]
    pub version: RunVersion,

    /// Optional mtDNA copy number vector (one value per line).
    #[arg(long)]
    pub mtcopy: Option<PathBuf>,

    /// Optional heteroplasmy vector (one value per line).
    #[arg(long)]
    pub heteroplasmy: Option<PathBuf>,

    /// Optional mtDNA deletions vector (one value per line).
    #[arg(long)]
    pub mtdeletions: Option<PathBuf>,

    /// Optional proteomics ETC stoichiometry vector (one value per line).
    #[arg(long)]
    pub proteomics_etc: Option<PathBuf>,

    /// Optional proteomics ATP coupling vector (one value per line).
    #[arg(long)]
    pub proteomics_atp: Option<PathBuf>,

    /// H5AD gene symbol column in var (default: auto).
    #[arg(long)]
    pub gene_symbol_key: Option<String>,

    /// Gene symbol column for MTX features.tsv (1-based).
    #[arg(long)]
    pub gene_symbol_col: Option<usize>,

    /// Execution mode.
    #[arg(long, value_enum, default_value = "standalone")]
    pub run_mode: RunMode,

    /// Input format selector.
    #[arg(long, value_enum, default_value = "auto")]
    pub input_format: InputFormat,

    // Enable additive redox-proxy extension stage.
    #[arg(long, default_value_t = false)]
    pub redox: bool,
}

#[derive(clap::ValueEnum, Debug, Copy, Clone, PartialEq, Eq)]
pub enum RunVersion {
    V1,
    V2,
}

#[derive(clap::ValueEnum, Debug, Copy, Clone, PartialEq, Eq)]
pub enum RunMode {
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
            let agg_mode = match input_spec.mode {
                InputMode::Sample => AggregationMode::Sample,
                InputMode::Cluster => AggregationMode::Cluster,
                InputMode::Cell => AggregationMode::Cell,
            };

            let gene_symbol_col = args.gene_symbol_col.map(|v| v.saturating_sub(1));

            // Recognize .h5ad (and .h5ad.gz). Falling back on extension keeps
            // the CLI usable when the file is missing/inaccessible at this
            // point; deeper detection happens inside kira-scio.
            let is_h5ad = args
                .input
                .file_name()
                .and_then(|f| f.to_str())
                .map(|name| {
                    let lower = name.to_ascii_lowercase();
                    lower.ends_with(".h5ad") || lower.ends_with(".h5ad.gz")
                })
                .unwrap_or(false);

            let detected_input = if is_h5ad {
                None
            } else {
                Some(detect_input_format(&args.input, args.input_format)?)
            };
            let bd_input_path =
                if matches!(detected_input, Some(DetectedInputFormat::BDRhapsodyDense)) {
                    Some(resolve_bd_input_path(&args.input)?)
                } else {
                    None
                };
            let input_format_label = if is_h5ad {
                "h5ad"
            } else {
                match detected_input.expect("checked above") {
                    DetectedInputFormat::Tenx => "10x",
                    DetectedInputFormat::BDRhapsodyDense => "bd_rhapsody",
                }
            };
            let expression_source = match detected_input {
                Some(DetectedInputFormat::BDRhapsodyDense) => {
                    if bd_input_path
                        .as_ref()
                        .and_then(|p| p.file_name().and_then(|v| v.to_str()))
                        .map(|name| {
                            name.eq_ignore_ascii_case("raw_counts.tsv")
                                || name.eq_ignore_ascii_case("raw_counts.tsv.gz")
                                || name.ends_with("_raw_counts.tsv")
                                || name.ends_with("_raw_counts.tsv.gz")
                        })
                        .unwrap_or(false)
                    {
                        ExpressionSource::RawUmiCounts
                    } else {
                        ExpressionSource::NormalizedExpression
                    }
                }
                _ => ExpressionSource::RawUmiCounts,
            };
            let is_bd_rhapsody =
                matches!(detected_input, Some(DetectedInputFormat::BDRhapsodyDense));
            if is_bd_rhapsody {
                info!(
                    path = ?bd_input_path,
                    "Detected BD Rhapsody/raw-counts dense expression format"
                );
            }

            let mut cached_result: Option<(
                Vec<String>,
                Vec<String>,
                kira_mitoqc::cache::ExpressionSoAView<'static>,
                ConfigV1,
            )> = None;
            if cache_path.exists() {
                let (features, barcodes) = if is_h5ad {
                    load_h5ad_metadata_guarded(&args.input, args.gene_symbol_key.as_deref())?
                } else if is_bd_rhapsody {
                    load_bd_rhapsody_metadata(
                        bd_input_path
                            .as_deref()
                            .expect("resolved for BD Rhapsody input"),
                    )?
                } else {
                    load_mtx_metadata(&args.input, gene_symbol_col)?
                };
                let config = load_config_autodetect(&args.assets, &features)?;
                let view = mmap_expr_bin(&cache_path)?;
                let resolved = resolve_all_genesets(
                    &GeneIndex::try_from_feature_list(&features)?,
                    &config.geneset,
                );
                let expected_genes = expected_expression_gene_count(&resolved);
                if view.mode == requested_cache_mode && view.genes == expected_genes {
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
                        cached_genes = view.genes,
                        expected_genes,
                        "cache mismatch, rebuilding expression cache"
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
                } else if is_bd_rhapsody {
                    let bd = load_bd_rhapsody(
                        bd_input_path
                            .as_deref()
                            .expect("resolved for BD Rhapsody input"),
                    )?;
                    info!(
                        features = bd.features.len(),
                        barcodes = bd.barcodes.len(),
                        nnz = bd.matrix.nnz(),
                        "Loaded BD Rhapsody dense input"
                    );
                    if args.run_mode == RunMode::Pipeline {
                        let cache_parent = if args.input.is_dir() {
                            args.input.clone()
                        } else {
                            args.input
                                .parent()
                                .map(|p| p.to_path_buf())
                                .unwrap_or_else(|| PathBuf::from("."))
                        };
                        let shared_path = cache_parent.join("kira-organelle.bin");
                        info!(
                            path = ?shared_path,
                            "Pipeline mode: reading shared cache generated by kira-organelle"
                        );
                        let shared = mmap_organelle_bin(&shared_path)?;
                        let mtx = shared.to_mtx_input();
                        (mtx.matrix, mtx.features, mtx.barcodes)
                    } else {
                        (bd.matrix, bd.features, bd.barcodes)
                    }
                } else {
                    if args.run_mode == RunMode::Pipeline {
                        let discovery = discover_dataset_files(&args.input)?;
                        let shared_name =
                            resolve_shared_cache_filename(discovery.prefix.as_deref());
                        let shared_path = discovery.input_dir.join(shared_name);
                        info!(
                            path = ?shared_path,
                            "Pipeline mode: reading shared cache generated by kira-organelle"
                        );
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

                let cluster_map = load_cluster_map_for_mode(&args, is_h5ad, &barcodes, agg_mode)?;

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
            let mut primitives = compute_primitives(&view, &resolved);
            if is_bd_rhapsody {
                let cluster_map_for_fraction =
                    load_cluster_map_for_mode(&args, is_h5ad, &barcodes, agg_mode)?;
                let mito_symbols: std::collections::BTreeSet<String> = geneset
                    .all_mtdna()
                    .into_iter()
                    .map(|s| s.to_string())
                    .collect();
                let mito_fraction = compute_mito_fraction_from_file(
                    &args.input,
                    &mito_symbols,
                    agg_mode,
                    cluster_map_for_fraction.as_ref(),
                )?;
                if mito_fraction.len() == primitives.mtdna_mean.len() {
                    primitives.mtdna_mean = mito_fraction;
                    info!(
                        expression_type = "normalized",
                        "Applied BD Rhapsody mitochondrial fraction semantics"
                    );
                } else {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "BD Rhapsody mito fraction length does not match computed sample count",
                    )
                    .into());
                }
            }

            info!("Computing proxies");
            let proxies = compute_proxies_v1(&primitives, &resolved, &config.refs)?;

            info!("Computing axes and decay");
            let scored = score_profile_v1(&proxies, &config.weights);

            // SoA-space index: maps gene symbols to their row in the
            // prepared SoA (not to be confused with `gene_index`, which
            // indexes the input feature list).
            let soa_index = SoaIndex::from_resolved(&resolved);

            let metabolic_metrics =
                compute_metabolic_metrics(&view, &soa_index, &scored.axes.ros);

            info!("Computing redox extension stage");
            let redox_metrics = if args.redox {
                Some(compute_redox_metrics(
                    &view,
                    &soa_index,
                    &scored.axes,
                    &proxies,
                )?)
            } else {
                None
            };

            info!("Classifying failure modes");
            let states = classify_v1_with_redox(
                &scored.axes,
                &scored.decay,
                &config.refs,
                redox_metrics.as_ref(),
            );

            info!("Computing explainability");
            let explain = explain_v1(
                &proxies,
                &scored.axes,
                &scored.decay,
                &states,
                &config.weights,
            );

            info!("Assembling profiles");
            let mut profiles =
                assemble_profiles_v1(&states, &scored.decay, &scored.axes, &proxies, &explain);

            if let Some(redox) = redox_metrics.as_ref() {
                for i in 0..profiles.len() {
                    if redox.low_confidence.get(i).copied().unwrap_or(false) {
                        profiles[i].interpretation.push(
                            "LOW_CONFIDENCE: redox proxy panel coverage is limited".to_string(),
                        );
                    }
                }
            }

            info!(path = ?args.out, "Writing outputs");
            write_json(&args.out, &profiles)?;
            write_axes_tsv(&args.out, &barcodes, &scored.axes)?;
            write_decay_tsv(&args.out, &scored.decay)?;
            write_proxies_tsv(&args.out, &proxies)?;
            if let Some(redox) = redox_metrics.as_ref() {
                write_redox_metrics_tsv(&args.out, &barcodes, redox)?;
            }

            if args.run_mode == RunMode::Pipeline {
                let shared_cache_name = if is_h5ad {
                    "kira-organelle.bin".to_string()
                } else {
                    discover_dataset_files(&args.input)
                        .map(|d| resolve_shared_cache_filename(d.prefix.as_deref()).to_string())
                        .unwrap_or_else(|_| "kira-organelle.bin".to_string())
                };
                let metabolic_summary = build_metabolic_summary(
                    &metabolic_metrics,
                    &barcodes,
                    matches!(agg_mode, AggregationMode::Cluster),
                );

                write_summary_json_with_redox(
                    &args.out,
                    &profiles,
                    input_format_label,
                    expression_source,
                    redox_metrics.as_ref(),
                    Some(&metabolic_summary),
                )?;
                write_mito_metrics_tsv(&args.out, &barcodes, &profiles, &metabolic_metrics)?;
                write_pipeline_step_json(&args.out, &shared_cache_name)?;
            }

            if args.version == RunVersion::V2 {
                info!("Computing v2 proxies and axes");
                let refs_v2_path = args.assets.join("refs_v2.toml");
                let weights_v2_path = args.assets.join("weights_v2.toml");
                let refs_v2 = if refs_v2_path.is_file() {
                    load_refs_v2(&refs_v2_path)?
                } else {
                    info!(
                        asset = "refs_v2.toml",
                        "asset file not found; using embedded default"
                    );
                    load_refs_v2_embedded()?
                };
                let weights_v2 = if weights_v2_path.is_file() {
                    load_weights_v2(&weights_v2_path)?
                } else {
                    info!(
                        asset = "weights_v2.toml",
                        "asset file not found; using embedded default"
                    );
                    load_weights_v2_embedded()?
                };

                let omics = OptionalOmicsInputs {
                    mt_dna_copy_number: read_vec_file(args.mtcopy.as_ref())?,
                    heteroplasmy: read_vec_file(args.heteroplasmy.as_ref())?,
                    mt_dna_deletions: read_vec_file(args.mtdeletions.as_ref())?,
                    proteomics_etc: read_vec_file(args.proteomics_etc.as_ref())?,
                    proteomics_atp: read_vec_file(args.proteomics_atp.as_ref())?,
                };

                // `proxies` is no longer needed past this point — move it into
                // compute_proxies_v2 to avoid cloning all 8 per-sample vectors.
                let proxies_v2 = compute_proxies_v2(&primitives, proxies, &refs_v2, &omics);
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
    let human_asset = assets_dir.join("geneset_v1.toml");
    let mouse_asset = assets_dir.join("geneset_mouse_v1.toml");
    let use_embedded = !human_asset.is_file() || !mouse_asset.is_file();
    let human = if use_embedded {
        info!(
            assets = %assets_dir.display(),
            "assets directory is incomplete; using embedded defaults"
        );
        ConfigV1::load_embedded_with_geneset("geneset_v1.toml")?
    } else {
        ConfigV1::load_from_assets_dir_with_geneset(assets_dir, "geneset_v1.toml")?
    };
    let mouse = if use_embedded {
        ConfigV1::load_embedded_with_geneset("geneset_mouse_v1.toml")?
    } else {
        ConfigV1::load_from_assets_dir_with_geneset(assets_dir, "geneset_mouse_v1.toml")?
    };

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

fn expected_expression_gene_count(resolved: &kira_mitoqc::input::ResolvedGeneSets) -> usize {
    resolved.mtdna_complex_i.genes.len()
        + resolved.mtdna_complex_iii.genes.len()
        + resolved.mtdna_complex_iv.genes.len()
        + resolved.mtdna_complex_v.genes.len()
        + resolved.nuclear_oxphos_complex_i.genes.len()
        + resolved.nuclear_oxphos_complex_ii.genes.len()
        + resolved.nuclear_oxphos_complex_iii.genes.len()
        + resolved.nuclear_oxphos_complex_iv.genes.len()
        + resolved.nuclear_oxphos_complex_v.genes.len()
        + resolved.ros.genes.len()
        + resolved.mitophagy.genes.len()
        + resolved.fusion.genes.len()
        + resolved.fission.genes.len()
        + resolved.biogenesis.genes.len()
        + OXPHOS_PANEL.len()
        + GLYCOLYSIS_PANEL.len()
        + FAO_PANEL.len()
        + ROS_PANEL.len()
        + BIOGENESIS_PANEL.len()
}

pub fn geneset_overlap_hits(
    features: &[String],
    geneset: &kira_mitoqc::core::types::GeneSet,
) -> usize {
    let mut feature_set: rustc_hash::FxHashSet<&str> =
        rustc_hash::FxHashSet::with_capacity_and_hasher(features.len(), Default::default());
    feature_set.extend(features.iter().map(String::as_str));

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
    use std::io::{BufRead, BufReader};
    let Some(path) = path else {
        return Ok(None);
    };
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut values = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
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

fn load_cluster_map_for_mode(
    args: &RunArgs,
    is_h5ad: bool,
    barcodes: &[String],
    agg_mode: AggregationMode,
) -> Result<Option<kira_mitoqc::data::ClusterMap>, Box<dyn std::error::Error>> {
    if agg_mode != AggregationMode::Cluster {
        return Ok(None);
    }
    let path = args.clusters.as_ref().ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "--clusters is required for cluster mode",
        )
    })?;
    if is_h5ad {
        let column = parse_h5ad_cluster_arg(path)?;
        Ok(Some(load_h5ad_clusters_guarded(&args.input, &column)?))
    } else {
        Ok(Some(load_cluster_map(path, barcodes)?))
    }
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
