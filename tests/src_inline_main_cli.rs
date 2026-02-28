use clap::{CommandFactory, Parser};

#[path = "../src/main.rs"]
mod bin_main;

#[test]
fn run_mode_default_is_standalone() {
    let cli = bin_main::Cli::parse_from(["kira-mitoqc", "run", "--input", "in", "--out", "out"]);
    let bin_main::Commands::Run(args) = cli.command;
    assert_eq!(args.run_mode, bin_main::RunMode::Standalone);
}

#[test]
fn run_mode_pipeline_is_accepted() {
    let cli = bin_main::Cli::parse_from([
        "kira-mitoqc",
        "run",
        "--input",
        "in",
        "--out",
        "out",
        "--run-mode",
        "pipeline",
    ]);
    let bin_main::Commands::Run(args) = cli.command;
    assert_eq!(args.run_mode, bin_main::RunMode::Pipeline);
}

#[test]
fn input_format_default_is_auto() {
    let cli = bin_main::Cli::parse_from(["kira-mitoqc", "run", "--input", "in", "--out", "out"]);
    let bin_main::Commands::Run(args) = cli.command;
    assert_eq!(args.input_format, kira_mitoqc::input::InputFormat::Auto);
}

#[test]
fn input_format_bd_rhapsody_is_accepted() {
    let cli = bin_main::Cli::parse_from([
        "kira-mitoqc",
        "run",
        "--input",
        "in.tsv",
        "--out",
        "out",
        "--input-format",
        "bd-rhapsody",
    ]);
    let bin_main::Commands::Run(args) = cli.command;
    assert_eq!(
        args.input_format,
        kira_mitoqc::input::InputFormat::BdRhapsody
    );
}

#[test]
fn clap_schema_is_valid() {
    bin_main::Cli::command().debug_assert();
}

#[test]
fn geneset_overlap_prefers_mouse_symbols() {
    let human = kira_mitoqc::config::ConfigV1::load_from_assets_dir_with_geneset(
        std::path::Path::new("assets"),
        "geneset_v1.toml",
    )
    .expect("human config");
    let mouse = kira_mitoqc::config::ConfigV1::load_from_assets_dir_with_geneset(
        std::path::Path::new("assets"),
        "geneset_mouse_v1.toml",
    )
    .expect("mouse config");
    let features = vec![
        "mt-Nd1".to_string(),
        "Atp5f1a".to_string(),
        "Sod2".to_string(),
    ];
    let h = bin_main::geneset_overlap_hits(&features, &human.geneset);
    let m = bin_main::geneset_overlap_hits(&features, &mouse.geneset);
    assert!(m > h);
}
