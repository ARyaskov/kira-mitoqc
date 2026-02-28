//! Input specification for different aggregation modes.

use std::path::PathBuf;

use clap::ValueEnum;

/// Input aggregation mode.
#[derive(Debug, Copy, Clone, ValueEnum, PartialEq, Eq)]
pub enum InputMode {
    Sample,
    Cluster,
    Cell,
}

/// Input matrix format.
#[derive(Debug, Copy, Clone, ValueEnum, PartialEq, Eq)]
pub enum InputFormat {
    Auto,
    Tenx,
    BdRhapsody,
}

/// Input specification (path + mode).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputSpec {
    pub input_path: PathBuf,
    pub mode: InputMode,
}
