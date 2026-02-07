# kira-mitoqc

Deterministic mitochondrial QC scoring for single-cell expression matrices. This repo is a staged implementation focused on stable gene sets, fixed weights, and auditable rules (no ML, no training).

## Requirements

- Rust 1.95+ (Edition 2024)
- Optional: HDF5 library for H5AD input support (`--features h5ad`)
  - macOS: `brew install hdf5`


## Installation

Install from [crates.io::kira-mitoqc](https://crates.io/crates/kira-mitoqc) (Rust 1.95+ / Windows / Linux / macOS):

```bash
cargo install kira-mitoqc
```

Or

Build from source (Rust 1.95+):

```bash
cargo build --release
```


## Layout

- `assets/` contains the versioned TOML configs (`geneset_v1.toml`, `weights_v1.toml`, `refs_v1.toml`).
- `assets/` also includes `geneset_mouse_v1.toml` for mouse datasets.
- `src/` contains the CLI and library modules.

## CLI

```bash
cargo run -- run --input ./data --out ./out --mode sample --assets ./assets
```

H5AD input (feature-gated):

```bash
cargo run --features h5ad -- run --input ./data.h5ad --out ./out --mode sample --assets ./assets
```

Run mode:

```bash
# default
cargo run -- run --input ./data --out ./out --mode sample --run-mode standalone

# pipeline shared-cache mode
cargo run -- run --input ./data --out ./out --mode sample --run-mode pipeline
```

In `pipeline` mode, `kira-mitoqc` creates a shared cache file in the input directory:
- `kira-organelle.bin` for standard 10x names.
- `<PREFIX>.kira-organelle.bin` for prefixed datasets like `GSM123_matrix.mtx`.

After writing, it reopens this cache via mmap and continues downstream computation from the mmap-backed cache data.

`CACHE_FILE.md` is the canonical format specification for this shared cache.

Geneset selection:
- `kira-mitoqc` auto-detects human vs mouse geneset from input feature symbol overlap.
- It logs the selected geneset at run start (`geneset_v1.toml` or `geneset_mouse_v1.toml`).
