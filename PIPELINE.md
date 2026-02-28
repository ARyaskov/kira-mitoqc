# kira-mitoqc: Runtime Pipeline (Current Implementation)

This document describes what the utility actually does at runtime in the current codebase.

## 0) Entry point and key flags

```bash
kira-mitoqc run --input <path> --out <dir> [options]
```

Main behavior switches:

- `--mode sample|cell|cluster` (aggregation target)
- `--run-mode standalone|pipeline` (extra pipeline contract outputs)
- `--input-format auto|10x|bd-rhapsody` (format selection)
- `--cache <dir>` (location of `expr.bin`)
- `--version v1|v2` (optional additive v2 branch)
- `--redox` (optional additive redox branch)

## 1) Input path validation and format resolution

The run starts with path existence check.

Input format resolution rules:

- `.h5ad` extension => `h5ad`
- Otherwise:
  - explicit `--input-format` override is honored
  - `auto` detection:
    - directory:
      - if BD-style raw-counts file exists (`raw_counts.tsv(.gz)` or `*_raw_counts.tsv(.gz)`) => `bd_rhapsody`
      - else => `10x`
    - file:
      - BD/raw-counts-like name => `bd_rhapsody`
      - else heuristic parse of first non-empty lines for BD signature
      - fallback => `10x`

For detected `bd_rhapsody`, the effective input file is resolved with:

- direct file path, or
- directory lookup for `raw_counts.tsv(.gz)` / `*_raw_counts.tsv(.gz)`

## 2) Expression semantics tagging

The run sets an explicit expression source tag:

- `raw_umi_counts` for:
  - `10x`
  - `h5ad`
  - `bd_rhapsody` when file name indicates raw counts (`raw_counts.tsv*`)
- `normalized_expression` for other BD dense inputs

This tag is later written into `summary.json` in pipeline mode.

## 3) Cache handshake (`expr.bin`)

The expression cache path is:

- `<cache_dir>/expr.bin`, where `cache_dir = --cache` or `--out`

If cache exists:

1. Read metadata (`features`, `barcodes`) from the selected input format.
2. Load config (`assets`) with geneset autodetection (human vs mouse by overlap).
3. mmap `expr.bin`.
4. Reuse only if cache aggregation mode matches requested `--mode`.

If mode mismatches, cache is rebuilt.

## 4) Matrix loading branch

### 4.1 H5AD

- Loaded via guarded H5AD loader.
- `run-mode pipeline` is rejected for H5AD.

### 4.2 BD Rhapsody / raw-counts dense text

- Streaming parser supports:
  - comments beginning with `#`
  - both orientations:
    - gene-major (`gene<TAB>cell1...`)
    - cell-major (`barcode<TAB>gene1...`, raw-counts style)
  - gzipped `.gz` files
- Deterministic checks/warnings:
  - strict column count validation
  - duplicate genes warning
  - duplicate barcodes warning
  - non-finite values replaced by `0.0` (warning)

Pipeline mode is supported for this branch:

- shared cache `kira-organelle.bin` is written in input directory (or file parent dir),
- then reopened via mmap.

### 4.3 10x MTX

- `standalone`: load MatrixMarket + features/genes + barcodes.
- `pipeline`:
  - discover optional prefix,
  - write shared cache:
    - `kira-organelle.bin` or `<PREFIX>.kira-organelle.bin`,
  - reopen via mmap.

## 5) Expression preparation

Common across input branches:

1. Load config and resolve gene sets.
2. Optional cluster map loading for `--mode cluster`:
   - TSV barcode->cluster for non-H5AD
   - `--clusters obs:<column>` for H5AD
3. Aggregate matrix by `sample|cell|cluster`.
4. Build canonical `ExpressionSoA` in fixed gene order.
5. Write `expr.bin` with mode metadata.
6. mmap `expr.bin` for downstream compute.

## 6) Primitive signals and BD mtDNA override

Primitive signals are computed from SoA (SIMD dispatch with scalar fallback).

For `bd_rhapsody`, an explicit mtDNA fraction override is applied:

- `fraction = sum(mito_gene_expression) / sum(all_gene_expression)`
- computed for the active aggregation mode
- replaces `primitives.mtdna_mean`

This ensures no implicit library-size normalization is introduced for BD inputs.

## 7) v1 proxies, scoring, classification, explainability

Core v1 sequence:

1. `compute_proxies_v1` (with QC checks on gene-set coverage)
2. `score_profile_v1` (axes + decay)
3. optional redox stage (`--redox`)
4. `classify_v1_with_redox`
5. `explain_v1`
6. `assemble_profiles_v1`

### Redox additive stage (`--redox`)

Inputs:

- expression view, gene index, v1 axes, v1 normalized proxies
- optional panel resources from `resources/mitochondria/redox/*.tsv`

Outputs:

- `mito_oxidative_stress_index`
- `redox_buffering_capacity`
- `mito_redox_mismatch`
- `mitochondrial_stress_adaptation_score`
- `redox_regime` (`Baseline`, `CompensatedOxidativeStress`, `UnbufferedOxidativeStress`, `RedoxOverload`)
- low-confidence flag when panel coverage is limited

Classification priority with redox enabled:

1. `RedoxOverload`
2. `UnbufferedOxidativeStress`
3. `CompensatedOxidativeStress`
4. legacy v1 state rules

When low-confidence redox coverage is detected, profile interpretation gets:

- `LOW_CONFIDENCE: redox proxy panel coverage is limited`

## 8) Output writing

Always written:

- `mitochondrial_profile.json`
- `axes.tsv`
- `decay.tsv`
- `proxies.tsv`

Written only with `--redox`:

- `mitochondrial_redox_metrics.tsv`

Pipeline-mode outputs (`--run-mode pipeline`):

- `summary.json`
  - includes `input.input_format`
  - includes `input.expression_type`
  - includes additive `redox` section when redox is enabled
- `mito_metrics.tsv`
- `pipeline_step.json`

## 9) Optional v2 additive branch (`--version v2`)

Runs after v1 outputs are prepared.

Uses:

- `refs_v2.toml`
- `weights_v2.toml`
- optional omics vectors (`--mtcopy`, `--heteroplasmy`, `--mtdeletions`, `--proteomics-etc`, `--proteomics-atp`)

Writes:

- `mitochondrial_profile.v2.json`

## 10) End-to-end flow map

### 10x directory

`matrix.mtx + features/genes + barcodes` ->
detect -> load -> aggregate -> `expr.bin` (mmap) ->
primitives -> proxies -> scores -> classify -> explain ->
`mitochondrial_profile.json` + TSV outputs
[`+ summary/metrics/manifest` in pipeline mode]
[`+ v2 json` when `--version v2`]

### BD dense/raw-counts file or directory

dense file resolve -> parse -> aggregate -> `expr.bin` (mmap) ->
primitives + mtDNA fraction override ->
proxies -> scores -> classify -> explain ->
same outputs as above
[`+ summary/metrics/manifest + shared cache` in pipeline mode]
[`+ redox TSV` when `--redox`]
[`+ v2 json` when `--version v2`]

### H5AD

h5ad load -> aggregate -> `expr.bin` (mmap) ->
same v1/v2 chain

Note: H5AD pipeline mode is currently not supported.

