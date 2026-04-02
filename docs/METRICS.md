# Mitochondrial-Metabolic Extension Metrics

This document defines additive mitochondrial-metabolic metrics implemented in `kira-mitoqc`.

Scope:
- single-sample compatible scRNA-seq workflows
- `--mode sample|cell|cluster`
- deterministic additive outputs in pipeline mode (`mito_metrics.tsv`, `summary.json`)

## Caveats

These metrics are **transcriptional approximations** of metabolic programs.

- They are not direct biochemical measurements.
- `OSL` is not a direct ROS readout.
- They should be interpreted with existing `kira-mitoqc` v1/v2 outputs and, when available, alongside `kira-energetics` and `kira-organelle`.

## Panel version and genes

Panel version constant:
- `MITO_METABOLIC_PANEL_V1`

Panels (filtered deterministically to genes present in annotation):

1. OXPHOS core ETC panel
- `NDUFA1..NDUFA13`, `NDUFB1..NDUFB11`
- `SDHA`, `SDHB`, `SDHC`, `SDHD`
- `UQCRC1`, `UQCRC2`
- `COX4I1`, `COX5A`, `COX6C`
- `ATP5F1A`, `ATP5F1B`, `ATP5F1C`

2. Glycolysis
- `HK1`, `HK2`, `PFKP`, `PFKM`, `ALDOA`, `GAPDH`, `PGK1`, `ENO1`, `PKM`, `LDHA`

3. Fatty-acid oxidation (FAO)
- `CPT1A`, `CPT2`, `ACADVL`, `ACADM`, `HADHA`, `HADHB`

4. ROS / oxidative stress response
- `SOD1`, `SOD2`, `GPX1`, `GPX4`
- `PRDX1..PRDX6`
- `TXN`, `TXNRD1`
- `NFE2L2` (alias-aware for `NRF2`)

5. Mitochondrial biogenesis / dynamics
- `PPARGC1A`, `TFAM`, `NRF1`, `OPA1`, `MFN1`, `MFN2`, `DNM1L` (alias-aware for `DRP1`)

## Notation

Let `E[g,c]` be per-sample expression for gene `g` and sample/cell/cluster `c` using `log1p(max(raw_expr, 0))` within this stage.

For panel `P`, values in `c` are `v_i = E[g_i,c]` for resolved panel genes.

Constants:
- `MIN_GENES = 3`
- trimmed fraction `t = 0.1`
- robust z denominator epsilon `eps = 1e-6`

## Core panel score

For each panel `P` and sample `c`:

1. If resolved gene count `< MIN_GENES` -> panel core is `NaN`.
2. Otherwise compute 10% trimmed mean:

`TM(P,c) = mean(v_sorted[k : n-k])`, where `k = floor(t*n)`.

Core outputs:
- `oxphos_core = TM(OXPHOS,c)`
- `gly_core = TM(Glycolysis,c)`
- `fao_core = TM(FAO,c)`
- `ros_core = TM(ROS,c)`
- `bio_core = TM(Biogenesis,c)`

## Robust z-score normalization

For any core series `S(c)` over all output samples:

- `m = median(S)` over finite values
- `d = median(|S - m|)` (MAD) over finite values
- If `d == 0`, z-scores for finite elements are `0`
- Else:

`Z(c) = (S(c) - m) / (1.4826*d + eps)`

## Derived metrics

Let `Zox`, `Zgly`, `Zfao`, `Zros`, `Zbio` be robust z-scores of the corresponding core series.

Let `Zmito` be robust z-score of existing mitochondrial stress axis (`axes.ros`).

### MRI (Metabolic Rigidity Index)

- `D(c) = max(Zox, Zgly, Zfao)`
- `VarZ(c) = var_pop({Zox, Zgly, Zfao})`
- `MRI(c) = D(c) - sqrt(VarZ(c))`

Interpretation: high dominant program with low cross-program dispersion suggests rigidity.

### OSL (ROS / Oxidative Stress Load)

- `OSL(c) = 0.6*Zros + 0.4*Zmito`
- fallback: if `Zmito` unavailable, `OSL(c) = Zros`

Interpretation: transcriptional oxidative stress load proxy.

### ESS (Energetic Strain Score)

- `Supply(c) = max(Zox, Zgly)`
- `ESS(c) = max(0, Zros - Supply(c))`

Interpretation: oxidative stress not matched by OXPHOS/glycolysis supply program.

### MCB (Mitochondrial Compensation Balance)

- `MCB(c) = Zbio - Zros`

Interpretation: positive values suggest adaptive compensation; negative suggests insufficient compensation.

### OGI (OXPHOS-Glycolysis Imbalance)

- `OGI(c) = Zox - Zgly`

Interpretation: positive favors OXPHOS-like dominance, negative favors glycolytic dominance.

## Thresholds and flags

Constants:
- `metabolic_rigid_high`: `MRI >= 2.0`
- `ros_high`: `OSL >= 2.0`
- `energetic_strain_high`: `ESS >= 1.5`
- `compensation_failure`: `MCB <= -1.5`

Additional internal OGI dominance cutoffs:
- OXPHOS-dominant: `OGI >= 1.5`
- Glycolysis-dominant: `OGI <= -1.5`

NaN policy:
- panel-missing dependent scores become `NaN`
- all boolean flags become `false`
- panel missingness counters are reported in `summary.json` under `mitochondrial_metabolic.missingness`

## Output contract additions

## `mito_metrics.tsv` (pipeline mode)

Added columns:
- `oxphos_core`, `gly_core`, `fao_core`, `ros_core`, `bio_core`
- `MRI`, `OSL`, `ESS`, `MCB`, `OGI`
- `metabolic_rigid_high`, `ros_high`, `energetic_strain_high`, `compensation_failure`

## `summary.json` (pipeline mode)

Added block:

- `mitochondrial_metabolic`
  - `panel_version`
  - `thresholds`
  - `global_stats` (`median`, `MAD` for MRI/OSL/ESS/MCB/OGI)
  - `cluster_stats` (if available; per-cluster med/p10/p90 + flag fractions)
  - `top_clusters_by_mri`
  - `top_clusters_by_ess`
  - `missingness`
