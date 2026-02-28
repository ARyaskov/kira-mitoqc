# kira-mitoqc Metrics Specification

This document defines the canonical metrics produced by `kira-mitoqc`, including formulas, units, directionality, and JSON contracts.

Scope:
- v1 core scoring (`mitochondrial_profile.json`)
- optional redox extension (`--redox`)
- pipeline JSON artifacts (`summary.json`, `pipeline_step.json`)
- v2 bundle output (`mitochondrial_profile.v2.json`, when `--version v2`)

## Canonical Conventions

1. Determinism
- No stochastic steps.
- Fixed config-controlled constants and weights.

2. Sample axis semantics
- All metrics are computed per `sample`.
- `sample` means:
  - `--mode cell`: one cell/barcode
  - `--mode cluster`: one cluster aggregate
  - `--mode sample`: one dataset-wide aggregate

3. Naming
- Metric identifiers are stable snake_case strings.
- JSON uses these exact identifiers for proxy keys.

4. Value domains
- Normalized proxy scores: `[0, 1]` (clamped).
- Axis scores: expected `[0, 1]` under valid weights.
- `decay_score`: `[0, 1]` (clamped).
- `robustness_margin`: `[0, 1]`, defined as `1 - decay_score`.
- Some raw proxies are unbounded and may be negative (e.g., `ATP_coupling_loss_raw`, `NADH_imbalance_raw`).

5. Directionality
- Higher is worse for:
  - all normalized v1 proxies
  - all axis scores
  - `decay_score`
- Higher is better for:
  - `robustness_margin`
  - `redox_buffering_capacity` (optional redox stage)

## Notation

Let:
- `E[g]` = normalized expression of gene `g` (input scale).
- `mean(S)` = arithmetic mean over gene set `S`.
- `eps` = small positive constant (`refs_v1.toml`, default `1e-6`).
- `clamp01(x) = min(max(x, 0), 1)`.

Units:
- `expr` = expression units of input (typically log1p CPM/TPM contract; see refs normalization).
- `1` = dimensionless.

## Primitive Signals (Per Sample)

- `C_I = mean(mtDNA.complex_I)` (`expr`)
- `C_III = mean(mtDNA.complex_III)` (`expr`)
- `C_IV = mean(mtDNA.complex_IV)` (`expr`)
- `C_V = mean(mtDNA.complex_V)` (`expr`)
- `mtDNA_mean = mean(all mtDNA genes)` (`expr`)
- `nuclear_mean = mean(all nuclear_oxphos genes)` (`expr`)
- `ros_mean = mean(ros_detox genes)` (`expr`)
- `mitophagy_mean = mean(mitophagy genes)` (`expr`)
- `fusion_mean = mean(dynamics.fusion)` (`expr`)
- `fission_mean = mean(dynamics.fission)` (`expr`)
- `biogenesis_mean = mean(biogenesis genes)` (`expr`)
- `ATP_mt = mean(["MT-ATP6","MT-ATP8"])` (`expr`)
- `ATP_nu = E["ATP5F1A"]` (`expr`, or `0` if unresolved)
- `stoich_variance = var_pop([C_I, C_III, C_IV, C_V])` (`expr^2`)
  - implementation uses population variance: divide by `4`.

## v1 Proxy Metrics

Constants are loaded from `refs_v1.toml`.

### 1) `ETC_stoichiometry_loss`
- Raw: `stoich_variance` (`expr^2`)
- Normalized: `clamp01(stoich_variance / STOICH_REF)` (`1`)
- Meaning: ETC complex stoichiometric imbalance.

### 2) `mtDNA_expression_uncoupling`
- Raw: `abs(mtDNA_mean - nuclear_mean)` (`expr`)
- Normalized: `clamp01(raw / UNCOUPLING_REF)` (`1`)
- Meaning: mito-nuclear expression uncoupling.

### 3) `ROS_response_overdrive`
- Raw: `ros_mean` (`expr`)
- Normalized: `clamp01(raw / ROS_REF)` (`1`)
- Meaning: stress-response transcriptional overdrive.

### 4) `NADH_imbalance`
- Intermediate: `redox_proxy = C_I / (ros_mean + eps)` (`1`)
- Raw: `1 - redox_proxy` (`1`)
- Normalized: `clamp01(raw / REDOX_REF)` (`1`)
- Meaning: inferred redox pressure imbalance.

### 5) `ATP_coupling_loss`
- Intermediate: `ratio = ATP_mt / (ATP_nu + eps)` (`1`)
- Raw: `1 - ratio` (`1`)
- Normalized: `clamp01(raw / ATP_REF)` (`1`)
- Meaning: ATP coupling mismatch.

### 6) `mitophagy_excess`
- Raw: `mitophagy_mean` (`expr`)
- Normalized: `clamp01(raw / MITO_REF)` (`1`)
- Meaning: elevated mitophagy pressure.

### 7) `dynamics_imbalance`
- Intermediate: `ratio = fusion_mean / (fission_mean + eps)` (`1`)
- Raw: `abs(log2(ratio))` (`1`)
- Normalized: `clamp01(raw / DYN_REF)` (`1`)
- Meaning: fusion/fission disequilibrium.

### 8) `biogenesis_failure`
- Raw: `1 - biogenesis_mean` (`1` relative expression contrast)
- Normalized: `clamp01(raw / BIO_REF)` (`1`)
- Meaning: insufficient biogenesis compensation.

## Axis and Global Metrics

Weights are loaded from `weights_v1.toml`.

### Axis scores

`bioenergetics`:
`w1*ETC_stoichiometry_loss + w2*mtDNA_expression_uncoupling + w3*ATP_coupling_loss`

`ros`:
`w4*ROS_response_overdrive + w5*NADH_imbalance`

`dynamics`:
`w6*dynamics_imbalance + w7*mitophagy_excess`

`regulation`:
`w8*biogenesis_failure`

All axis metrics are dimensionless (`1`).

### Global metrics

`decay_score = Wb*bioenergetics + Wr*ros + Wd*dynamics + Wg*regulation`

`robustness_margin = clamp01(1 - decay_score)`

Both are dimensionless (`1`).

## Classification Metric

Thresholds come from `refs_v1.toml`.

Legacy v1 priority:
1. `ros > ROS_HIGH && bioenergetics < BIOENERGETICS_LOW` -> `ROS-dominant decay`
2. `bioenergetics > BIOENERGETICS_HIGH` -> `Bioenergetic collapse`
3. `dynamics > DYNAMICS_HIGH && regulation < REGULATION_LOW` -> `Mitophagy-locked depletion`
4. `bioenergetics > STRUCTURAL_BIO_MIN && dynamics > STRUCTURAL_DYN_MIN` -> `Structural fragmentation`
5. else -> `Compensated but fragile`

If redox extension is enabled, redox regimes have higher priority:
1. `RedoxOverload`
2. `UnbufferedOxidativeStress`
3. `CompensatedOxidativeStress`
4. fallback to legacy rules above.

## Optional Redox Metrics (`--redox`)

Output file: `mitochondrial_redox_metrics.tsv`.

Let:
- `oxidative_norm` = min-max normalized expression panel score (`[0,1]`)
- `buffering_norm` = min-max normalized buffering panel score (`[0,1]`)
- `ros_proxy = ROS_response_overdrive` (`[0,1]`)
- `nadh_proxy = NADH_imbalance` (`[0,1]`)
- `bio_fail = biogenesis_failure` (`[0,1]`)

Definitions:
- `mito_oxidative_stress_index = clamp01(0.65*oxidative_norm + 0.20*ros_proxy + 0.15*nadh_proxy)` (`1`, higher=worse)
- `redox_buffering_capacity = clamp01(0.70*buffering_norm + 0.20*(1-bio_fail) + 0.10*(1-ros_axis))` (`1`, higher=better)
- `mito_redox_mismatch = clamp(oxidative_stress - buffering_capacity, -1, 1)` (`1`)
- `mitochondrial_stress_adaptation_score = clamp01(0.45*bioenergetics + 0.25*regulation + 0.30*((mismatch+1)/2))` (`1`)

Regime classification:
- `mismatch >= 0.45 && oxidative >= 0.75` -> `RedoxOverload`
- else if `mismatch >= 0.20 && oxidative >= 0.55` -> `UnbufferedOxidativeStress`
- else if `oxidative >= 0.45 && buffering >= 0.45` -> `CompensatedOxidativeStress`
- else -> `Baseline`

Low-confidence flag:
- set when either oxidative or buffering panel coverage is `< 25%`.

## JSON Output Contract

## `mitochondrial_profile.json` (always emitted)

Top-level:
- JSON array of `MitoProfileV1`, length = number of output samples.

`MitoProfileV1` object:
- `mitochondrial_state`: string enum
  - `ROS-dominant decay`
  - `Bioenergetic collapse`
  - `Structural fragmentation`
  - `Mitophagy-locked depletion`
  - `Compensated but fragile`
  - `CompensatedOxidativeStress`
  - `UnbufferedOxidativeStress`
  - `RedoxOverload`
- `decay_score`: number (`[0,1]`)
- `robustness_margin`: number (`[0,1]`)
- `axes`: object
  - `bioenergetics`: number
  - `ros`: number
  - `dynamics`: number
  - `regulation`: number
- `proxies`: object
  - `normalized`: object map `<proxy_key -> number[]>`
  - `raw`: object map `<proxy_key -> number[]>`
  - v1 proxy keys:
    - `ETC_stoichiometry_loss`
    - `mtDNA_expression_uncoupling`
    - `ROS_response_overdrive`
    - `NADH_imbalance`
    - `ATP_coupling_loss`
    - `mitophagy_excess`
    - `dynamics_imbalance`
    - `biogenesis_failure`
  - In profile JSON, each proxy vector is length `1` (the profile sample value).
- `drivers`: array of objects
  - `key`: proxy key enum (same key space as above)
  - `axis`: enum: `Bioenergetics | Ros | Dynamics | Regulation`
  - `contribution`: number (weighted contribution to global score)
- `interpretation`: array of strings

Minimal example:

```json
[
  {
    "mitochondrial_state": "Compensated but fragile",
    "decay_score": 0.4123,
    "robustness_margin": 0.5877,
    "axes": {
      "bioenergetics": 0.44,
      "ros": 0.39,
      "dynamics": 0.33,
      "regulation": 0.21
    },
    "proxies": {
      "normalized": {
        "ETC_stoichiometry_loss": [0.31],
        "mtDNA_expression_uncoupling": [0.42],
        "ROS_response_overdrive": [0.37],
        "NADH_imbalance": [0.41],
        "ATP_coupling_loss": [0.36],
        "mitophagy_excess": [0.29],
        "dynamics_imbalance": [0.33],
        "biogenesis_failure": [0.21]
      },
      "raw": {
        "ETC_stoichiometry_loss": [0.077]
      }
    },
    "drivers": [
      {
        "key": "mtDNA_expression_uncoupling",
        "axis": "Bioenergetics",
        "contribution": 0.05145
      }
    ],
    "interpretation": [
      "Mitochondrial function maintained by compensatory mechanisms"
    ]
  }
]
```

## `summary.json` (pipeline mode)

Object contract:
- `tool`: `"kira-mitoqc"`
- `input`:
  - `mode`: `"pipeline"`
  - `n_samples`: integer
  - `input_format`: `"10x" | "bd_rhapsody" | "h5ad"`
  - `expression_type`: `"raw_umi_counts" | "normalized_expression"`
- `mitochondrial_state_distribution`: object map `<state -> fraction>` (rounded to 6 decimals)
- `decay`:
  - `decay_score_median`: number
  - `robustness_margin_median`: number
- `axes_median`:
  - `bioenergetics`, `ros`, `dynamics`, `regulation`: numbers
- optional `redox`:
  - `regime_fractions`: map for oxidative regimes
  - `mean_mito_redox_mismatch`: number
  - `high_redox_overload_fraction`: number

## `pipeline_step.json` (pipeline mode)

Object contract:
- `tool`: `"kira-mitoqc"`
- `mode`: `"pipeline"`
- `artifacts`:
  - `summary`: `"summary.json"`
  - `primary_metrics`: `"mito_metrics.tsv"`
  - `shared_cache`: `<kira-organelle.bin or prefixed variant>`
- `sample_metrics`:
  - `file`: `"mito_metrics.tsv"`
  - `id_column`: `"cell_id"`
  - `state_column`: `"mitochondrial_state"`
- `axes`: `["bioenergetics","ros","dynamics","regulation"]`

## v2 JSON (`mitochondrial_profile.v2.json`, optional)

Top-level: array of:
- `v1`: full `MitoProfileV1`
- `v2`:
  - `axes` / `decay_score` / `robustness_margin`
  - `proxies`:
    - `v1` (embedded v1 proxy block)
    - `v2_raw` map
    - `v2_normalized` map
  - `refs_version`: `"v2"`

v2-only proxy keys:
- `mtDNA_copy_number_instability`
- `mtDNA_heteroplasmy_burden`
- `mtDNA_deletion_burden`
- `proteomics_ETC_stoichiometry_loss`
- `proteomics_ATP_coupling_loss`
