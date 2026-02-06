# kira-mitoqc — Functional Proxy Formulas (v1)

This document is **normative**. Any deviation requires a version bump.

The goal is to deterministically convert normalized expression data into interpretable mitochondrial decay proxies.

No learning. No adaptive thresholds. No stochasticity.

---

## General definitions

Let:

* `E[g]` — normalized gene expression (log1p CPM / TPM)
* `mean(S)` — arithmetic mean over gene set `S`
* `ε` — small constant (`1e-6`, from `refs_v1.toml`)
* `clamp(x)` — `min(max(x, 0), 1)`

All reference constants (`*_REF`) are fixed and loaded from `refs_v1.toml`.

---

## 1. ETC stoichiometry loss

Measures imbalance between mitochondrial ETC complexes rather than absolute expression.

```
C_I   = mean(mtDNA.complex_I)
C_III = mean(mtDNA.complex_III)
C_IV  = mean(mtDNA.complex_IV)
C_V   = mean(mtDNA.complex_V)

stoich_variance = variance([C_I, C_III, C_IV, C_V])

ETC_stoichiometry_loss =
  clamp(stoich_variance / STOICH_REF)
```

Interpretation:
High value indicates disrupted electron flow consistency across the respiratory chain.

---

## 2. mtDNA ↔ nuclear expression uncoupling

Detects breakdown of mitochondrial–nuclear coordination.

```
mtDNA_expr   = mean(all mtDNA genes)
nuclear_expr = mean(all nuclear_oxphos genes)

uncoupling_raw =
  abs(mtDNA_expr - nuclear_expr)

mtDNA_expression_uncoupling =
  clamp(uncoupling_raw / UNCOUPLING_REF)
```

---

## 3. ROS response overdrive

High antioxidant expression is interpreted as stress saturation, not robustness.

```
ros_level = mean(ros_detox.genes)

ROS_response_overdrive =
  clamp(ros_level / ROS_REF)
```

---

## 4. NADH / redox imbalance (proxy)

Indirect redox pressure estimation using ETC–detox balance.

```
redox_proxy =
  C_I / (ros_level + ε)

NADH_imbalance =
  clamp(1 - redox_proxy / REDOX_REF)
```

---

## 5. ATP coupling loss

Detects decoupling between ETC activity and ATP synthesis.

```
ATP_mt = mean(["MT-ATP6", "MT-ATP8"])
ATP_nu = E["ATP5F1A"]

ATP_coupling_raw =
  ATP_mt / (ATP_nu + ε)

ATP_coupling_loss =
  clamp(1 - ATP_coupling_raw / ATP_REF)
```

---

## 6. Mitophagy excess

Represents sustained mitophagy pressure rather than healthy turnover.

```
mitophagy_signal = mean(mitophagy.genes)

mitophagy_excess =
  clamp(mitophagy_signal / MITO_REF)
```

---

## 7. Mitochondrial dynamics imbalance

Fusion–fission disequilibrium.

```
fusion  = mean(dynamics.fusion)
fission = mean(dynamics.fission)

dynamics_ratio =
  fusion / (fission + ε)

dynamics_imbalance =
  clamp(abs(log2(dynamics_ratio)) / DYN_REF)
```

---

## 8. Biogenesis failure

Failure of compensatory mitochondrial renewal.

```
biogenesis_signal = mean(biogenesis.genes)

biogenesis_failure =
  clamp(1 - biogenesis_signal / BIO_REF)
```

---

## 9. Axis aggregation

### Bioenergetics axis

```
bioenergetics =
  w1 * ETC_stoichiometry_loss +
  w2 * mtDNA_expression_uncoupling +
  w3 * ATP_coupling_loss
```

### ROS axis

```
ros =
  w4 * ROS_response_overdrive +
  w5 * NADH_imbalance
```

### Dynamics axis

```
dynamics =
  w6 * dynamics_imbalance +
  w7 * mitophagy_excess
```

### Regulation axis

```
regulation =
  w8 * biogenesis_failure
```

All axis weights are fixed and sum to 1 within each axis.

---

## 10. Global decay score

```
DecayScore =
  Wb * bioenergetics +
  Wr * ros +
  Wd * dynamics +
  Wg * regulation
```

Global weights are fixed and sum to 1.

---

## 11. Failure mode classification (v1 rules)

Deterministic rule-based classification.

```
if ros > 0.75 and bioenergetics < 0.5:
    state = "ROS-dominant decay"

elif bioenergetics > 0.8:
    state = "Bioenergetic collapse"

elif dynamics > 0.7 and regulation < 0.4:
    state = "Mitophagy-locked depletion"

elif bioenergetics > 0.6 and dynamics > 0.6:
    state = "Structural fragmentation"

else:
    state = "Compensated but fragile"
```

Thresholds are loaded from `refs_v1.toml` and must not be hard-coded.

---

## Design guarantees

* Fully deterministic
* Explainable by construction
* Stable across datasets
* SIMD- and mmap-friendly
* Extendable in v2 (proteomics, heteroplasmy) without breaking v1
