# symer — Symbolic Kernel for Two-Step Pairing Processes

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23282847.svg)](https://doi.org/10.5281/zenodo.23282847)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[中文说明](README.md)

**symer** is a configurable symbolic-dynamics kernel: one engine covers both
stepwise aggregation (AA+BB polycondensation, Flory/Carothers family) and
two-stage rumor spreading (modified Maki–Thompson / Daley–Kendall cascade).

Every channel rate of a chained two-step system carries a common **carrier
factor** — the fuel pool [A] (consumed) for aggregation, the spreader pool S
(catalytic) for rumors. The factor cancels **exactly** in the trajectory
equations dY/dX, leaving closed-form, parameter-free laws:

| family | carrier | trajectory equation | resonance |
|---|---|---|---|
| aggregation | [A] (consumed) | db1/db0 = −1 + (α/2)·b1/b0 | α = 2 (log tail) |
| rumor, merged ledger | S (catalytic) | dS/di = −1 + μ·S/i | μ = 1 (algebraic tail) |
| rumor, professional ledger | S (catalytic) | ds1/di = −1 + γ·s1/i, γ = α+μ | γ = 1 (log branch) |

**Frozen inventory law** (exact theorem): terminal spreader-pool composition
(x₁, x₂) = (μ, α)/(α+μ), independent of seed size — invertible linearly as
α̂ = μ·x̂₂/x̂₁. Verified SSA (N = 6×10⁴, 6 seeds/cell): max deviation 0.95%
(x₁ caliber; ≤0.32% for the x₂-share caliber, at the Monte-Carlo floor).
On real Twitter cascades (Higgs-Twitter, SNAP; 354,930 retweet events):
identity-check intercept 1.08 (theory +1), pooled α̂ = 2.21% (joint
cascade-level cluster-bootstrap 95% CI [1.63%, 2.95%]; per-cascade median 0
with bootstrap upper end 0.41%) — field arming channel ≈60× below recruitment
nominal, ≈40× after the Sec. 7.3 calibration correction.

## Install & run

```bash
pip install -r requirements.txt
python run_tests.py                                # 49 regression tests
python examples/demo_deviations.py                 # model-locking demo
RUMOR_STAGE=sym python tests/test_rumor_gate1.py   # symbolic block only
```

## Invert real cascades

```bash
python3 cli/x1_invert.py --activity higgs-activity_time.txt.gz --rtype RT \
    --out demo --theta-hours 6 2 24
```

Outputs per-cascade ledgers (x̂₁, x̂₂), the fit-free identity check
(regression intercept, theory +1), μ̂, α̂ = μ̂·x̂₂/x̂₁ with hierarchical
bootstrap bands. Data: [SNAP Higgs-Twitter](https://snap.stanford.edu/data/higgs-twitter.html).

## Layout

```
symer/symer/  spec.py (channels/conservation) · engine.py (trajectory +
            carrier cancellation + closed forms) · configs/polymer.py,
            rumor.py
tests/      16 polymer anchors (Carothers, Flory Đ=2, trajectory identities)
            + 33 rumor tests (closed forms, α-blindness, frozen inventory,
            ODE/SSA cross-validation, ledger conservation)
cli/        x1_invert.py — frozen-inventory inversion on cascade archives
examples/   demo_deviations.py
scripts/    v0.1.1: review-revision reproduction scripts (simulation
            validation, negative controls, network CI, g* phase analysis,
            Figs. 9–10) + shared engine w3_lib.py / inversion w5_invert.py
results/    reference JSON outputs (source of every number in the paper)
logs/       raw run logs (seeds 20261009–20261012)
```

Run (see `data_README.md` for the SNAP Higgs-Twitter download step):
`python3 scripts/w7_gstar.py` (no data needed) · `scripts/w7_sim_archives.py` ·
`scripts/w7_network_ci.py`.

## Citation

**Tang, Q. (2026). *symer*: symbolic kernel v0 for chained two-step pairing
processes (aggregation / rumor families). Zenodo.**
Base archive (v0.1.0): [10.5281/zenodo.23282848](https://doi.org/10.5281/zenodo.23282848);
v0.1.3 version DOI: see the Release page; concept DOI (all versions):
[10.5281/zenodo.23282847](https://doi.org/10.5281/zenodo.23282847).
Machine-readable metadata: `CITATION.cff`. License: MIT.
