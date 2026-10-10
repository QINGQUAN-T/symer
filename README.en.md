# symer — Symbolic Kernel for Two-Step Pairing Processes

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
α̂ = μ·x̂₂/x̂₁. Verified SSA (N = 6×10⁴): max deviation 0.61%. On real
Twitter cascades (Higgs-Twitter, SNAP; 354,930 retweet events): identity-check
intercept 1.08 (theory +1), α̂ ≤ 0.44% (95% CI) — field arming channel is two
orders of magnitude weaker than recruitment.

## Install & run

```bash
pip install -r requirements.txt
python run_tests.py                                # 48 regression tests
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
symer/      spec.py (channels/conservation) · engine.py (trajectory + carrier
            cancellation + closed forms) · configs/polymer.py, rumor.py
tests/      16 polymer anchors (Carothers, Flory Đ=2, trajectory identities)
            + 32 rumor tests (closed forms, α-blindness, frozen inventory,
            ODE/SSA cross-validation)
cli/        x1_invert.py — frozen-inventory inversion on cascade archives
examples/   demo_deviations.py
```

## Citation

See `CITATION.cff`. A Zenodo DOI is minted on the first GitHub release.
License: MIT.
