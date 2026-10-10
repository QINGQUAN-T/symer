# symer —— 链式两步配对过程的符号核 v0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23282847.svg)](https://doi.org/10.5281/zenodo.23282847)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

[English version](README.en.md)

**symer** is a configurable symbolic-dynamics kernel: one engine covers both
stepwise aggregation (AA+BB polycondensation, the Flory/Carothers family) and
two-stage rumor spreading (a modified Maki–Thompson/Daley–Kendall cascade).
All channel rates carry a common *carrier factor* that cancels **exactly** in
the trajectory equations, leaving closed-form, fit-free laws — e.g. the frozen
inventory law $(x_1,x_2)=(\mu,\alpha)/(\alpha+\mu)$, inverted linearly to
recover $\alpha$ from real cascade ledgers (Higgs-Twitter: $\hat\alpha\le0.44\%$,
95% CI). Companion paper: *Three laws of the two-stage rumor cascade*
(preprint; code archived on Zenodo, DOI badge above). 48 regression tests,
all passing.

可配置的符号动力学内核：同一套引擎同时覆盖**逐步聚合**（AA+BB 缩聚，
Flory/Carothers 家族）与**两阶段谣言传播**（Maki–Thompson/Daley–Kendall
改良）两类系统。可引用代码资产（伴随论文发布）。

## 核心思想

链式两步系统的通道速率都携带一个公共**载体因子**——聚合家族是燃料池
[A]（消耗型），谣言家族是传播者池 S（催化型，`partner_consumed=False`
的"催化家族接口"）。在轨迹方程 dY/dX 中该因子**完全消去**，留下只含
无量纲参数的闭式可解方程：

| 家族 | 载体 | 引擎 | 共振线 |
|---|---|---|---|
| 聚合 | [A]（消耗） | db1/db0 = −1 + (α/2)·b1/b0 | α = 2（对数尾） |
| 谣言·合并层 | S（催化） | dS/di = −1 + μ·S/i | μ = 1（代数尾） |
| 谣言·职业层 | S（催化） | ds1/di = −1 + γ·s1/i, γ = α+μ | γ = 1（对数支） |

## 目录

```
symer/
  spec.py            Channel / ReactionSystem：状态、通道、守恒律配置化
  engine.py          轨迹推导 + 载体消去校验 + 闭式求解 + 残差归零
  configs/polymer.py AA+BB 两步缩聚（ChemRxiv 原文模型）
  configs/rumor.py   两阶段谣言级联（v1.0-lock，W1 锁定版）
tests/
  test_polymer_anchors.py  T1 Carothers / T2 Flory Đ=2 / T3 原文轨迹式（16 项）
  test_rumor_gate1.py      R1–R8 引擎闭式解、α 盲性、冻结库存律、ODE/SSA 交叉（32 项）
examples/
  demo_deviations.py       模型锁定论证演示（备忘录 §五：S×S 出传通道 + 账本指纹）
run_tests.py         一键回归（48 项全部通过，2026-10-08）
```

## 运行

```bash
pip install -r requirements.txt
python run_tests.py          # 全量回归
RUMOR_STAGE=sym python tests/test_rumor_gate1.py   # 只跑符号段
```

## 主要结果（W1 闸门 1）

1. **聚合三锚点全部复现**：Carothers Xn=(1+r)/(1+r−2pB)；α=1 时
   1−pB=√b0 且 Đ=1+p→2（Flory 最可几）；轨迹式
   b1=2(b0^{α/2}−b0)/(2−α) 及 α=2 对数支；SI 数值锚点
   X̄w(r=0.8,pB=0.45)=3.016736、冻结态 (α=3,r=0.9)：b0f=0.05676、
   b1f=0.08648、Xn=19。
2. **谣言版引擎闭式可解**：ds1/di = −1 + (α+μ)s1/i，载体 S 完全消去，
   残差恒 0（含 γ=1 对数支）。
3. **合并层 α 盲性（符号证明）**：dS/dt = k1·S·i − k2·S² 与 α 和
   S1/S2 分工无关——边际统计量完全不携带 α。
4. **冻结库存律**：终态出传者池职业构成 (x1,x2) = (μ,α)/(α+μ)，
   与种子规模无关；α = μ·x2/x1（线性反演）。
5. **饱和二分性**：μ<1 有限时刻饱和（硬冻结）；μ≥1 渐近饱和，
   S/i → 1/(μ−1)（代数尾）；μ=1 为对数临界。
6. **随机仿真交叉**：Gillespie SSA（N=6×10⁴）终态职业构成与库存律
   偏差 0.35%，满足 W2 ±2% 闸门的前置预期。

模型锁定细节与偏离说明见《W1_模型锁定与推导备忘录.md》。

---

## v0.1 新增：真实级联的 x₁ 反演 CLI（W5）

`cli/x1_invert.py` —— 在转发级联档案上执行冻结库存反演：

```bash
python3 cli/x1_invert.py --activity <activity_time.txt.gz> --rtype RT \
    --out demo --theta-hours 6 2 24
```

输出：逐级联账本（x̂₁, x̂₂）、免拟合恒等式检验（轨迹回归截距，理论 +1）、
μ̂ 与 α̂ = μ̂·x̂₂/x̂₁ 及层级 bootstrap 置信带。已在 Higgs-Twitter
（SNAP，354,930 转推事件）上验证：截距中位 1.08（θ=6h，N=298 级联），
α̂ ≤ 0.44%（95% CI）。映射约定与误差预算见论文 §VII.C 与 W5 报告。

## 复现与数据

- **论文数值**：论文中每一个数字都可由发布脚本复现（全部随机种子、参数
  网格、原始日志与结构化输出随 Zenodo 归档包发布：
  [10.5281/zenodo.23282848](https://doi.org/10.5281/zenodo.23282848)，即
  Release v0.1.0 的归档）。
- **真实级联数据**：`cli/x1_invert.py` 默认读取
  [SNAP Higgs-Twitter](https://snap.stanford.edu/data/higgs-twitter.html)
  的 `higgs-activity_time.txt.gz`（354,930 转推事件，298 条多用户级联）。
- **引用**：见 `CITATION.cff`。引用本版本：
  **Tang, Q. (2026). symer: 符号核 v0（聚合/谣言两家族的统一符号引擎）. Zenodo. [10.5281/zenodo.23282848](https://doi.org/10.5281/zenodo.23282848)**；全版本概念 DOI：[10.5281/zenodo.23282847](https://doi.org/10.5281/zenodo.23282847)。

## 许可

MIT © 2026，见 [LICENSE](LICENSE)。
