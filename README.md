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
recover $\alpha$ from real cascade ledgers (Higgs-Twitter pooled:
$\hat\alpha=2.21\%$, joint cascade-level cluster-bootstrap 95% CI
$[1.63,2.95]\%$; per-cascade median $0$ with bootstrap upper end $0.41\%$;
$\approx40\times$ below recruitment after calibration). Companion paper:
*Terminal-state theorems for a two-stage rumor model: closed-form
trajectories, frozen inventories, and architectural blindness* (preprint;
code archived on Zenodo, DOI badge above). 49 regression tests, all passing.

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
symer/symer/         符号核 Python 包（同目录下 import symer）
  spec.py            Channel / ReactionSystem：状态、通道、守恒律配置化
  engine.py          轨迹推导 + 载体消去校验 + 闭式求解 + 残差归零
  configs/polymer.py AA+BB 两步缩聚（ChemRxiv 原文模型）
  configs/rumor.py   两阶段谣言级联（v1.0-lock，W1 锁定版）
tests/
  test_polymer_anchors.py  T1 Carothers / T2 Flory Đ=2 / T3 原文轨迹式（16 项）
  test_rumor_gate1.py      R1–R8 引擎闭式解、α 盲性、冻结库存律、ODE/SSA 交叉（33 项）
examples/
  demo_deviations.py       模型锁定论证演示（备忘录 §五：S×S 出传通道 + 账本指纹）
scripts/             论文级可复现脚本与图级数据生成器（v0.1.1–0.1.2，见下表）
results/             图级源数据 JSON 与图 PDF（results/figs/）
logs/                原始运行日志（相对路径记录）
cli/                 x1_invert.py 真实级联反演命令行
run_tests.py         一键回归（49 项全部通过，2026-10-18）
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
5. **饱和与尾部分类**：对一切 μ>0 均为渐近饱和、不存在有限时刻饱和
   （硬冻结已被解析排除，见论文定理 1 与附录 A）；μ>1 代数尾 S/i → 1/(μ−1)；
   μ=1 对数临界；μ<1 幂律尾 i∼t^{−1/μ}。（早期文档中"μ<1 有限时刻饱和"
   的说法不准确，已在 v0.1.2 更正。）
6. **随机仿真交叉**：Gillespie SSA（N=6×10⁴，6 seeds/格点）终态职业
   构成与库存律最大偏差 0.95%（x₁ 口径；x₂ 份额口径 ≤0.32%，处于蒙特
   卡洛噪声底），满足 W2 ±2% 闸门的前置预期；种子规模扫描
   （ε=60→6000）x₁ 偏差单调收缩 0.95%→0.28%。

模型锁定细节与偏离说明见论文正文 §2 与附录 D（W1 备忘录内容已并入
论文，不再随包单独发布）。

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
池化 α̂ = 2.21%（级联级簇 bootstrap 95% CI [1.63%, 2.95%]；每级联估计
中位 0，上端 0.41%）。映射约定与误差预算见论文 §7.2–7.3 与
results/w7_measure.json。

## v0.1.1 新增：审稿修订复现脚本（W7）

响应审稿意见新增的全部数值结果，现由 `scripts/` 下的可复现脚本生成：

| 脚本 | 内容 | 论文对应 |
|---|---|---|
| `w7_sim_archives.py` | 档案模拟 SSA + 冻结库存反演 | §Simulation validation |
| `w7_nullb.py` | Null B 机制区分（耦合升级 vs 非耦合再爆发） | §7.4 |
| `w7_network_ci.py` | 网络双账本 CI 重叠 | §Network model, Fig. 7 |
| `w7_negative_controls.py` | 负对照（失配/打乱） | §Negative controls |
| `w7_gstar.py` | g* 相图与分支检查 | §S×R, Table 1 |
| `w7_measure.py` | 实测级联统计量 | §Data |
| `w7_figs.py` | 论文图 3–10 统一渲染（图级数据在 results/） | Figs. 3–10 |
| `w3_lib.py` / `w5_invert.py` | 共享引擎与反演管线 | — |

运行方式（先按 [data_README.md](data_README.md) 下载 Higgs-Twitter 数据到 `data/`）：

```bash
python3 scripts/w7_gstar.py            # 无需外部数据
python3 scripts/w7_sim_archives.py     # 需 SSA，约数分钟（--quick 可冒烟测试）
python3 scripts/w7_network_ci.py       # 合成网络，无需外部数据
python3 scripts/w7_nullb.py            # Null B 机制区分
python3 scripts/w7_measure.py          # 真实数据测量统计（需 data/）
python3 scripts/w7_negative_controls.py  # 真实数据负控制（需 data/）
python3 scripts/w3_figs_data.py        # 图 3/4 数据（闭式解 + SSA）
python3 scripts/w4_lone_spreader.py    # 图 8 数据（批量向量化 SSA）
python3 scripts/w6_network_ledger.py   # 图 5 数据（网络 CTMC）
python3 scripts/w6_failure_boundary.py # 图 6 数据（触达边界 + HMF + scale-free）
python3 scripts/w7_figs.py             # 全部论文图渲染至 results/figs/
```

各脚本的参考输出（论文数字来源）在 `results/*.json`，图级源数据与图
PDF 分别在 `results/*.json` 与 `results/figs/`，原始运行日志在 `logs/`
（日志记录相对路径；重跑时用 `python3 scripts/<脚本> 2>&1 | tee logs/<名字>.log`
即可再生成）。

## 复现与数据

- **论文数值**：论文中每一个数字都可由发布脚本复现（全部随机种子、参数
  网格、原始日志与结构化输出随 Zenodo 归档包发布：
  v0.1.0 基础归档 [10.5281/zenodo.23282848](https://doi.org/10.5281/zenodo.23282848)；
  v0.1.1 首轮审稿修订归档；v0.1.2 二轮修订归档（缺陷脚本修复、图级数据补齐、
  负控制扩充 C4），DOI 见本版本 Release 页，挂在概念 DOI
  [10.5281/zenodo.23282847](https://doi.org/10.5281/zenodo.23282847) 下）。
- **真实级联数据**：`cli/x1_invert.py` 默认读取
  [SNAP Higgs-Twitter](https://snap.stanford.edu/data/higgs-twitter.html)
  的 `higgs-activity_time.txt.gz`（354,930 转推事件，298 条多用户级联）。
- **引用**：见 `CITATION.cff`。引用最新版本（v0.1.2）请用其版本 DOI
  （见 Release 页）；引用任意版本可用概念 DOI：[10.5281/zenodo.23282847](https://doi.org/10.5281/zenodo.23282847)。

## 许可

MIT © 2026，见 [LICENSE](LICENSE)。
