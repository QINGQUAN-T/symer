"""聚合物家族回归锚点（计划书 W1 闸门 1 的前三项）。

T3  轨迹式：引擎自动推导 db1/db0 = -1 + (alpha/2) b1/b0，闭式解等于原文 (6)，
    残差恒 0；alpha=2 共振分支 -b0 ln b0。
T1  Carothers：Xn = (1+r)/(1+r-2pB)；alpha=1 时 1-pB = sqrt(b0)，r=1 退化为
    Xn = 1/(1-p)；冻结线复现（alpha=1: b0f=(1-r)^2, b1f=2r(1-r)；alpha=3, r=0.9
    数值锚点 b0f=0.05676, b1f=0.08648, Xn=19）。
T2  Flory 分散度：alpha=1, r=1 时 Đ = 1+p（SI §7 强测试 T1），p->1 极限 Đ = 2；
    数值锚点 r=0.8, pB=0.45: Xw = 3.016736（SI 双路径校验值）。
"""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from symer import derive_trajectory, solve_trajectory, verify_solution, carrier_cancellation_check
from symer.configs.polymer import build_polymer, xn_all, xw_all, freezing_line

PASS = []


def check(name, cond, detail=""):
    PASS.append((name, bool(cond)))
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")


def test_trajectory_anchor():
    """T3 原文轨迹式 (6)。"""
    sysm = build_polymer()
    b0, b1 = sysm.states["b0"], sysm.states["b1"]
    alpha = sysm.params["alpha"]

    ok_carrier, slope = carrier_cancellation_check(sysm, b0, b1)
    check("T3a 载体 [A] 完全消去", ok_carrier)

    expected = -1 + alpha * b1 / (2 * b0)
    check("T3b 轨迹方程 = -1 + (alpha/2) b1/b0", sp.simplify(slope - expected) == 0,
          f"slope = {slope}")

    sol = solve_trajectory(slope, b0, b1, ic_driver=1, ic_tracked=0)
    paper = 2 * (b0**(alpha / 2) - b0) / (2 - alpha)
    # sympy 返回 exp(alpha*log(b0)/2) 形式：残差恒 0（T3d）+ 初条件 + 数值采样
    ic_ok = sp.simplify(sol.subs(b0, 1)) == 0
    numeric_ok = all(
        abs(float(sol.subs({b0: bt, alpha: av}) - paper.subs({b0: bt, alpha: av}))) < 1e-12
        for av in (0.5, 1.7, 3.3) for bt in (0.05, 0.4, 0.9)
    )
    check("T3c 闭式解 == 原文(6)（IC + 数值采样 + 残差）", ic_ok and numeric_ok,
          f"sol = {sp.simplify(sol)}")

    res, ok = verify_solution(slope, b0, b1, sol)
    check("T3d 残差恒 0", ok)

    # alpha=2 共振分支：-b0 ln b0
    sys2 = build_polymer(alpha=sp.Integer(2))
    slope2 = derive_trajectory(sys2, b0, b1)
    sol2 = solve_trajectory(slope2, b0, b1, 1, 0)
    check("T3e alpha=2 对数尾 == -b0 ln b0",
          sp.simplify(sol2 - (-b0 * sp.log(b0))) == 0, f"sol2 = {sp.simplify(sol2)}")
    res2, ok2 = verify_solution(slope2, b0, b1, sol2)
    check("T3f alpha=2 残差恒 0", ok2)


def test_carothers_anchor():
    """T1 Carothers 方程 + 冻结线复现。"""
    b0, b1 = sp.symbols("b0 b1", nonnegative=True)
    p, r = sp.symbols("p r", positive=True)

    # Carothers 恒等式
    check("T1a Xn=(1+r)/(1+r-2pB) 恒等（r=1 退化为 1/(1-p)）",
          sp.simplify(xn_all(p, 1) - 1 / (1 - p)) == 0)

    # alpha=1: 1 - pB = sqrt(b0)
    sys1 = build_polymer(alpha=sp.Integer(1))
    slope1 = derive_trajectory(sys1, b0, b1)
    sol1 = solve_trajectory(slope1, b0, b1, 1, 0)  # b1 = 2(sqrt(b0)-b0)
    pB = 1 - b0 - sol1 / 2
    check("T1b alpha=1 时 pB = 1 - sqrt(b0)", sp.simplify(pB - (1 - sp.sqrt(b0))) == 0)

    # 冻结线：b1f = 2(1-r-b0f)；alpha=1 解析: b0f=(1-r)^2, b1f=2r(1-r)
    bl = freezing_line(b0, b1, r)
    check("T1c 冻结线配置 2b0+b1=2(1-r)", bl == sp.Eq(2 * b0 + b1, 2 * (1 - r)))
    b0f = (1 - r) ** 2
    b1f = 2 * (1 - r) - 2 * b0f  # 冻结线：2*b0f + b1f = 2(1-r)
    check("T1d alpha=1 冻结: b1f = 2r(1-r)", sp.simplify(b1f - 2 * r * (1 - r)) == 0)

    # 数值锚点（SI §8 验证记录）：alpha=3, r=0.9
    u = sp.symbols("u", positive=True)
    root = sp.nsolve(u**3 - 2 * u**2 + sp.Rational(1, 10), u, 0.24)  # b0f = u^2
    b0f_n = float(root) ** 2
    b1f_n = 2 * (0.1 - b0f_n)
    check("T1e alpha=3,r=0.9 冻结 b0f≈0.05676", abs(b0f_n - 0.05676) < 5e-5, f"{b0f_n:.6f}")
    check("T1f alpha=3,r=0.9 冻结 b1f≈0.08648", abs(b1f_n - 0.08648) < 5e-5, f"{b1f_n:.6f}")
    xn_n = float(xn_all(sp.Rational(9, 10), sp.Rational(9, 10)))
    check("T1g alpha=3,r=0.9 终态 Xn=19", abs(xn_n - 19.0) < 1e-9, f"{xn_n}")


def test_flory_dispersity_anchor():
    """T2 Flory Đ = 2（经 SI §7.7 重均公式 + 轨迹式）。"""
    b0, b1 = sp.symbols("b0 b1", nonnegative=True)
    p = sp.Symbol("p", positive=True)

    # alpha=1, r=1: b1 = 2(sqrt(b0)-b0), pB = p, b0 = (1-p)^2
    # （约定 p ∈ (0,1)：sqrt(b0) = 1-p，直接代入避免 Abs）
    b1_expr = 2 * (1 - p) - 2 * (1 - p) ** 2
    xw = xw_all(p, b1_expr, 1)
    xn = xn_all(p, 1)
    Đ = sp.simplify(xw / xn)
    check("T2a alpha=1,r=1: Đ = 1+p（SI 强测试 T1）", sp.simplify(Đ - (1 + p)) == 0,
          f"Đ = {Đ}")
    check("T2b p->1 极限 Đ = 2（Flory 最可几分布）",
          sp.limit(Đ, p, 1) == 2, f"limit = {sp.limit(Đ, p, 1)}")

    # SI 数值锚点：r=0.8, pB=0.45, alpha=1 -> b0=0.3025, b1=0.495, Xw=3.016736
    xw_n = float(xw_all(sp.Rational(45, 100), sp.Rational(495, 1000), sp.Rational(8, 10)))
    check("T2c SI 数值锚点 Xw(r=0.8,pB=0.45)=3.016736",
          abs(xw_n - 3.016736) < 1e-6, f"{xw_n:.6f}")


if __name__ == "__main__":
    test_trajectory_anchor()
    test_carothers_anchor()
    test_flory_dispersity_anchor()
    n_fail = sum(1 for _, ok in PASS if not ok)
    print(f"\n聚合物锚点：{len(PASS) - n_fail}/{len(PASS)} 通过")
    sys.exit(1 if n_fail else 0)
