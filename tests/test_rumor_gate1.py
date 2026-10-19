"""谣言配置 Gate-1 验证（W1 闸门 1 的谣言侧）。

R1  引擎推导：ds1/di = -1 + (alpha+mu) s1/i，载体 S=s1+s2 完全消去；
R2  合并层 alpha 盲性（符号证明）：dS/dt = k1 S i - k2 S^2，不含 alpha；
R3  合并层轨迹：dS/di = -1 + mu S/i（与聚合物引擎同族，alpha/2 <-> mu）；
R4  职业层闭式解 + 残差恒 0（含 gamma=1 对数分支）；
R5  冻结库存律（符号）：x1 = mu/(alpha+mu) 与种子无关；x2 = alpha/(alpha+mu)；
R6  数值交叉：ODE 积分 vs 闭式解 vs 库存律（mu<1 幂律尾渐近饱和）；
R7  mu>=1 代数尾：S/i -> 1/(mu-1)，库存律仍收敛；
R8  Gillespie SSA 抽查：终态出传者职业构成 vs (mu, alpha)/(alpha+mu)。
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from symer import derive_trajectory, solve_trajectory, verify_solution, carrier_cancellation_check
from symer.configs.rumor import build_rumor, career_slope, pooled_slope

PASS = []


def check(name, cond, detail=""):
    PASS.append((name, bool(cond)))
    print(f"[{'PASS' if cond else 'FAIL'}] {name} {detail}")


# ----------------------------------------------------------------------
def test_engine_and_blindness():
    """R1-R3：引擎推导、载体消去、合并层 alpha 盲性（符号级）。"""
    sysm = build_rumor()
    i, s1, s2, x = (sysm.states[n] for n in ("i", "s1", "s2", "x"))
    alpha, mu, k1 = (sysm.params[n] for n in ("alpha", "mu", "k1"))
    k2 = mu * k1

    ok_carrier, slope = carrier_cancellation_check(sysm, i, s1)
    check("R1a 载体 S=s1+s2 完全消去", ok_carrier)
    gamma = alpha + mu
    expected = -1 + gamma * s1 / i
    check("R1b 职业层轨迹方程 = -1 + (alpha+mu) s1/i",
          sp.simplify(slope - expected) == 0, f"slope = {sp.simplify(slope)}")

    # 合并层 alpha 盲性：S = s1+s2 的方程与 alpha 和分工无关
    rhs = sysm.rhs()
    dSdt = sp.simplify(rhs[s1] + rhs[s2])
    pooled_expected = k1 * (s1 + s2) * i - k2 * (s1 + s2) ** 2
    check("R2a 合并层 dS/dt = k1 S i - k2 S^2", sp.simplify(dSdt - pooled_expected) == 0)
    check("R2b alpha 盲（d(dS/dt)/d alpha = 0）", not dSdt.has(alpha))
    # 分工无关性：对 s1 的任意分工扰动，dS/dt 不变
    dSdt_s2 = sp.simplify(dSdt.subs(s2, sp.Symbol("S") - s1))
    check("R2c 合并层与 S1/S2 分工无关", not dSdt_s2.has(s1))

    # 合并层轨迹方程与聚合物引擎同构：dS/di = -1 + mu S/i
    dSdi = sp.simplify(dSdt / rhs[i])
    check("R3 合并层轨迹 = -1 + mu S/i（聚合物家族同构，alpha/2 <-> mu）",
          sp.simplify(dSdi - (-1 + mu * (s1 + s2) / i)) == 0, f"dS/di = {dSdi}")


def test_career_closed_form():
    """R4：职业层闭式解 + 残差（gamma≠1 幂律支 与 gamma=1 对数支）。"""
    i, s1 = sp.symbols("i s1", nonnegative=True)
    alpha = sp.Symbol("alpha", positive=True)
    mu = sp.Symbol("mu", positive=True)
    gamma = alpha + mu
    eps = sp.Symbol("epsilon", positive=True)
    i0 = 1 - eps

    slope = -1 + gamma * s1 / i
    sol = solve_trajectory(slope, i, s1, ic_driver=i0, ic_tracked=eps)
    expected = (eps + i0 / (1 - gamma)) / i0**gamma * i**gamma - i / (1 - gamma)
    sym_ok = sp.simplify(sp.simplify(sol - expected)) == 0
    num_ok = all(
        abs(complex(sol.subs({i: iv, alpha: av, mu: mv, eps: ev})).real
            - complex(expected.subs({i: iv, alpha: av, mu: mv, eps: ev})).real) < 1e-10
        for av, mv in ((0.8, 0.3), (2.0, 0.5), (1.5, 0.25))  # 避开 gamma=1 共振点（R4c/D 单测）
        for ev in (1e-3,) for iv in (0.9, 0.3, 0.05)
    )
    check("R4a 闭式解 = C1 i^gamma - i/(1-gamma)（符号或数值采样）",
          sym_ok or num_ok, f"sol = {sp.simplify(sol)}")
    _, ok = verify_solution(slope, i, s1, sol)
    check("R4b 残差恒 0", ok)

    # gamma=1 对数分支
    slope1 = -1 + s1 / i
    sol1 = solve_trajectory(slope1, i, s1, ic_driver=i0, ic_tracked=eps)
    check("R4c gamma=1 对数支 s1 = i(C1 - ln i)",
          sp.simplify(sol1 - i * (eps / i0 + sp.log(i0) - sp.log(i))) == 0,
          f"sol1 = {sp.simplify(sol1)}")
    _, ok1 = verify_solution(slope1, i, s1, sol1)
    check("R4d gamma=1 残差恒 0", ok1)


def test_inventory_law_symbolic():
    """R5：冻结库存律 x1 = mu/gamma（种子无关），x2 = alpha/gamma。"""
    alpha = sp.Symbol("alpha", positive=True)
    mu = sp.Symbol("mu", positive=True)
    gamma = alpha + mu
    eps = sp.Symbol("epsilon", positive=True)
    i0 = 1 - eps
    i = sp.Symbol("i", positive=True)

    s1 = (eps + i0 / (1 - gamma)) / i0**gamma * i**gamma - i / (1 - gamma)
    # x1 = mu \int_0^{i0} s1(i)/i di：初等原函数 F(ig)，取端点差（gamma>0 收敛于 0 处）
    ig = sp.Symbol("ig", positive=True)
    F = sp.simplify(((eps + i0 / (1 - gamma)) / i0**gamma) * ig**gamma / gamma
                    - ig / (1 - gamma))
    x1 = sp.simplify(mu * (F.subs(ig, i0) - sp.limit(F, ig, 0, "+")))
    check("R5a x1 = mu/(alpha+mu)", sp.simplify(x1 - mu / gamma) == 0, f"x1 = {x1}")
    check("R5b 库存律与种子 epsilon 无关", not sp.simplify(x1).has(eps))
    # gamma=1 极限连续性：x1 -> mu
    x1_gamma1 = sp.limit(x1, alpha, 1 - mu)
    check("R5c gamma=1 极限 x1 -> mu（对数支连续）",
          sp.simplify(x1_gamma1 - mu) == 0, f"limit = {x1_gamma1}")


# ----------------------------------------------------------------------
def _rhs_num(t, y, alpha, mu):
    i, s1, s2, xx = y
    S = s1 + s2
    return [-S * i,            # di/dt   (k1 = 1)
            S * i - (alpha + mu) * S * s1,   # ds1/dt
            alpha * S * s1 - mu * S * s2,    # ds2/dt
            mu * S * S]        # dx/dt


def _closed_s1(i, alpha, mu, eps):
    i0 = 1 - eps
    gamma = alpha + mu
    if abs(gamma - 1) < 1e-12:
        return i * (eps / i0 + np.log(i0) - np.log(i))
    return (eps + i0 / (1 - gamma)) / i0**gamma * i**gamma - i / (1 - gamma)


def run_mean_field(alpha, mu, eps, t_max=2e4):
    """数值积分平均场 ODE，返回终端量与轨迹采样。"""
    i0, s10 = 1 - eps, eps

    def event_saturation(t, y, *a):
        return y[0] - 1e-8
    event_saturation.terminal = True
    event_saturation.direction = -1

    sol = solve_ivp(_rhs_num, (0, t_max), [i0, s10, 0.0, 0.0],
                    args=(alpha, mu), rtol=1e-11, atol=1e-14,
                    events=event_saturation, dense_output=True)
    yf = sol.y[:, -1]
    hit = len(sol.t_events[0]) > 0
    return sol, hit, yf


def test_numeric_crosscheck():
    """R6：mu<1 渐近饱和（无有限时刻冻结，定理 1）+ 库存律数值验证 + 闭式解沿轨迹吻合。"""
    cases = [(0.8, 0.3, 1e-3, 1.1), (0.6, 0.4, 1e-3, 1.0), (2.0, 0.5, 1e-2, 2.5)]
    for alpha, mu, eps, gamma in cases:
        sol, hit, yf = run_mean_field(alpha, mu, eps)
        i_f, s1_f, s2_f, x_f = yf
        x1_expect = mu / gamma
        x2_expect = alpha / gamma
        tag = f"(alpha={alpha}, mu={mu}, eps={eps})"
        if mu < 1:
            check(f"R6a{tag} mu<1 渐近饱和（数值 i 低于 1e-8；无有限时刻冻结）", hit and i_f < 2e-8,
                  f"t_end={sol.t[-1]:.3f}")
        # 终端账本：出传者池分解 x1, x2——用闭式积分流验证：
        # 数值上以 -mu*s1/i 与 -mu*s2/i 对 i 积分（伴随方程）得到 x1, x2
        def rhs_ledger(t, y, *a):
            base = _rhs_num(t, y[:4], alpha, mu)
            i, s1, s2, xx = y[:4]
            S = s1 + s2
            dx1di = 0.0 if i <= 0 or S <= 0 else mu * s1 / i
            dx2di = 0.0 if i <= 0 or S <= 0 else mu * s2 / i
            # 随 t 正向积流出：dx1/dt = (dx1/di)(di/dt)
            return base + [-dx1di * base[0], -dx2di * base[0]]
        solL = solve_ivp(rhs_ledger, (0, sol.t[-1]), [1 - eps, eps, 0.0, 0.0, 0.0, 0.0],
                         rtol=1e-10, atol=1e-14)
        x1_n, x2_n = solL.y[4, -1], solL.y[5, -1]
        check(f"R6b{tag} 冻结库存 x1=mu/gamma",
              abs(x1_n - x1_expect) < 1e-5, f"x1={x1_n:.8f} vs {x1_expect:.8f}")
        # x2 账本在截断时刻 ie 缺失幂律尾：mu∫_0^{ie}(S/i)di = C_S ie^mu - mu*ie/(1-mu)
        # （S(i) = C_S i^mu - i/(1-mu)），修正后应精确回到定律值
        ie = solL.y[0, -1]
        C_S = (eps + (1 - eps) / (1 - mu)) / (1 - eps) ** mu
        tail = C_S * ie ** mu - mu * ie / (1 - mu) if mu < 1 else mu / (mu - 1) * ie
        check(f"R6c{tag} 冻结库存 x2=alpha/gamma（含解析尾修正）",
              abs(x2_n + tail - x2_expect) < 1e-6,
              f"x2={x2_n:.8f}+tail{tail:.2e} vs {x2_expect:.8f}")
        # 闭式解沿轨迹吻合
        ts = np.linspace(sol.t[0], sol.t[-1] * 0.98, 200)
        yy = sol.sol(ts)
        err = 0.0
        for t, y in zip(ts, yy.T):
            if y[0] > 1e-8:
                err = max(err, abs(y[1] - _closed_s1(y[0], alpha, mu, eps)))
        check(f"R6d{tag} 闭式解沿轨迹最大偏差 < 1e-7", err < 1e-7, f"max|ds1|={err:.2e}")


def test_algebraic_tail():
    """R7：mu>=1 无限时饱和（代数尾），账本积分仍收敛于库存律。"""
    alpha, mu, eps = 0.5, 2.0, 1e-3
    gamma = alpha + mu
    sol, hit, yf = run_mean_field(alpha, mu, eps, t_max=5e4)
    i_f, s1_f, s2_f, x_f = yf
    S_f = s1_f + s2_f
    ratio_S = S_f / i_f if i_f > 0 else np.inf
    ratio_s1 = s1_f / i_f if i_f > 0 else np.inf
    check("R7a mu>1 不发生有限时刻饱和", not hit, f"i(t_max)={i_f:.2e}")
    check("R7b 代数尾 S/i -> 1/(mu-1)=1", abs(ratio_S - 1.0) < 2e-3, f"S/i={ratio_S:.6f}")
    check("R7c 代数尾 s1/i -> 1/(gamma-1)", abs(ratio_s1 - 1 / (gamma - 1)) < 5e-3,
          f"s1/i={ratio_s1:.6f} vs {1/(gamma-1):.6f}")
    # 账本积分到 t_max：尾部缺失通量 O(i(t_max))，x1(t_max) -> mu/gamma
    def rhs_ledger(t, y):
        base = _rhs_num(t, y[:4], alpha, mu)
        i, s1, s2 = y[0], y[1], y[2]
        S = s1 + s2
        dx1di = mu * s1 / i if i > 1e-300 and S > 0 else 0.0
        dx2di = mu * s2 / i if i > 1e-300 and S > 0 else 0.0
        return base + [-dx1di * base[0], -dx2di * base[0]]
    solL = solve_ivp(rhs_ledger, (0, sol.t[-1]),
                     [1 - eps, eps, 0.0, 0.0, 0.0, 0.0],
                     rtol=1e-10, atol=1e-14)
    x1_n, x2_n = solL.y[4, -1], solL.y[5, -1]
    check("R7d 账本 x1 -> mu/gamma（尾部通量可忽略）",
          abs(x1_n - mu / gamma) < 2e-3, f"x1={x1_n:.6f} vs {mu/gamma:.6f}")
    check("R7e 账本 x2 -> alpha/gamma",
          abs(x2_n - alpha / gamma) < 2e-3, f"x2={x2_n:.6f} vs {alpha/gamma:.6f}")


# ----------------------------------------------------------------------
def gillespie(N, alpha, mu, eps_count, rng):
    """直接法 SSA：返回 (i_f, x1, x2)（计数）。"""
    I, S1, S2, R = N - eps_count, eps_count, 0, 0
    x1 = x2 = 0
    k1 = 1.0
    k2 = mu * k1
    while S1 + S2 > 0:
        a1 = k1 * (S1 + S2) * I / N                      # REC
        a2 = alpha * k1 * (S1 * (S1 - 1) + S1 * S2) / N  # ARM（同态对自校正）
        a3 = k2 * S1 * (S1 + S2 - 1) / N                 # STIF hearer=S1
        a4 = k2 * S2 * (S1 + S2 - 1) / N                 # STIF hearer=S2
        a = a1 + a2 + a3 + a4
        if a <= 0:
            break
        tau = rng.exponential(1 / a)
        u = rng.random() * a
        if u < a1:
            I -= 1; S1 += 1
        elif u < a1 + a2:
            S1 -= 1; S2 += 1
        elif u < a1 + a2 + a3:
            S1 -= 1; R += 1; x1 += 1
        else:
            S2 -= 1; R += 1; x2 += 1
    return I, x1, x2


def test_ssa_sanity():
    """R8：SSA 终态出传者职业构成 vs (mu, alpha)/(alpha+mu)。"""
    alpha, mu = 1.0, 0.5
    gamma = alpha + mu
    N, seeds = 60000, 4
    eps_count = 60
    x1_tot = x2_tot = i_tot = 0
    rng = np.random.default_rng(20260924)
    for _ in range(seeds):
        I, x1, x2 = gillespie(N, alpha, mu, eps_count, rng)
        x1_tot += x1; x2_tot += x2; i_tot += I
    frac2 = x2_tot / (x1_tot + x2_tot)
    law = alpha / gamma
    i_f_frac = i_tot / (N * seeds)
    check("R8a SSA 职业构成 x2/(x1+x2) ≈ alpha/(alpha+mu)（±3%）",
          abs(frac2 - law) < 0.03, f"x2份={frac2:.4f} vs 律={law:.4f}")
    check("R8b SSA 终态无知者占比小（mu<1 近饱和）", i_f_frac < 0.10,
          f"i_f/N = {i_f_frac:.4f}")
    check("R8b2 SSA 终态账本守恒 x1+x2=N-1（定理 8，有限 N 恰余一个传播者）",
          (x1_tot + x2_tot) == (N - 1) * seeds,
          f"x1+x2={x1_tot + x2_tot} vs (N-1)*seeds={(N - 1) * seeds}")
    print(f"      SSA 详情: N={N}, seeds={seeds}, i_f/N={i_f_frac:.4f}, "
          f"x1={x1_tot}, x2={x2_tot}")


if __name__ == "__main__":
    import os
    stage = os.environ.get("RUMOR_STAGE", "all")
    if stage in ("sym", "all"):
        test_engine_and_blindness()
        test_career_closed_form()
        test_inventory_law_symbolic()
    if stage in ("num", "all"):
        test_numeric_crosscheck()
        test_algebraic_tail()
    if stage in ("ssa", "all"):
        test_ssa_sanity()
    n_fail = sum(1 for _, ok in PASS if not ok)
    print(f"\n谣言 Gate-1[{stage}]：{len(PASS) - n_fail}/{len(PASS)} 通过")
    sys.exit(1 if n_fail else 0)
