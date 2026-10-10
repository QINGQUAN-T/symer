"""符号引擎：轨迹方程推导、闭式求解、残差校验、载体消去验证。

数学内核（两个家族共用的"动力学引擎"范式）
------------------------------------------
设驱动变量 X（单调耗尽，如聚合 b0、谣言无知者 i）与被跟踪变量 Y
（中间态，如 B* 库存 b1、新传播者 s1）。若所有通道速率都携带同一公共
载体因子 C（聚合 [A]；谣言 S1+S2），则

    dY/dX = (sum stoich_Y * rate) / (sum stoich_X * rate)

中 C 完全消去，轨迹方程只剩无量纲参数（聚合 alpha；谣言 gamma = alpha + mu）。
这就是闸门 1 要求的"闭式可解"的来源，也是本文与聚合论文同构的枢纽。
"""
from dataclasses import dataclass

import sympy as sp


@dataclass
class TrajectoryResult:
    driver: sp.Symbol
    tracked: sp.Symbol
    slope: sp.Expr          # d(tracked)/d(driver) 的右端 f(tracked, driver)
    solution: sp.Expr       # tracked(driver) 闭式解（含初条件）
    residual: sp.Expr       # 代回 ODE 后的残差（应为 0）
    carrier_free: bool      # 载体因子是否完全消去


# ----------------------------------------------------------------------
def derive_trajectory(system, driver, tracked):
    """从通道清单自动推导 d(tracked)/d(driver)，并校验载体消去。"""
    rhs = system.rhs()
    num = sp.simplify(rhs[tracked])
    den = sp.simplify(rhs[driver])
    slope = sp.cancel(sp.together(num / den))
    slope = sp.simplify(slope)
    return slope


def carrier_cancellation_check(system, driver, tracked):
    """载体因子必须从轨迹方程斜率中完全消失。"""
    slope = derive_trajectory(system, driver, tracked)
    return not slope.has(system.carrier), slope


# ----------------------------------------------------------------------
def solve_trajectory(slope, driver, tracked, ic_driver, ic_tracked):
    """求解 dY/dX = f(Y, X)，初条件 Y(X0)=Y0，返回闭式解表达式。

    策略：先用 sympy.dsolve 求线性一阶 ODE 的通解，再由初条件定常数；
    若 dsolve 失败，退回到幂律拟设 Y = C1 * X**m + C2 * X。
    所有符号参数保留为符号（共振/对数分支由调用方以参数特殊值另行验证）。
    """
    x, y = driver, tracked
    t = sp.Dummy("t")
    Y = sp.Function(y.name)(t)
    ode = sp.Eq(Y.diff(t), slope.subs({x: t, y: Y}))
    try:
        gensol = sp.dsolve(ode, Y)
        rhs_gen = gensol.rhs
        # dsolve 对含符号指数的线性 ODE 可能返回 Piecewise（一般支 + 共振支）；
        # 取一般支，共振/对数分支由调用方以参数特殊值另行验证
        if isinstance(rhs_gen, sp.Piecewise):
            rhs_gen = rhs_gen.args[0][0]
        # 提取任意常数并按初条件求解
        consts = sorted(rhs_gen.free_symbols - {t}
                        - set(slope.free_symbols - {y, x}),
                        key=lambda s: s.name)
        # 统一处理：把通解在 t=x0 处的值等于 y0 解常数
        unknowns = [c for c in rhs_gen.free_symbols
                    if c not in slope.free_symbols or c == y]
        unknowns = [c for c in unknowns if c != t and c != x and c != y]
        if not unknowns:
            unknowns = consts or [sp.Dummy("C1")]
        # 通解对本家族恒为常数的线性函数：按系数直接解常数，避免
        # sp.solve 在符号指数下产生 Piecewise 条件分支
        C = unknowns[0]
        collected = sp.collect(sp.expand(rhs_gen), C, evaluate=False)
        coeff = collected.get(C, sp.Integer(0))
        rest = collected.get(sp.Integer(1), sp.Integer(0))
        C_val = sp.simplify((ic_tracked - rest.subs(t, ic_driver))
                            / coeff.subs(t, ic_driver))
        sol = rhs_gen.subs(C, C_val)
        sol = sp.simplify(sp.expand(sp.powsimp(sp.simplify(sol.subs(t, x)),
                                               force=True)))
        return sol
    except Exception:
        # 幂律拟设：Y = C1 * X**m + C2 * X（覆盖本家族全部非共振情形）
        m, C1, C2 = sp.Dummy("m"), sp.Dummy("C1"), sp.Dummy("C2")
        trial = C1 * x**m + C2 * x
        eq = sp.simplify(trial.diff(x) - slope.subs(y, trial))
        poly = sp.Poly(sp.together(eq).as_numer_denom()[0], x)
        coeff_eqs = [sp.Eq(c, 0) for c in poly.all_coeffs()]
        sol_m = sp.solve(coeff_eqs, m, dict=True)
        best = None
        for cand in sol_m:
            mm = cand[m]
            if mm == 0:  # 平凡拟设不采用
                continue
            eqs = [sp.Eq((trial.diff(x) - slope.subs(y, trial)).subs(m, mm), 0)]
            sol_c = sp.solve(eqs, [C1, C2], dict=True)
            if sol_c:
                best = trial.subs(m, mm).subs(sol_c[0])
                break
        if best is None:
            raise
        # 用初条件重新定常数（拟设解可能未消耗自由度）
        consts = [c for c in best.free_symbols if c not in slope.free_symbols]
        subs = sp.solve(sp.Eq(best.subs(x, ic_driver), ic_tracked), consts, dict=True)
        if subs:
            best = best.subs(subs[0])
        return sp.simplify(best)


# ----------------------------------------------------------------------
def verify_solution(slope, driver, tracked, solution):
    """残差校验：把闭式解代回轨迹方程，应恒等于 0。"""
    x, y = driver, tracked
    lhs = sp.diff(solution, x)
    rhs = slope.subs(y, solution)
    residual = sp.simplify(lhs - rhs)
    return residual, sp.simplify(residual) == 0
