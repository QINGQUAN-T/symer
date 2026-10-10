"""偏离论证演示（配套《W1 备忘录》§五，2026-10-08 确认版）。

Demo A（偏离一）：计划书原通道（S2×S1-only 出传）的三个失败 vs 锁定版（任意 S×S hearer-stifling）
  A1 引擎结构：锁定版 ds1/di = -1 + (alpha+mu) s1/i（载体 S 消去）；
     计划书版 ds1/di = -1 + alpha s1/i + mu*s1*s2/(S*i)（含分工比，(i,s1) 不封闭）
  A2 alpha 是否泄漏进合并层：锁定版 S(t) 对 alpha=0.5/3 最大偏差 1.2e-11（alpha 盲）；
     计划书版 0.573（alpha 泄漏 57%）
  A3 饱和后全 R：锁定版 s2 -> 0（2e-5 为缓慢衰减尾，x -> 1）；
     计划书版 s2 永久滞留 0.580（僵尸传播者，终态统计无法定义）

Demo B（偏离二）：快照退化 vs 账本解析（锁定版，饱和判据 i<1e-10）
  快照 s2/(s1+s2)：alpha=0 时 0，任意 alpha>0 时 ->1（无法分辨 alpha 大小，且活体池 -> 0）
  账本 x2/(x1+x2) = alpha/(alpha+mu)：随 alpha 单调迁移，与种子 eps 无关 —— 可反演

实测输出（2026-10-08，sympy/scipy 环境，μ=0.5, ε=1e-3）：
  Demo A2: locked 1.225e-11 / plan 5.727e-01
  Demo A3: locked s1=2.0e-10 s2=2.0e-5 x=0.999980（衰减中）/ plan s2=0.580330 x=0.209835
  Demo B:  alpha=0   -> 快照 0.000000, 账本 0.000000
           alpha=0.5 -> 快照 0.999885, 账本 0.499990 (律 0.500000)
           alpha=1.0 -> 快照 0.999990, 账本 0.666660 (律 0.666667)
           alpha=3.0 -> 快照 0.999998, 账本 0.857140 (律 0.857143)
"""
import numpy as np
from scipy.integrate import solve_ivp

MU, EPS = 0.5, 1e-3
T_END = 4e5


def rhs_locked(t, y, alpha):
    """锁定版 v1.0：招募 k1=1；武装 alpha；任意 S×S hearer-出传 mu。"""
    i, s1, s2, x = y
    S = s1 + s2
    return [-S * i,
            S * i - (alpha + MU) * s1 * S,
            alpha * s1 * S - MU * s2 * S,
            MU * S * (s1 + s2)]


def rhs_plan(t, y, alpha):
    """计划书原通道：出传仅 S2×S1（S1 单亡或双亡，(i,s1,s2) 方程相同）。"""
    i, s1, s2, x = y
    S = s1 + s2
    return [-S * i,
            S * i - alpha * s1 * S - MU * s1 * s2,
            alpha * s1 * S - MU * s1 * s2,
            MU * s1 * s2]


def run(rhs, alpha, t_end=T_END):
    def ev(t, y, *args):
        return y[0] - 1e-10
    ev.terminal, ev.direction = True, -1
    s = solve_ivp(rhs, (0, t_end), [1 - EPS, EPS, 0, 0], args=(alpha,),
                  rtol=1e-11, atol=1e-16, dense_output=True, events=ev)
    return s


if __name__ == "__main__":
    print("=" * 72)
    print("Demo A1  引擎结构（符号层面）")
    print("-" * 72)
    print("锁定版:  ds1/di = [S*i - (alpha+mu)*S*s1] / (-S*i) = -1 + (alpha+mu)*s1/i")
    print("         -> S 完全消去，常系数线性 ODE，闭式解 s1 = C1*i^gamma - i/(1-gamma)")
    print("计划书版: ds1/di = -1 + alpha*s1/i + mu*s1*s2/(S*i)")
    print("         -> 末项含分工比 s2/S：二维非线性系统，(i,s1) 不封闭，无闭式轨迹")
    print()
    print("=" * 72)
    print("Demo A2  alpha 是否泄漏进合并层（mu=0.5, eps=1e-3）")
    print("-" * 72)
    ts = np.linspace(0.5, 40, 300)
    for name, rhs in (("锁定版", rhs_locked), ("计划书版", rhs_plan)):
        a = run(rhs, 0.5)
        b = run(rhs, 3.0)
        Sa = a.sol(ts)[1] + a.sol(ts)[2]
        Sb = b.sol(ts)[1] + b.sol(ts)[2]
        d = np.max(np.abs(Sa - Sb)) / np.max(Sa)
        tag = "alpha 盲" if d < 1e-6 else "alpha 泄漏进边际动力学"
        print(f"  {name}: max|S(t;alpha=0.5)-S(t;alpha=3)|/max S = {d:.3e}   <- {tag}")
    print()
    print("=" * 72)
    print("Demo A3  饱和后是否全 R（t=4e5, i 已 <1e-10）")
    print("-" * 72)
    for name, rhs in (("锁定版", rhs_locked), ("计划书版", rhs_plan)):
        s = run(rhs, 3.0)
        i, s1, s2, x = s.y[:, -1]
        note = ("s2 为缓慢衰减尾，x->1（全 R）" if name == "锁定版"
                else "s2 永久滞留（僵尸传播者，非全 R）")
        print(f"  {name}: i={i:.2e}  s1={s1:.2e}  s2={s2:.6f}  x={x:.6f}   <- {note}")
    print()
    print("=" * 72)
    print("Demo B   快照 vs 账本（锁定版，饱和判据 i<1e-10）")
    print("-" * 72)
    print(f"  {'alpha':>6} | {'s1(饱和)':>10} {'s2(饱和)':>10} {'快照s2/(s1+s2)':>14} | "
          f"{'x1':>9} {'x2':>9} {'账本x2/(x1+x2)':>14} {'定律alpha/(alpha+mu)':>20}")
    for alpha in (0.0, 0.5, 1.0, 3.0):
        s = run(rhs_locked, alpha)
        i, s1, s2, x = s.y[:, -1]
        snap = s2 / (s1 + s2) if s1 + s2 > 0 else float("nan")

        def rhs_l(t, y):
            b = rhs_locked(t, y[:4], alpha)
            return b + [MU * y[1] * (y[1] + y[2]), MU * y[2] * (y[1] + y[2])]

        sl = solve_ivp(rhs_l, (0, s.t[-1]), [1 - EPS, EPS, 0, 0, 0, 0],
                       rtol=1e-10, atol=1e-16)
        x1, x2 = sl.y[4, -1], sl.y[5, -1]
        led = x2 / (x1 + x2)
        law = alpha / (alpha + MU)
        print(f"  {alpha:>6.1f} | {s1:>10.2e} {s2:>10.2e} {snap:>14.6f} | "
              f"{x1:>9.5f} {x2:>9.5f} {led:>14.6f} {law:>20.6f}")
    print("""  判读：快照在 alpha=0 时为 0、对任意 alpha>0 都趋于 1 —— 只能区分"有无预热"，
        无法分辨 alpha 大小（且活体池本身 -> 0，无测量窗口）；
        账本 (mu,alpha)/(alpha+mu) 随 alpha 单调迁移，且与种子 eps 无关 —— 可反演。""")
