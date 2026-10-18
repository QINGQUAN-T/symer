# -*- coding: utf-8 -*-
"""W7-R4 g*(mu) 精确门限 + g>mu 物理分支验证（审稿意见 #5/#3）。

均值场 S×R 扰动模型（alpha=1, k1=1, eps=1e-3）:
  di/dt = -iS;  ds1/dt = S i - (1+mu) s1 S - g s1 r;  ds2/dt = S s1 - mu s2 S - g s2 r
账本: x_k = ∫ mu s_k S dt, 沿轨道 dx_k/di = -mu s_k / i  (精确 ODE 积分)
律 (7): x1_law = mu/(1+mu), share2_law = 1/(1+mu)
门限 g*: |dev| 首次超过 2% 的 g (brentq)。
g>mu 分支: 验证闭式 S(i;g) 在 [i_f, i0] 上为正且 i_f(g)>0 存在。
输出: tmp/w7_gstar.json
"""
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/results/w7_gstar.json'
EPS, I0 = 0.001, 0.999


def integrate(mu, g):
    """高精度积分扰动模型(内嵌账本累积器), 返回 (i_f, x1, x2, share2)。
    账本口径与 SSA 一致: x_k = ∫ (mu s_k S + g s_k R) dt。"""
    def rhs(t, y):
        ii, s1, s2, r, X1, X2 = y
        S_ = s1 + s2
        return [-ii * S_,
                S_ * ii - (1.0 + mu) * S_ * s1 - g * s1 * r,
                S_ * s1 - mu * S_ * s2 - g * s2 * r,
                mu * S_ * S_ + g * S_ * r,
                mu * s1 * S_ + g * s1 * r,
                mu * s2 * S_ + g * s2 * r]
    y0 = [I0, EPS, 0.0, 0.0, 0.0, 0.0]
    s = solve_ivp(rhs, (0, 5e5), y0, method="LSODA", rtol=1e-12, atol=1e-15,
                  dense_output=False)
    i_f = float(s.y[0][-1])
    x1, x2 = float(s.y[4][-1]), float(s.y[5][-1])
    tot = x1 + x2
    return i_f, x1, x2, (x2 / tot if tot > 0 else 0.0)


def dev_share(mu, g):
    i_f, x1, x2, sh = integrate(mu, g)
    law = 1.0 / (1.0 + mu)
    return (sh - law) / law * 100.0


def dev_x1(mu, g):
    i_f, x1, x2, sh = integrate(mu, g)
    law = mu / (1.0 + mu)
    return (x1 - law) / law * 100.0


def gate(mu, f):
    """g* = |dev| 首次超过 2% 的 g（偏差可非单调, 密集扫描+局部求根）。"""
    gs = np.geomspace(1e-4, 1.5, 120)
    vals = [abs(f(mu, g)) for g in gs]
    for k in range(1, len(gs)):
        if vals[k - 1] < 2.0 <= vals[k]:
            return float(brentq(lambda gg: abs(f(mu, gg)) - 2.0,
                                gs[k - 1], gs[k], xtol=1e-5))
    return None                          # 扫描域内不过门


def main():
    res = {'gates': [], 'branch_check': []}
    for mu in (0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
        g1 = gate(mu, dev_x1)
        g2 = gate(mu, dev_share)
        res['gates'].append(dict(mu=mu, g_x1=g1, g_share=g2))
        print(f"mu={mu}: g*(x1)={g1}, g*(share2)={g2}")

    print('== g>mu branch check ==')
    for mu, g in [(0.5, 0.6), (0.5, 1.0), (0.5, 2.0), (2.0, 2.5), (2.0, 4.0)]:
        i_f, x1, x2, sh = integrate(mu, g)
        # 闭式在 [i_f, I0] 上的最小值 (数值格点)
        Csg = None
        # S(i;g) = C i^(mu-g) + g/(g-mu) - (1+g)/(1+g-mu) i, C 由 S(I0)=EPS 定
        gg = g
        C = (EPS - gg / (gg - mu) + (1 + gg) * I0 / (1 + gg - mu)) / I0 ** (mu - gg)
        xs = np.geomspace(max(i_f, 1e-12), I0, 4000)
        Ss = C * xs ** (mu - gg) + gg / (gg - mu) - (1 + gg) * xs / (1 + gg - mu)
        res['branch_check'].append(dict(mu=mu, g=g, i_f=i_f,
                                        S_min_on_bracket=float(Ss.min()),
                                        positive=bool(Ss.min() > -1e-9)))
        print(f"mu={mu} g={g}: i_f={i_f:.5f}  S_min[i_f,I0]={Ss.min():.4f} "
              f"({'物理支成立' if Ss.min() > -1e-9 else '负值!'})")
    json.dump(res, open(OUT, 'w'), indent=1)
    print('->', OUT)


if __name__ == '__main__':
    main()
