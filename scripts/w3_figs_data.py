# -*- coding: utf-8 -*-
"""W3 图级源数据存档（审稿意见 M6.4）: 论文图 3 与图 4 + 附录 C 小 g 标度。

图 3 (Law I):
  (a) 尾部塌缩: μ∈{0.5,1,2} 均值场 ODE 轨迹按式(6)预测的尾部因子缩放
      (常数由初值闭式给出, 免拟合), 平台值 -> 1;
  (b) 相平面 S/i vs i: 三分支渐近。
图 4 (Law II 指纹校准):
  Gillespie SSA N=6e4, 9 个 (μ,α) 参数点 × 3 个种子规模 eps∈{60,600,6000}
  × 6 seeds; 终态 x1/N 与 x2/(x1+x2) 对库存律 (7) 的偏差; 最大 |偏差|。
附录 C (Moderate #1) 小 g 标度:
  闭式 (9) 求根: i_f(g) ~ g^{1/μ} (μ<1); i_f ~ (μ−1)g/μ (μ>1)。

输出: results/w3_law1_data.json, results/w3_law2_data.json,
      results/w3_smallg_data.json
随机种子: seed=20261018 (SSA), 记录于 JSON。
"""
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w3_lib import gillespie

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EPS, S0, I0 = 1e-3, 1e-3, 0.999
K1 = 1.0


# ---------------- 图 3: 尾部塌缩 + 相平面 ----------------
def traj(mu, t_max=2e5):
    """均值场轨迹 (α=0: di/dt=-iS; dS/dt=iS-μS²)。返回 (t, i, S)。"""
    def rhs(t, y):
        i, S = y
        return [-i * S, i * S - mu * S * S]
    s = solve_ivp(rhs, (0, t_max), [I0, S0], method="LSODA",
                  rtol=1e-11, atol=1e-15, dense_output=True, max_step=t_max / 2e4)
    tg = np.geomspace(1e-2, min(t_max, s.t[-1]), 400)
    y = s.sol(tg)
    return tg, y[0], y[1]


def law1_data():
    out = {'meta': {'model': 'mean-field merged (alpha-blind), k1=1, i0=0.999, S0=1e-3',
                    'integrator': 'LSODA rtol=1e-11 atol=1e-15', 'seed': None}}
    panels = []
    for mu in (0.5, 1.0, 2.0):
        t, i, S = traj(mu)
        if mu < 1:
            C = (S0 + i0_term(mu)) / I0 ** mu      # S = C i^μ − i/(1−μ)
            fac = (mu * K1 * C) ** (1.0 / mu)
            scaled = i * fac * t ** (1.0 / mu)
            const = {'C': C, 'factor': fac, 'formula': 'i*(mu*k1*C)^(1/mu)*t^(1/mu)'}
        elif mu == 1:
            # S = i(C − ln i), C = S0/i0 + ln i0; 尾部 i k1 t ln t → 1
            C = S0 / I0 + np.log(I0)
            scaled = i * K1 * t * np.log(np.maximum(t, 1.0) * np.e)
            const = {'C': C, 'factor': 1.0, 'formula': 'i*k1*t*ln(e t)'}
        else:
            scaled = i * K1 * t                    # i k1 t → 1
            const = {'C': None, 'factor': 1.0, 'formula': 'i*k1*t'}
        # 平台值: 取最后 50 点均值与极差
        tail = scaled[-50:]
        # 相平面
        ratio = S / np.maximum(i, 1e-300)
        panels.append(dict(
            mu=mu, t=t.tolist(), i=i.tolist(), S_over_i=ratio.tolist(),
            scaled=scaled.tolist(), const=const,
            plateau_mean=float(np.mean(tail)), plateau_min=float(np.min(tail)),
            plateau_max=float(np.max(tail)),
            asymptote={'mu<1': 'S/i = C i^(mu-1)', 'mu=1': 'S/i = C - ln i',
                       'mu>1': 'S/i -> 1/(mu-1)'}[('mu<1' if mu < 1 else
                                                   'mu=1' if mu == 1 else 'mu>1')]))
        print(f"mu={mu}: plateau={panels[-1]['plateau_mean']:.4f} "
              f"[{panels[-1]['plateau_min']:.4f},{panels[-1]['plateau_max']:.4f}]")
    out['panels'] = panels
    json.dump(out, open(f'{ROOT}/results/w3_law1_data.json', 'w'), indent=1)
    print('-> results/w3_law1_data.json')


def i0_term(mu):
    return I0 / (1.0 - mu) if mu < 1 else 0.0


# ---------------- 图 4: 指纹校准 SSA ----------------
def law2_data():
    rng = np.random.default_rng(20261018)
    N = 60000
    mus = (0.5, 1.0, 2.0)
    alphas = (0.5, 1.0, 2.0)
    epss = (60, 600, 6000)
    seeds = 6
    pts = []
    for mu in mus:
        for al in alphas:
            g = al + mu
            for eps in epss:
                x1s, x2s, ifs = [], [], []
                for _ in range(seeds):
                    ic, x1, x2 = gillespie(N, al, mu, eps, rng)
                    x1s.append(x1 / N); x2s.append(x2 / N); ifs.append(ic / N)
                x1m, x2m = float(np.mean(x1s)), float(np.mean(x2s))
                share2 = x2m / (x1m + x2m)
                pts.append(dict(
                    mu=mu, alpha=al, eps=eps, N=N, seeds=seeds,
                    x1_mean=x1m, x2_mean=x2m, i_f_mean=float(np.mean(ifs)),
                    x1_law=mu / g, share2=share2, share2_law=al / g,
                    dev_x1_pct=(x1m - mu / g) / (mu / g) * 100.0,
                    dev_share_pct=(share2 - al / g) / (al / g) * 100.0))
                print(f"mu={mu} a={al} eps={eps}: "
                      f"dev_x1={pts[-1]['dev_x1_pct']:+.3f}% "
                      f"dev_share={pts[-1]['dev_share_pct']:+.3f}%")
    mx = max(max(abs(p['dev_x1_pct']), abs(p['dev_share_pct'])) for p in pts)
    out = {'meta': {'engine': 'w3_lib.gillespie (direct method, S-1 self-correction)',
                    'N': N, 'seeds_per_cell': seeds, 'seed': 20261018},
           'points': pts, 'max_abs_dev_pct': mx,
           'gate2_pct': 2.0}
    json.dump(out, open(f'{ROOT}/results/w3_law2_data.json', 'w'), indent=1)
    print(f'-> results/w3_law2_data.json  max|dev|={mx:.3f}% (闸门 ±2%)')


# ---------------- 附录 C: 小 g 标度 ----------------
def smallg_data():
    """闭式 (9): S(i;g) = C i^(mu-g) + g/(g-mu) − (1+g)/(1+g-mu) i,
    C 由 S(I0)=EPS 定; i_f = S 的正根 (g 小于门限时 i_f>0)。"""
    def i_f_of(mu, g):
        if abs(g - mu) < 1e-12 or abs(1.0 + g - mu) < 1e-12:
            return np.nan
        # S(i) = C i^(mu-g) + A + B i ; S(I0)=EPS
        A = g / (g - mu)
        B = -(1.0 + g) / (1.0 + g - mu)
        C = (EPS - A - B * I0) / I0 ** (mu - g)
        Sf = lambda i: C * i ** (mu - g) + A + B * i
        # 在 (0, I0) 上找符号变化并求根
        grid = np.geomspace(1e-12, I0, 4000)
        vals = Sf(grid)
        sign = np.where(vals > 0)[0]
        if len(sign) == 0:
            return 0.0
        k = sign[0]
        if k == 0:
            return float(grid[0])
        lo, hi = grid[k - 1], grid[k]
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if Sf(mid) > 0:
                hi = mid
            else:
                lo = mid
        return float(0.5 * (lo + hi))

    out = {'meta': {'formula': 'eq.(9) closed form, root i_f(g); k1=1, i0=0.999, eps=1e-3'},
           'cases': []}
    for mu, gmax, scaling in ((0.5, 0.05, 'g^(1/mu)'),
                              (1.5, 0.5, '(mu-1)g/mu'),
                              (2.0, 0.5, '(mu-1)g/mu')):
        gs = np.geomspace(1e-4, gmax, 25)
        ifs = [i_f_of(mu, g) for g in gs]
        ifs = np.array(ifs, float)
        ok = np.isfinite(ifs) & (ifs > 0)
        x, y = np.log(gs[ok]), np.log(ifs[ok])
        slope, ic = np.polyfit(x, y, 1)
        pred = 1.0 / mu if mu < 1 else 1.0
        out['cases'].append(dict(
            mu=mu, gs=gs.tolist(), i_f=ifs.tolist(), scaling=scaling,
            loglog_slope=float(slope), predicted_exponent=float(pred),
            amp=float(np.exp(ic)),
            pred_amp=(mu - 1.0) / mu if mu > 1 else None))
        print(f"mu={mu}: slope={slope:.4f} (pred {pred:.4f}), amp={np.exp(ic):.4f}")
    json.dump(out, open(f'{ROOT}/results/w3_smallg_data.json', 'w'), indent=1)
    print('-> results/w3_smallg_data.json')


if __name__ == '__main__':
    law1_data()
    law2_data()
    smallg_data()
