# -*- coding: utf-8 -*-
"""W7-R2b Null B: 非耦合再爆发 vs 职业耦合升级（审稿意见 #8 的机制区分）。

合成 alpha=0 级联 + 退休后独立再爆发（速率 lambda，与接触过程无耦合）:
  - 终端账本反演 alpha_hat 是否"通过"（预期: alpha_hat ≈ lambda·tau_eff，即
    终端账本无法区分两种机制 —— 可识别性的诚实边界）
  - 时序判别量: 再爆发时刻 t_b 的 S(t_b)/S_max 分布
    (升级耦合: 再爆发与 S(t) 活动相关 -> 偏向大 S;
     非耦合: 与 S(t) 无关 -> S(t_b)/S_max 近似活动加权均匀)
判别统计量 D = mean[S(t_b)/S_max] - mean_t[S(t)/S_max] (耦合>0, 非耦合≈0)
输出: results/w7_nullb_results.json（图 10 相关统计见 w7_figs.py；v4 稿 §7.4）
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w7_sim_archives import ssa_events, time_warp
from w5_invert import analyze_cascade

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/results/w7_nullb_results.json'


def coupled_bursts(events, beta=0.5, theta=2.0):
    """真升级档案（w7_sim_archives 主路径）: 直接返回 warped 事件。"""
    return time_warp(events, beta)


def uncoupled_bursts(events, stifled, rng, lam, beta=0.5):
    """alpha=0 级联 + 退休用户以速率 lam 独立再爆发（指数等待, 最多3次）。"""
    ev = [(u, t, 'REC') for u, t, _ in events]
    last = {}
    for u, t, _ in ev:
        last[u] = max(t, last.get(u, -np.inf))
    T_hor = 1.5 * max(t for _, t, _ in ev)   # 视界: 超出级联时间尺度的再爆发不记录
    for u in stifled:
        t = last.get(u)
        if t is None or not np.isfinite(t):
            continue
        for _ in range(3):
            t = t + rng.exponential(1.0 / lam)
            if t > T_hor:
                break
            ev.append((u, t, 'UNC'))
    return time_warp(ev, beta)


def S_profile(events, bins=200):
    """由事件重建 S(t)/S_max (a posteriori, 停留窗 theta=2)。"""
    theta = 2.0
    users = {}
    for u, t, _ in events:
        users.setdefault(u, []).append(t)
    first = np.array([min(ts) for ts in users.values()])
    last = np.array([max(ts) for ts in users.values()]) + theta
    P = len(users)
    t0, t1 = first.min(), max(last.max(), max(t for _, t, _ in events))
    grid = np.linspace(t0, t1, bins)
    S = np.array([((first <= t) & (last > t)).sum() for t in grid], float)
    return grid, S / max(S.max(), 1e-9)


def discrim(events, bursts_idx):
    """判别量: 再爆发时刻密度与 S(t)/S_max 轮廓的 Pearson 相关 rho。
    耦合升级(速率∝S): rho>0; 非耦合(时间均匀): rho≈<=0。"""
    grid, S = S_profile(events)
    if len(bursts_idx) == 0 or S.max() <= 0:
        return 0.0, float(S.mean())
    hist, _ = np.histogram([t for _, t in bursts_idx], bins=grid)
    if hist.std() < 1e-12:
        return 0.0, float(S.mean())
    rho = float(np.corrcoef(hist, S[:-1])[0, 1])
    return rho, float(S.mean())


def main():
    rng = np.random.default_rng(20261013)
    res = {'coupled': [], 'uncoupled': []}
    for m in range(250):
        N = int(np.clip(np.round(np.exp(rng.uniform(np.log(120), np.log(2000)))), 120, 2000))
        ev, st = ssa_events(N, 0.0, 1.37, 1, rng)          # alpha=0
        # -- Null B: 注入非耦合再爆发 lam=5e-4 (单位 k1=1，与代码调用一致)
        ev_unc = uncoupled_bursts(ev, st, rng, lam=5e-4)
        arch = [(u, t) for u, t, _ in ev_unc]
        c = analyze_cascade(arch, 2.0, 16)
        if c is not None:
            bursts = [(u, t) for u, t, k in ev_unc if k == 'UNC']
            D, base = discrim(ev_unc, bursts)
            res['uncoupled'].append(dict(ahat=c['alpha'], mu=c['mu'], D=D, base=base))
        # -- 真升级对照: alpha=0.02
        ev2, st2 = ssa_events(N, 0.02, 1.37, 1, rng)
        ev2w = coupled_bursts(ev2)
        arch2 = [(u, t) for u, t, _ in ev2w]
        c2 = analyze_cascade(arch2, 2.0, 16)
        if c2 is not None:
            bursts2 = [(u, t) for u, t, k in ev2w if k == 'ARM']
            D2, base2 = discrim(ev2w, bursts2)
            res['coupled'].append(dict(ahat=c2['alpha'], mu=c2['mu'], D=D2, base=base2))
    for k in res:
        a = res[k]
        print(k, 'n=', len(a),
              '| ahat_med=%.4f' % np.median([x['ahat'] for x in a]),
              '| rho_med=%.4f' % np.median([x['D'] for x in a]))
    json.dump(res, open(OUT, 'w'), indent=1)
    print('->', os.path.relpath(OUT, ROOT))


if __name__ == '__main__':
    main()
