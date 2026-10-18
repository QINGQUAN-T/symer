# -*- coding: utf-8 -*-
"""W7-R1 模拟档案端到端验证（审稿意见 #7）。

已知 (alpha, mu) 的两阶段均匀混合 SSA（事件率与 w3_lib.gillespie 逐项一致），
逐事件记录 招募/升级，生成"平台式转发档案"（沉默=无事件；再激活=追加事件），
然后送入与真实数据完全相同的反演管线（w5_invert.analyze_cascade / pooled）。

档案退化（平台真实性）：
  D1 删失 c: 只保留 t <= c * span 的事件（档案在级联仍活跃时截断）
  D2 事件丢失 p: 每条事件以概率 p 独立缺失
零假设（负控制·合成侧）：alpha=0 档案 -> alpha_hat 的假阳性率。

输出: results/w7_simarchive_results.json
      （图 10 由 scripts/w7_figs.py 渲染至 results/figs/fig10-simarchive.pdf）
"""
import json, sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w5_invert import analyze_cascade, pooled  # 复用同一反演管线

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_JSON = ROOT + '/results/w7_simarchive_results.json'
# 图由 scripts/w7_figs.py 统一渲染: results/figs/fig10-simarchive.pdf（v4 稿图 10）


def ssa_events(N, alpha, mu, eps_count, rng, k1=1.0):
    """个体级两阶段 SSA（事件率与 w3_lib.gillespie 逐项一致）。
    返回 (archive_events, stifled_ids)；archive: [(user_id, t, kind)] kind=REC|ARM。
    平台语义: 招募=账号首转(首爆发), 升级=沉默后再转(再爆发)。"""
    ic, s1, s2 = N - eps_count, eps_count, 0
    i_set = set(range(eps_count, N))
    s1_set = set(range(eps_count))
    s2_set = set()
    users = list(range(N))
    rng.shuffle(users)
    events, stifled = [], []
    t = 0.0
    while True:
        S = s1 + s2
        a1 = k1 * ic * S / N
        a2 = alpha * k1 * s1 * (S - 1) / N
        a3 = mu * k1 * s1 * (S - 1) / N
        a4 = mu * k1 * s2 * (S - 1) / N
        a = a1 + a2 + a3 + a4
        if a <= 1e-12:
            break
        t += rng.exponential(1.0 / a)
        u = rng.random() * a
        if u < a1:                                   # 招募: I -> S1
            j = rng.choice(tuple(i_set))
            i_set.discard(j); s1_set.add(j)
            ic -= 1; s1 += 1
            events.append((users[j], t, 'REC'))
        elif u < a1 + a2:                            # 升级: S1 -> S2
            j = rng.choice(tuple(s1_set))
            s1_set.discard(j); s2_set.add(j)
            s1 -= 1; s2 += 1
            events.append((users[j], t, 'ARM'))
        elif u < a1 + a2 + a3:                       # S1 退休
            j = rng.choice(tuple(s1_set))
            s1_set.discard(j)
            s1 -= 1; stifled.append(users[j])
        else:                                        # S2 退休
            j = rng.choice(tuple(s2_set))
            s2_set.discard(j)
            s2 -= 1; stifled.append(users[j])
    return events, stifled


def time_warp(ev, beta=0.5):
    """平台式阵发时序: 单调重参数化 t' = T^(1-b) t^b (人类活动重尾)。
    轨道 (i,S) 不变, 仅改变分箱采样与会话间隔结构。"""
    T = max(t for _, t, _ in ev)
    if T <= 0:
        return ev
    return [(u, (t / T) ** beta * T, k) for u, t, k in ev]


def to_archive(events, rng, censor=1.0, drop=0.0):
    """真值事件 -> 平台档案: 删失 + 事件丢失; 返回 [(u,t), ...]（无 kind）。"""
    ev = [(u, t) for u, t, _ in events]
    if censor < 1.0:
        tmax = max(t for _, t in ev)
        cut = censor * tmax
        ev = [(u, t) for u, t in ev if t <= cut]
    if drop > 0:
        ev = [(u, t) for u, t in ev if rng.random() >= drop]
    return ev


def run_cell(mu, alpha, theta, M, rng, censor=1.0, drop=0.0, Nmin=120, Nmax=2000, beta=0.5):
    rows = []
    for m in range(M):
        N = int(np.clip(np.round(np.exp(rng.uniform(np.log(Nmin), np.log(Nmax)))), Nmin, Nmax))
        eps = 1
        events, stifled = ssa_events(N, alpha, mu, eps, rng)
        arch = to_archive(time_warp(events, beta), rng, censor, drop)
        c = analyze_cascade(arch, theta, 10)
        if c is not None:
            c['true_alpha'] = alpha; c['true_mu'] = mu
            c['true_x2'] = None
            rows.append(c)
    summ = pooled(rows, 800, np.random.default_rng(12345)) if rows else {}
    ic_med = float(np.median([r['intercept'] for r in rows])) if rows else float('nan')
    r2_med = float(np.median([r['r2'] for r in rows])) if rows else float('nan')
    ahat = np.array([r['alpha'] for r in rows]) if rows else np.array([np.nan])
    rel = (ahat - alpha) / alpha if alpha > 0 else ahat
    return dict(mu=mu, alpha=alpha, theta=theta, censor=censor, drop=drop,
                n_ok=len(rows), n_try=M, intercept_med=ic_med, r2_med=r2_med,
                ahat_med=float(np.nanmedian(ahat)),
                ahat_lo=float(np.nanpercentile(ahat, 25)),
                ahat_hi=float(np.nanpercentile(ahat, 75)),
                relbias_med=float(np.nanmedian(rel)) if alpha > 0 else None,
                x2_law=alpha / (alpha + mu),
                x2_pooled=(sum(r['n_bursts'][2] for r in rows) /
                           max(1, sum(r['stifled'] for r in rows))) if rows else None)


def main(quick=False):
    rng_master = np.random.default_rng(2026)
    results = {'cells': [], 'null': None, 'degradation': []}
    # --quick: 冒烟测试模式（缩减网格/种子数，仅验证脚本可运行；结果写 results/ 但
    # 与论文数字不对应，正式复现请不带该参数全量运行）
    N_REP_MAIN = 8 if quick else 200
    N_REP_CAL = 8 if quick else 60
    N_REP_NULL = 10 if quick else 400
    N_REP_DEG = 8 if quick else 200

    # --- θ 校准: 单细胞快速扫 θ ---
    print('== theta calibration ==')
    cal = []
    for th in [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0]:
        c = run_cell(1.37, 0.02, th, N_REP_CAL, rng_master)
        cal.append(c)
        print(f"  theta={th}: x2_pooled={c['x2_pooled']}, ahat_med={c['ahat_med']:.4f}, "
              f"ic={c['intercept_med']:.2f}, r2={c['r2_med']:.2f}, n_ok={c['n_ok']}")
    results['theta_calibration'] = cal
    THETAS = [0.5, 2.0, 6.0]           # 短/中/长 驻留窗（对应真实 2/6/24h 的定性角色）

    # --- 主网格: mu x alpha x theta ---
    print('== main grid ==')
    MUS = [0.5, 1.37, 2.0]
    ALPHAS = [0.0, 0.005, 0.01, 0.02, 0.04, 0.08]
    for mu in MUS:
        for al in ALPHAS:
            for th in THETAS:
                c = run_cell(mu, al, th, N_REP_MAIN, rng_master)
                results['cells'].append(c)
                print(f"  mu={mu} a={al} th={th}: ahat={c['ahat_med']:.4f} "
                      f"[{c['ahat_lo']:.4f},{c['ahat_hi']:.4f}] ic={c['intercept_med']:.2f} n={c['n_ok']}")

    # --- alpha=0 零假设: 假阳性率（α̂ > 0.44% 判为"检出"）---
    print('== alpha=0 null ==')
    null_ahats = []
    for mu in MUS:
        for th in THETAS[:1]:
            c = run_cell(mu, 0.0, th, N_REP_NULL, rng_master)
            null_ahats.append(c)
            print(f"  mu={mu} th={th}: ahat_med={c['ahat_med']:.5f}")
    results['null'] = null_ahats

    # --- 档案退化: 删失 / 事件丢失 (固定 mu=1.37, alpha=0.02, th=1.0) ---
    print('== degradation ==')
    for cen in [1.0, 0.9, 0.75]:
        c = run_cell(1.37, 0.02, 1.0, N_REP_DEG, rng_master, censor=cen)
        results['degradation'].append(c)
        print(f"  censor={cen}: ahat={c['ahat_med']:.4f} n_ok={c['n_ok']}/{c['n_try']}")
    for dr in [0.0, 0.15, 0.3]:
        if dr == 0.0:
            continue
        c = run_cell(1.37, 0.02, 1.0, N_REP_DEG, rng_master, drop=dr)
        results['degradation'].append(c)
        print(f"  drop={dr}: ahat={c['ahat_med']:.4f} n_ok={c['n_ok']}/{c['n_try']}")

    json.dump(results, open(OUT_JSON, 'w'), indent=1)
    print('->', os.path.relpath(OUT_JSON, ROOT))


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true', help='冒烟测试：缩减网格与种子数')
    a = ap.parse_args()
    main(quick=a.quick)
