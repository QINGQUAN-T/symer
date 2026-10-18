# -*- coding: utf-8 -*-
"""W7-R3 网络账本偏差的置信区间与参数依赖（审稿意见 #6）。

A) 指数置信区间: Poisson <k>∈{4,6,8,12,20,40}, mu=0.5, alpha=1, N=2e4。
   每点 B=10 个独立图实现 × 3 SSA 种子 -> 每实现账本偏差 ->
   log-log 幂律拟合 + bootstrap 指数 CI。
B) 参数依赖 (<k>=8): mu∈{0.5,2}×alpha∈{0.5,1,3}; N∈{5e3,2e4}。
输出: results/w7_network_ci.json
      （图 7 由 scripts/w7_figs.py 渲染至 results/figs/fig7-networkci.pdf；v4 稿图 7）
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w3_lib import NetSim, build_graph

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/results/w7_network_ci.json'


def ledger_dev(adj, n, alpha, mu, eps, rng, runs=3):
    devs, share2s = [], []
    for _ in range(runs):
        sim = NetSim(adj, alpha, mu, eps, rng)
        r = sim.run()
        s2 = r['x2'] / max(1, r['x1'] + r['x2'])
        law = alpha / (alpha + mu)
        devs.append(abs(s2 - law) / law * 100)
        share2s.append(s2)
    return float(np.mean(devs)), float(np.mean(share2s))


def main():
    rng = np.random.default_rng(20261013)
    res = {'exp': [], 'scan': []}

    print('== A) exponent CI (Poisson, mu=0.5, alpha=1, N=2e4) ==')
    ks = (4, 6, 8, 12, 20, 40)
    devs_by_k = {}
    for zt in ks:
        per_real = []
        for b in range(10):
            adj, z, n = build_graph('poisson', 20000, mean_k=zt, rng=rng)
            eps = max(1, round(0.001 * n))
            d, s2 = ledger_dev(adj, n, 1.0, 0.5, eps, rng, runs=3)
            per_real.append(d)
        devs_by_k[zt] = per_real
        row = dict(k=zt, dev_mean=float(np.mean(per_real)),
                   dev_sd=float(np.std(per_real)),
                   dev_med=float(np.median(per_real)), n_real=len(per_real))
        res['exp'].append(row)
        print(f"  <k>={zt}: dev={row['dev_mean']:.2f}% ± {row['dev_sd']:.2f}")

    # 幂律拟合 + bootstrap CI（对实现重采样）
    ks_a = np.array([r['k'] for r in res['exp']], float)
    dm = np.array([r['dev_mean'] for r in res['exp']])
    A = np.column_stack([np.ones_like(np.log(ks_a)), np.log(ks_a)])
    coef, *_ = np.linalg.lstsq(A, np.log(dm), rcond=None)
    exp_hat = -coef[1]
    boot = []
    for _ in range(4000):
        idx = rng.integers(0, len(ks_a), len(ks_a))
        kb, db = ks_a[idx], dm[idx]
        if np.any(db <= 0):
            continue
        cb, *_ = np.linalg.lstsq(np.column_stack([np.ones_like(np.log(kb)), np.log(kb)]),
                                 np.log(db), rcond=None)
        boot.append(-cb[1])
    res['fit'] = dict(exponent=float(exp_hat),
                      lo=float(np.percentile(boot, 2.5)),
                      hi=float(np.percentile(boot, 97.5)))
    print(f"  exponent = -{exp_hat:.3f} [{-res['fit']['hi']:.3f}, {-res['fit']['lo']:.3f}] (95%)")

    print('== B) parameter scan (<k>=8) ==')
    for label, alpha, mu, N in [('mu0.5', 1.0, 0.5, 20000), ('mu1', 1.0, 1.0, 20000),
                                ('mu2', 1.0, 2.0, 20000), ('a0.5', 0.5, 0.5, 20000),
                                ('a3', 3.0, 0.5, 20000), ('N5e3', 1.0, 0.5, 5000)]:
        adj, z, n = build_graph('poisson', N, mean_k=8, rng=rng)
        eps = max(1, round(0.001 * n))
        per_real = [ledger_dev(adj, n, alpha, mu, eps, rng, runs=2)[0]
                    for _ in range(4)]
        res['scan'].append(dict(label=label, alpha=alpha, mu=mu, N=N, z=z,
                                dev_mean=float(np.mean(per_real)),
                                dev_sd=float(np.std(per_real))))
        print(f"  {label}: dev={np.mean(per_real):.2f}% ± {np.std(per_real):.2f}")

    json.dump(res, open(OUT, 'w'), indent=1)
    print('->', os.path.relpath(OUT, ROOT))


if __name__ == '__main__':
    main()
