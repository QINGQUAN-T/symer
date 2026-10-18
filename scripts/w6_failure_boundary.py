# -*- coding: utf-8 -*-
"""W6 图级源数据存档（审稿意见 M6.4）: 论文图 6（网络失效边界）。

(a) 触达失效边界: Poisson <k>∈{4,6,8,12,20,40}, alpha=1, mu=0.5, N=2e4,
    8 SSA 种子/点; R_f/N = (x1+x2)/N 与 HMF（w3_lib.hmf_final, 同度分布）;
    完全混合极限 R_f=1-eps 处的水平线。
(b) 账本偏差 vs <k>（log-log）: Poisson 同扫描 + scale-free（tau=2.5, 3.0,
    匹配 <k>=8）diamonds。
种子协议: 本图 8 种子/点（图 7 为 10 图实现×3 种子，见 App. D 种子说明）。
输出: results/w6_failure_boundary_data.json
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w3_lib import NetSim, build_graph, hmf_final, poisson_pmf

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N = 20000
ALPHA, MU = 1.0, 0.5
KS = (4, 6, 8, 12, 20, 40)
N_SEEDS = 8


def reach_run(adj, n, rng):
    sim = NetSim(adj, ALPHA, MU, max(1, n // 1000), rng)
    r = sim.run()
    return (r['x1'] + r['x2']) / n, r['x2'] / max(1, r['x1'] + r['x2'])


def main():
    rng = np.random.default_rng(20261018)
    out = {'meta': dict(N=N, alpha=ALPHA, mu=MU, ks=list(KS), n_seeds=N_SEEDS,
                        seed=20261018),
           'reach': [], 'ledger_dev_poisson': [], 'ledger_dev_sf': [],
           'hmf': []}
    law_share = ALPHA / (ALPHA + MU)
    for zt in KS:
        rf, devs = [], []
        for s in range(N_SEEDS):
            adj, z, n = build_graph('poisson', N, mean_k=zt, rng=rng)
            r, share = reach_run(adj, n, rng)
            rf.append(r)
            devs.append((share - law_share) / law_share * 100)
        h = hmf_final(poisson_pmf(zt), ALPHA, MU)
        h_share = h['x2'] / max(1e-12, h['x1'] + h['x2'])
        out['reach'].append(dict(k=zt, rf_N_mean=float(np.mean(rf)),
                                 rf_N_sd=float(np.std(rf)),
                                 rf_N_all=[float(x) for x in rf],
                                 hmf_reach=float(h['x1'] + h['x2']),
                                 hmf_i_f=float(h['i_f']),
                                 hmf_share=float(h_share),
                                 hmf_share_dev_pct=float(
                                     (h_share - law_share) / law_share * 100)))
        out['ledger_dev_poisson'].append(dict(
            k=zt, dev_mean=float(np.mean(devs)), dev_sd=float(np.std(devs)),
            dev_all=[float(x) for x in devs]))
        print(f"<k>={zt}: Rf/N={np.mean(rf):.4f}±{np.std(rf):.4f} "
              f"HMF_reach={h['x1'] + h['x2']:.4f} dev={np.mean(devs):.2f}%±{np.std(devs):.2f}")
    # scale-free, 匹配 <k>=8
    K_TARGET = 8.0
    for tau in (2.5, 3.0):
        # 匹配 <k>≈8: 纯幂律矩 mean = kmin*(tau-1)/(tau-2) -> kmin
        kmin = max(2, int(round(K_TARGET * (tau - 2) / (tau - 1))))
        devs, zs, emp_pmf = [], [], None
        for s in range(N_SEEDS):
            adj, z, n = build_graph('sf', N, tau=tau, rng=rng, kmin=kmin)
            r, share = reach_run(adj, n, rng)
            devs.append((share - law_share) / law_share * 100)
            zs.append(z)
            # App. D: scale-free 的 HMF 用实证度直方图（理论 pmf 重尾使
            # ODE 类数爆炸）; 取最后一个图实现的度分布
            emp_pmf = np.bincount([len(a) for a in adj],
                                  minlength=max(len(a) for a in adj) + 1
                                  ).astype(float) / n
        h = hmf_final(emp_pmf, ALPHA, MU)
        h_share = h['x2'] / max(1e-12, h['x1'] + h['x2'])
        out['ledger_dev_sf'].append(dict(tau=tau, k_matched=float(np.mean(zs)),
                                         dev_mean=float(np.mean(devs)),
                                         dev_sd=float(np.std(devs)),
                                         dev_all=[float(x) for x in devs],
                                         hmf_share=float(h_share),
                                         hmf_dev_pct=float(
                                             (h_share - law_share) / law_share * 100)))
        print(f"sf tau={tau}: <k>={np.mean(zs):.2f} dev={np.mean(devs):.2f}%")
    json.dump(out, open(f'{ROOT}/results/w6_failure_boundary_data.json', 'w'),
              indent=1)
    print('-> results/w6_failure_boundary_data.json')


if __name__ == '__main__':
    main()
