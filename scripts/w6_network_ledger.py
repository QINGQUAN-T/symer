# -*- coding: utf-8 -*-
"""W6 图级源数据存档（审稿意见 M6.4）: 论文图 5（网络双账本二分性）。

Poisson 配置模型 <k>=8, N=2e4, mu=0.5, alpha∈{0.5,3}:
  - 相同图实现上比较两个 alpha 的 E[S_t] 轨迹（合并层 α 盲, 定理 6）:
    95% CI 带的重叠时间占比;
  - 终态无知者占比 i_f/N 两组对比;
  - 职业层 x2 份额 vs 均值场律 alpha/(alpha+mu)（偏差 = 网络失效）。
种子: 4 图实现 × 4 SSA 种子 = 16 runs/alpha（同一组图跨 alpha 复用）。
输出: results/w6_network_ledger_data.json
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w3_lib import NetSim, build_graph

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N = 20000
MEAN_K = 8
MU = 0.5
ALPHAS = (0.5, 3.0)
N_GRAPH, N_SEED = 4, 4
T_GRID = np.geomspace(0.5, 3e4, 160)


def run_traj(adj, n, alpha, rng, t_grid):
    sim = NetSim(adj, alpha, MU, max(1, n // 1000), rng)
    return sim.run(t_grid=t_grid)


def main():
    rng = np.random.default_rng(20261018)
    out = {'meta': dict(N=N, mean_k=MEAN_K, mu=MU, alphas=list(ALPHAS),
                        n_graph=N_GRAPH, n_seed=N_SEED, seed=20261018,
                        engine='w3_lib.NetSim event-driven CTMC'),
           'graphs': [], 'traj': {}, 'terminal': {}}
    for g in range(N_GRAPH):
        adj, z, n = build_graph('poisson', N, mean_k=MEAN_K, rng=rng)
        out['graphs'].append(dict(graph_id=g, mean_k_actual=z, N_eff=n))
        for a in ALPHAS:
            for s in range(N_SEED):
                r = run_traj(adj, n, a, rng, T_GRID)
                key = f'a{a:g}'
                out['traj'].setdefault(key, []).append(
                    dict(graph_id=g, seed=s, t=T_GRID.tolist(),
                         S=r['S_traj'].tolist(), i=r['i_traj'].tolist()))
                out['terminal'].setdefault(key, []).append(
                    dict(graph_id=g, seed=s, i_f=r['i_f'], x1=r['x1'], x2=r['x2']))
                print(f"graph{g} seed{s} alpha={a}: i_f/N={r['i_f']/n:.4f} "
                      f"x2_share={r['x2']/max(1,r['x1']+r['x2']):.4f}")
    # 汇总: CI 带重叠 + 终态统计
    law = {a: a / (a + MU) for a in ALPHAS}
    summ = {}
    for a in ALPHAS:
        term = out['terminal'][f'a{a:g}']
        n_eff = [g['N_eff'] for g in out['graphs']]
        ifs = np.array([t['i_f'] / n_eff[t['graph_id']] for t in term])
        shares = np.array([t['x2'] / max(1, t['x1'] + t['x2']) for t in term])
        summ[f'a{a:g}'] = dict(
            i_f_N_med=float(np.median(ifs)),
            x2_share_med=float(np.median(shares)),
            x2_share_mean=float(shares.mean()),
            law_share=law[a],
            x2_share_dev_pct=float((shares.mean() - law[a]) / law[a] * 100))
    # CI 带重叠: 两个 alpha 的 E[S_t] 的 95% 带重叠网格点占比
    S_mean, S_lo, S_hi = {}, {}, {}
    for a in ALPHAS:
        arr = np.array([r['S'] for r in out['traj'][f'a{a:g}']])
        S_mean[a] = arr.mean(0)
        S_lo[a] = np.percentile(arr, 2.5, axis=0)
        S_hi[a] = np.percentile(arr, 97.5, axis=0)
    lo = np.maximum(S_lo[0.5], S_lo[3.0])
    hi = np.minimum(S_hi[0.5], S_hi[3.0])
    overlap = float(np.mean(hi >= lo))
    summ['ci_overlap_frac'] = overlap
    out['summary'] = summ
    json.dump(out, open(f'{ROOT}/results/w6_network_ledger_data.json', 'w'), indent=1)
    print('-> results/w6_network_ledger_data.json')
    print(json.dumps(summ, indent=1))


if __name__ == '__main__':
    main()
