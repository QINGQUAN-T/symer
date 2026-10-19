# -*- coding: utf-8 -*-
"""W7 补充：模拟档案种子间变异扫描（审稿 v5 意见 C15/C16）。

对主网格 (mu, alpha) × {alpha} 在匹配阈值 theta=2（k1=1 单位）用 10 个独立
生成器种子重跑，每 (seed, branch) 60 个档案；报告每分支跨种子的
alpha_hat 中位数与 [min, max]，以及 alpha=0 分支的假阳性率。

输出: results/w7_seed_variation.json（每 seed 增量写入，可中断续用）
"""
import json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w7_sim_archives import run_cell

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = ROOT + '/results/w7_seed_variation.json'
MUS = [0.5, 1.37, 2.0]
ALPHAS = [0.0, 0.005, 0.01, 0.02, 0.04, 0.08]
THETA = 2.0
SEEDS = list(range(20261, 20271))
N_REP = 60

def load():
    if os.path.exists(OUT):
        return json.load(open(OUT))
    return {'theta': THETA, 'n_rep': N_REP, 'seeds_done': [], 'cells': {}}

def main():
    res = load()
    for seed in SEEDS:
        if seed in res['seeds_done']:
            continue
        rng = np.random.default_rng(seed)
        for mu in MUS:
            for al in ALPHAS:
                c = run_cell(mu, al, THETA, N_REP, rng)
                key = f'mu={mu}|a={al}'
                res['cells'].setdefault(key, []).append(
                    dict(seed=seed, ahat_med=c['ahat_med'], n_ok=c['n_ok'],
                         intercept_med=c['intercept_med']))
                print(f"seed={seed} mu={mu} a={al}: ahat_med={c['ahat_med']:.4f} n_ok={c['n_ok']}",
                      flush=True)
        res['seeds_done'].append(seed)
        # 增量汇总
        summary = {}
        for key, rows in res['cells'].items():
            meds = [r['ahat_med'] for r in rows if r['ahat_med'] is not None]
            if meds:
                summary[key] = dict(n_seeds=len(meds),
                                    med=float(np.median(meds)),
                                    lo=float(np.min(meds)), hi=float(np.max(meds)))
        res['summary'] = summary
        json.dump(res, open(OUT, 'w'), indent=1)
        print(f'-- seed {seed} done ({len(res["seeds_done"])}/{len(SEEDS)}) --', flush=True)
    print('ALL DONE', flush=True)

if __name__ == '__main__':
    t0 = time.time()
    main()
    print(f'total {time.time()-t0:.0f}s', flush=True)
