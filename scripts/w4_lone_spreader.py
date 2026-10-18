# -*- coding: utf-8 -*-
"""W4 图级源数据存档（审稿意见 M6.4）: 论文图 8（孤传播者，定理 8）。

36 个 (N, μ) 构型（N∈{1e2,3e2,1e3,3e3,1e4,3e4} × μ∈{0.3,0.5,1,2,3,5}, α=0）:
  - P(i_f>0) 的测量值与 95% 上界（rule of three: 3/M）
  - 终态 S_f = 1（E[S_f] 精确为 1）
  - 账本修正 1 − R_f/N ≈ 1/N
批量向量化 lockstep Gillespie（与 w3_lib.gillespie 逐事件同分布：每副本
独立直接法, 仅在副本间交错推进；随机流 seed=20261018）。

输出: results/w4_lone_spreader_data.json
"""
import json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NS = (100, 300, 1000, 3000, 10000, 30000)
MUS = (0.3, 0.5, 1.0, 2.0, 3.0, 5.0)
ALPHA = 0.0


def batch_gillespie(N, mu, eps_count, M, rng):
    """M 个独立副本的 lockstep 直接法 SSA。返回终态数组 (i_f, x1, x2, S_f)。

    与 w3_lib.gillespie 的单副本事件率逐项一致:
      a1 = ic*S/N; a2 = alpha*s1*(S-1)/N; a3 = mu*s1*(S-1)/N; a4 = mu*s2*(S-1)/N。
    lockstep 仅改变副本间的推进交错，不改变每副本的事件分布。"""
    ic = np.full(M, N - eps_count, dtype=np.int64)
    s1 = np.full(M, eps_count, dtype=np.int64)
    s2 = np.zeros(M, dtype=np.int64)
    x1 = np.zeros(M, dtype=np.int64)
    x2 = np.zeros(M, dtype=np.int64)
    alive = np.ones(M, dtype=bool)
    while alive.any():
        idx = np.nonzero(alive)[0]
        S = (s1[idx] + s2[idx]).astype(np.float64)
        icf = ic[idx].astype(np.float64)
        s1f = s1[idx].astype(np.float64)
        s2f = s2[idx].astype(np.float64)
        a1 = icf * S / N
        a2 = ALPHA * s1f * np.maximum(S - 1.0, 0.0) / N
        a3 = mu * s1f * np.maximum(S - 1.0, 0.0) / N
        a4 = mu * s2f * np.maximum(S - 1.0, 0.0) / N
        a = a1 + a2 + a3 + a4
        # a<=0 的副本已吸收（S_f=1: 末位传播者无 stifling 通道）→ 移出
        stuck = a <= 1e-12
        if stuck.any():
            alive[idx[stuck]] = False
        act = ~stuck
        idx_a = idx[act]
        if len(idx_a) == 0:
            continue
        ua = rng.random(len(idx_a)) * a[act]
        a1a, a2a, a3a = a1[act], a2[act], a3[act]
        rec = ua < a1a
        arm = ~rec & (ua < a1a + a2a)
        st1 = ~rec & ~arm & (ua < a1a + a2a + a3a)
        st2 = ~rec & ~arm & ~st1
        i_rec = idx_a[rec]; i_arm = idx_a[arm]; i_s1 = idx_a[st1]; i_s2 = idx_a[st2]
        ic[i_rec] -= 1; s1[i_rec] += 1
        s1[i_arm] -= 1; s2[i_arm] += 1
        s1[i_s1] -= 1; x1[i_s1] += 1
        s2[i_s2] -= 1; x2[i_s2] += 1
        alive[(s1 + s2) <= 0] = False
    return ic, x1, x2, s1 + s2


def main():
    rng = np.random.default_rng(20261018)
    out = {'meta': {'engine': 'lockstep batch direct-method SSA (per-replica rates '
                    'identical to w3_lib.gillespie)', 'alpha': ALPHA,
                    'seed': 20261018,
                    'replicates': 'M=4000 for N<=3000; 1000 for N=1e4; 400 for N=3e4'},
           'cells': []}
    for N in NS:
        M = 4000 if N <= 3000 else (1000 if N == 10000 else 400)
        eps = max(1, N // 1000)
        for mu in MUS:
            ic, x1, x2, sf = batch_gillespie(N, mu, eps, M, rng)
            p_if = float(np.mean(ic > 0))
            ub = 3.0 / M                        # rule of three
            one_minus_Rf = 1.0 - (x1 + x2) / N  # 账本修正 ≈ 1/N
            out['cells'].append(dict(
                N=N, mu=mu, M=M, eps=eps,
                P_if_pos=p_if, upper_bound_95=ub,
                E_Sf=float(np.mean(sf)), max_Sf=int(np.max(sf)),
                one_minus_Rf=float(np.mean(one_minus_Rf)), one_over_N=1.0 / N))
            c = out['cells'][-1]
            print(f"N={N:>6} mu={mu:>3}: P(i_f>0)={p_if} ub95={ub:.4f} "
                  f"E[Sf]={c['E_Sf']:.4f} 1-Rf/N={c['one_minus_Rf']:.5f} vs 1/N={1/N:.5f}")
    json.dump(out, open(f'{ROOT}/results/w4_lone_spreader_data.json', 'w'), indent=1)
    print('-> results/w4_lone_spreader_data.json')


if __name__ == '__main__':
    main()
