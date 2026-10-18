# -*- coding: utf-8 -*-
"""W7-R2 负控制与稳健性（审稿意见 #8）。真实 Higgs-Twitter RT 数据，θ=6h。

  C0 基线       : 原样管线（与 W5 演示同参数）
  C1 时间戳打乱 : 级联内事件时间随机重排（破坏用户-时间线耦合）
  C2 源账户打乱 : 事件重新随机指派到等尺寸的随机源（破坏级联内聚性）
  C4 保活跃度打乱 : 级联内各用户事件数保持不变、时间戳在级联内重排
                  （保留用户活跃度异质性，只破坏用户-时间耦合；审稿意见
                   Moderate #5 的中度对照。C0-C2 顺序与种子不变，结果与
                   v0.1.1 逐位一致）
输出: results/w7_negctl_results.json（v0.1.2：路径相对化，种子固定）
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w5_invert import load_rt, analyze_cascade

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/results/w7_negctl_results.json'
DATA = f'{ROOT}/data/higgs-activity.txt.gz'
THETA = 6.0 * 3600.0          # 6h，秒


def shuffle_activity_preserving(evs, rng):
    """C4: 保每用户事件数、级联内重排时间戳（边际时间线与活跃度分布不变）。"""
    users = {}
    for u, t in evs:
        users.setdefault(u, []).append(t)
    us = list(users.keys())
    counts = [len(users[u]) for u in us]
    ts_all = sorted(t for _, t in evs)
    perm = rng.permutation(len(ts_all))
    ts_sh = [ts_all[i] for i in perm]
    new_evs, pos = [], 0
    for u, c in zip(us, counts):
        for t in ts_sh[pos:pos + c]:
            new_evs.append((u, t))
        pos += c
    return new_evs


def summarize(casc, rng, shuffle_time=False, shuffle_source=False,
              shuffle_activity=False):
    if shuffle_source:
        keys = list(casc.keys())
        sizes = [len(casc[v]) for v in keys]
        all_ev = [e for evs in casc.values() for e in evs]
        order = rng.permutation(len(all_ev))
        new, pos = {}, 0
        rng.shuffle(keys)
        for v, n in zip(keys, sizes):
            new[v] = [all_ev[i] for i in order[pos:pos + n]]
            pos += n
        casc = new
    rows = []
    for v, evs in casc.items():
        if len(evs) < 30:
            continue
        ts = [t for _, t in evs]
        span = max(ts) - min(ts)
        if span < 2 * THETA:
            continue
        if max(ts) > T_END - 24 * 3600.0:
            continue                      # 删失保护：级联已休眠（与 W5 相同）
        ev = [(u, t) for u, t in evs]
        if shuffle_activity:
            ev = shuffle_activity_preserving(ev, rng)
        if shuffle_time:
            ts_sh = list(rng.permutation(ts))
            ev = [(u, t2) for (u, _), t2 in zip(ev, ts_sh)]
        c = analyze_cascade(ev, THETA, 16)
        if c is not None:
            rows.append(c)
    a = np.array([r['alpha'] for r in rows]) if rows else np.array([np.nan])
    return dict(n=len(rows),
                mu_med=float(np.nanmedian([r['mu'] for r in rows])) if rows else None,
                ic_med=float(np.nanmedian([r['intercept'] for r in rows])) if rows else None,
                r2_med=float(np.nanmedian([r['r2'] for r in rows])) if rows else None,
                ahat_med=float(np.nanmedian(a)),
                ahat_lo=float(np.nanpercentile(a, 25)) if rows else None,
                ahat_hi=float(np.nanpercentile(a, 75)) if rows else None,
                x2_pooled=(sum(r['n_bursts'][2] for r in rows) /
                           max(1, sum(r['stifled'] for r in rows))) if rows else None)


def main():
    casc = load_rt(DATA, 'RT')
    global T_END
    T_END = max(t for evs in casc.values() for _, t in evs)
    casc = {v: evs for v, evs in casc.items() if len(evs) >= 30}
    print('cascades >=30 users:', len(casc))
    rng = np.random.default_rng(20261013)
    res = {}
    for label, kw in [('C0_baseline', {}),
                      ('C1_time_shuffle', dict(shuffle_time=True)),
                      ('C2_source_shuffle', dict(shuffle_source=True)),
                      ('C4_activity_shuffle', dict(shuffle_activity=True))]:
        s = summarize(casc, rng, **kw)
        res[label] = s
        print(label, json.dumps(s, ensure_ascii=False))
    json.dump(res, open(OUT, 'w'), indent=1)
    print('->', os.path.relpath(OUT, ROOT))


if __name__ == '__main__':
    main()
