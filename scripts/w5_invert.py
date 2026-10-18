"""W5 数据演示：真实转发级联上的 x1 反演 CLI.

模型 ↔ 档案映射（W4 论文 §VII）：
  级联        = 同一源用户 v 的全部转推事件 (u, v, t)   [Higgs activity_time, RT 行]
  招募        = 用户 u 对该源内容的首次转推             -> i 减一, S 加一
  再激活(ARM) = u 在间隔 > theta 后再次转推同源内容     -> r 减一, s2 加一
  出传(R)     = u 的末次转推（a posteriori 判定）        -> S 减一
  孤传播者    = 级联末次事件的那个用户（命运删失）       -> 账本池 = P - 1（定理 6 的 −1/N）
免拟合检验：轨迹回归 ΔS/Δi = a + b·(S/i)（y 取 −dS/di 约定），截距 a 理论值 +1
（合并层 dS/di = −1 + μS/i），斜率 b = −μ；反演：mu_hat = −b；
alpha_hat = mu_hat * x2/x1（x1,x2 在"已出传"池上计）。
误差带：级联为单元的层级 bootstrap。
"""
import argparse
import gzip
import json

import numpy as np


def load_rt(path, rtype="RT"):
    """读 activity_time 文件，返回 {source: [(actor, t), ...]}（仅 RT）。"""
    casc = {}
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            p = ln.split()
            if len(p) != 4 or p[3] != rtype:
                continue
            u, v, t = int(p[0]), int(p[1]), int(p[2])
            casc.setdefault(v, []).append((u, t))
    return casc


def session_bursts(times, theta):
    """间隔 > theta 记为新一轮爆发。"""
    ts = np.sort(np.asarray(times))
    return 1 + int((np.diff(ts) > theta).sum())


def analyze_cascade(events, theta, n_bins):
    """单级联：账本 + 轨迹回归。返回 dict 或 None。"""
    events = sorted(events, key=lambda e: e[1])
    users = {}
    for u, t in events:
        users.setdefault(u, []).append(t)
    P = len(users)
    if P < 5:
        return None
    last = {u: max(ts) for u, ts in users.items()}
    first = {u: min(ts) for u, ts in users.items()}
    # 孤传播者：末次事件最晚者命运删失，账本池 = P-1
    immortal = max(users, key=lambda u: (last[u], first[u]))
    ledger = {1: 0, 2: 0}
    for u, ts in users.items():
        if u == immortal:
            continue
        ledger[1 if session_bursts(ts, theta) == 1 else 2] += 1
    stifled = ledger[1] + ledger[2]
    if stifled < 5 or ledger[1] == 0:
        return None
    # a posteriori 轨迹（用户级计数）：i(t)=P-已招募; S(t)=first<=t<last
    fs = np.array([first[u] for u in users])
    ls = np.array([last[u] for u in users]) + theta  # 停留窗口 tau=theta：末事件后仍可被接触
    t_lo, t_hi = min(e[1] for e in events), max(e[1] for e in events)
    edges = np.linspace(t_lo, t_hi + 1e-9, n_bins + 1)
    pts_x, pts_y = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        i0 = int(P - (fs <= a).sum())          # a 时刻未招募
        i1 = int(P - (fs < b).sum())           # b 时刻未招募
        s0 = int(((fs <= a) & (ls > a)).sum())  # a 时刻活跃
        s1 = int(((fs <= b) & (ls > b)).sum())  # b 时刻活跃
        di, dS = i0 - i1, s1 - s0
        if di < 2:
            continue
        x_mid = ((s0 + s1) / 2.0) / max((i0 + i1) / 2.0, 1e-9)
        pts_x.append(x_mid)
        pts_y.append(dS / di)                  # = -dS/di = 1 - mu*(S/i)
    n_pts = len(pts_y)
    if n_pts < 6:
        return None
    x = np.array(pts_x)
    y = np.array(pts_y)
    X = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    resid = y - X @ coef
    r2 = 1 - (resid @ resid) / ((y - y.mean()) @ (y - y.mean()) + 1e-300)
    if b >= 0:                                  # mu_hat = -b 须为正
        return None
    mu = -b
    x1f, x2f = ledger[1] / stifled, ledger[2] / stifled
    return dict(n_users=P, stifled=stifled, n_bursts=ledger, n_pts=n_pts,
                intercept=a, mu=mu, r2=float(r2), x1=x1f, x2=x2f,
                alpha=mu * x2f / x1f, x_mean=float(np.mean(x)),
                px=x.tolist(), py=y.tolist())


def pooled(casc_list, boot, rng):
    mus = np.array([c["mu"] for c in casc_list])
    alphas = np.array([c["alpha"] for c in casc_list])
    x2s = np.array([c["x2"] for c in casc_list])
    outs = {}
    for name, arr in (("mu", mus), ("alpha", alphas), ("x2", x2s)):
        bs = np.array([np.median(rng.choice(arr, len(arr))) for _ in range(boot)])
        outs[name] = dict(median=float(np.median(arr)),
                          lo=float(np.percentile(bs, 2.5)),
                          hi=float(np.percentile(bs, 97.5)),
                          n=int(len(arr)))
    return outs


def main():
    ap = argparse.ArgumentParser(description="x1 反演 CLI（转发级联档案）")
    ap.add_argument("--activity", required=True, help="Higgs activity_time(.gz)：u v t TYPE")
    ap.add_argument("--rtype", default="RT", choices=["RT", "RE", "MT"])
    ap.add_argument("--out", default="w5_demo")
    ap.add_argument("--theta-hours", type=float, nargs="+", default=[6.0],
                    help="会话化间隔（小时），可多个做敏感性")
    ap.add_argument("--min-users", type=int, default=30)
    ap.add_argument("--bins", type=int, default=16)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261009)
    ap.add_argument("--dataset", default="higgs-twitter (SNAP)")
    a = ap.parse_args()

    casc = load_rt(a.activity, a.rtype)
    t_end = max(t for evs in casc.values() for _, t in evs)
    rng = np.random.default_rng(a.seed)
    results = {"dataset": a.dataset, "rt_events": sum(len(v) for v in casc.values()),
               "sources": len(casc), "t_end_utc": int(t_end),
               "dormancy_hours": 24.0, "min_users": a.min_users,
               "bins": a.bins, "boot": a.boot, "theta": []}
    for th_h in a.theta_hours:
        th = th_h * 3600.0
        rows = []
        for v, evs in casc.items():
            if len(evs) < a.min_users:
                continue
            span = max(t for _, t in evs) - min(t for _, t in evs)
            if span < 2 * th:
                continue
            if max(t for _, t in evs) > t_end - 24 * 3600.0:
                continue  # 删失保护：级联已休眠
            c = analyze_cascade(evs, th, a.bins)
            if c is not None:
                c["source"] = int(v)
                rows.append(c)
        summ = pooled(rows, a.boot, np.random.default_rng(a.seed)) if rows else {}
        results["theta"].append(dict(theta_h=th_h, n_cascades=len(rows),
                                     pooled=summ, rows=rows))
        print(f"theta={th_h}h: 级联 {len(rows)} | mu_med={summ['mu']['median']:.3f} "
              f"[{summ['mu']['lo']:.3f},{summ['mu']['hi']:.3f}] | "
              f"alpha_med={summ['alpha']['median']:.3f} "
              f"[{summ['alpha']['lo']:.3f},{summ['alpha']['hi']:.3f}] | "
              f"x2={summ['x2']['median']:.3f}")
    json.dump(results, open(a.out + "_results.json", "w"), indent=1)
    print("->", a.out + "_results.json")


if __name__ == "__main__":
    main()
