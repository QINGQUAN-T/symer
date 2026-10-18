"""W3 共享库 v2：构型模型 CTMC + HMF 退火 ODE + well-mixed Gillespie.

约定（与 W1/W2 一致的无量纲化，网络采用"每边速率"制）：
  REC : S-I 边，速率 1/边                    -> I 端 -> S1
  ARM : 有向对 (u∈池 -> v=S1)，速率 alpha/对 -> v -> S2
  STIF: 有向对 (u∈池 -> v∈池)，速率 mu/对    -> v -> R（hearer=v 出传）
  与 MF 同轴对比时时间乘 <k>（每 spreader 接触率 = <k>·i 对齐 MF 的 k1·i）。
  一条 S1-S1 无向边产生 4 个事件项（双方 STIF + 双方 ARM）；S1-S2 产生 3 项；
  S2-S2 产生 2 项；S-I 产生 1 项（REC）。

α-盲可集总性（W3 定理的构造性体现）：
  REC/STIF/边型计数只依赖池 membership（S1/S2 对池完全对称），ARM 只改标签
  不动池与边型 —— 故 (i, S, r) 与全部边型计数过程对任意静态图严格 α-盲。
"""
import numpy as np

I, S1, S2, R = 0, 1, 2, 3
POOL = (S1, S2)


def build_graph(kind, N, mean_k=None, tau=None, rng=None, kmin=2):
    """返回 (adj list[np.array], 实际平均度, N_eff)。孤立点剔除。"""
    import networkx as nx
    rng = rng or np.random.default_rng()
    if kind == "poisson":
        seq = np.maximum(rng.poisson(mean_k, N), 0)
    elif kind == "sf":
        ks = np.arange(kmin, 10 * N)
        w = ks.astype(float) ** (-tau)
        seq = rng.choice(ks, size=N, p=w / w.sum())
    else:
        raise ValueError(kind)
    if seq.sum() % 2 == 1:
        seq = seq.copy()
        seq[-1] += 1
    G = nx.configuration_model(seq.tolist(), seed=int(rng.integers(1 << 30)))
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    G.remove_nodes_from(list(nx.isolates(G)))
    G = nx.convert_node_labels_to_integers(G)
    N_eff = G.number_of_nodes()
    adj = [np.fromiter(G.neighbors(v), dtype=np.int64, count=G.degree(v))
           for v in range(N_eff)]
    return adj, 2.0 * G.number_of_edges() / N_eff, N_eff


def build_graph_pmf(ks, p, N, rng):
    """按给定度分布 pmf 生成构型模型（截断 pmf，用于 SF+指数截断）。"""
    import networkx as nx
    seq = rng.choice(ks, size=N, p=p)
    if seq.sum() % 2 == 1:
        seq = seq.copy()
        seq[-1] += 1
    G = nx.configuration_model(seq.tolist(), seed=int(rng.integers(1 << 30)))
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    G.remove_nodes_from(list(nx.isolates(G)))
    G = nx.convert_node_labels_to_integers(G)
    N_eff = G.number_of_nodes()
    adj = [np.fromiter(G.neighbors(v), dtype=np.int64, count=G.degree(v))
           for v in range(N_eff)]
    return adj, 2.0 * G.number_of_edges() / N_eff, N_eff


class NetSim:
    """事件驱动 CTMC。桶 = swap-delete 数组 + 位置字典，O(1) 采样。"""

    def __init__(self, adj, alpha, mu, eps_count, rng):
        self.adj, self.alpha, self.mu, self.rng = adj, alpha, mu, rng
        self.N = len(adj)
        self.state = np.full(self.N, I, dtype=np.int8)
        self.ci, self.cs1, self.cs2 = self.N, 0, 0
        for v in self.rng.choice(self.N, size=eps_count, replace=False):
            self.state[v] = S1
        self.ci, self.cs1 = self.N - eps_count, eps_count
        self.bk = {"SI": [], "SS": [], "ARM": []}   # 编码见 _events_of
        self.bp = {k: {} for k in self.bk}
        for v in range(self.N):
            if self.state[v] != I:
                for w in self.adj[v]:
                    for key, code in self._events_of(v, w):
                        self._add(key, code)
        self.t = 0.0
        self.x1 = self.x2 = 0

    # ---- 事件枚举：有序视图 (u -> v)，v 为被作用端 ----
    def _events_of(self, u, v):
        su, sv = self.state[u], self.state[v]
        if su not in POOL:
            return ()
        if sv == I:
            return (("SI", u * self.N + v),)        # u 招募 v
        out = []
        if sv in POOL:
            out.append(("SS", v * self.N + u))      # hearer=v 出传
        if sv == S1:
            out.append(("ARM", v * self.N + u))     # v 被再激活
        return tuple(out)

    def _add(self, key, code):
        b, p = self.bk[key], self.bp[key]
        p[code] = len(b)
        b.append(code)

    def _remove(self, key, code):
        b, p = self.bk[key], self.bp[key]
        j = p.pop(code)
        last = b.pop()
        if last != code:
            b[j] = last
            p[last] = j

    def _set_state(self, v, new):
        old = self.state[v]
        if old == new:
            return
        for w in self.adj[v]:
            for key, code in self._events_of(v, w):
                self._remove(key, code)
            for key, code in self._events_of(w, v):
                self._remove(key, code)
        self.state[v] = new
        if old == I:
            self.ci -= 1
        if new == I:
            self.ci += 1
        if old == S1:
            self.cs1 -= 1
        if old == S2:
            self.cs2 -= 1
        if new == S1:
            self.cs1 += 1
        if new == S2:
            self.cs2 += 1
        if new == R:
            if old == S1:
                self.x1 += 1
            elif old == S2:
                self.x2 += 1
        for w in self.adj[v]:
            for key, code in self._events_of(v, w):
                self._add(key, code)
            for key, code in self._events_of(w, v):
                self._add(key, code)

    def _step(self):
        n_si, n_ss, n_arm = (len(self.bk["SI"]), len(self.bk["SS"]),
                             len(self.bk["ARM"]))
        a_rec, a_arm, a_stif = 1.0 * n_si, self.alpha * n_arm, self.mu * n_ss
        a = a_rec + a_arm + a_stif
        if a <= 0:
            return np.inf, None
        tau = self.rng.exponential(1.0 / a)
        u = self.rng.random() * a
        if u < a_rec:
            code = self.bk["SI"][self.rng.integers(n_si)]
            v = int(code % self.N)          # I 端
            return tau, (v, S1)
        if u < a_rec + a_arm:
            code = self.bk["ARM"][self.rng.integers(n_arm)]
            return tau, (int(code // self.N), S2)
        code = self.bk["SS"][self.rng.integers(n_ss)]
        return tau, (int(code // self.N), R)

    def run(self, t_max=np.inf, t_grid=None):
        """返回 dict 终态；t_grid 给定时同时返回 (i, S) 分段常数轨迹。"""
        tg = list(t_grid) if t_grid is not None else None
        ki = 0
        if tg is not None:
            out_i = np.empty(len(tg))
            out_S = np.empty(len(tg))
        while True:
            tau, act = self._step()
            t_next = self.t + tau
            if tg is not None:
                while ki < len(tg) and tg[ki] < min(t_next, t_max):
                    out_i[ki] = self.ci / self.N
                    out_S[ki] = (self.cs1 + self.cs2) / self.N
                    ki += 1
            if t_next > t_max or not np.isfinite(t_next):
                self.t = min(self.t + (tau if np.isfinite(tau) else 0), t_max)
                if tg is not None:
                    while ki < len(tg):
                        out_i[ki] = self.ci / self.N
                        out_S[ki] = (self.cs1 + self.cs2) / self.N
                        ki += 1
                break
            self.t = t_next
            v, new = act
            self._set_state(v, new)
        res = dict(t=self.t, i_f=self.ci, x1=self.x1, x2=self.x2)
        if tg is not None:
            res["i_traj"], res["S_traj"] = out_i, out_S
        return res


def simulate_network(adj, alpha, mu, eps_count, rng, t_grid=None, t_max=np.inf):
    return NetSim(adj, alpha, mu, eps_count, rng).run(t_max=t_max,
                                                      t_grid=t_grid)


# ----------------------------------------------------------------------
# HMF 退火 ODE（k 类；每边速率制，与网络同轴）
# ----------------------------------------------------------------------
def hmf_final(p_k, alpha, mu, t_max=5e4, n_store=3000):
    """返回 dict(i_f, x1, x2, t, S)。p_k 为类分布向量（下标=度）。"""
    from scipy.integrate import solve_ivp
    ks = np.nonzero(p_k > 1e-16)[0]
    pk = p_k[ks]
    nc = len(ks)
    kavg = float((ks * pk).sum())

    def rhs(t, y):
        ik, s1k, s2k = y[:nc], y[nc:2 * nc], y[2 * nc:]
        theta = float((ks * pk * (s1k + s2k)).sum()) / kavg
        return np.concatenate([
            -ks * ik * theta,
            ks * ik * theta - (alpha + mu) * ks * s1k * theta,
            alpha * ks * s1k * theta - mu * ks * s2k * theta])

    y0 = np.concatenate([np.full(nc, 1.0 - 1e-3), np.full(nc, 1e-3),
                         np.zeros(nc)])
    ts = np.linspace(0, t_max, n_store)
    sol = solve_ivp(rhs, (0, t_max), y0, method="LSODA", rtol=1e-9,
                    atol=1e-13, t_eval=ts)
    y = sol.y[:, -1]
    ys = sol.y
    s1t, s2t = ys[nc:2 * nc], ys[2 * nc:]
    Sk = s1t + s2t
    theta = ((ks * pk) @ Sk) / kavg
    x1 = mu * float(np.trapezoid(
        (s1t * (ks * pk)[:, None]).sum(axis=0) * theta, ts))
    x2 = mu * float(np.trapezoid(
        (s2t * (ks * pk)[:, None]).sum(axis=0) * theta, ts))
    return dict(i_f=float((pk * y[:nc]).sum()), x1=x1, x2=x2, t=ts,
                S=pk @ Sk)


def poisson_pmf(mean_k, kmax=200):
    from scipy.stats import poisson
    ks = np.arange(kmax + 1)
    p = poisson.pmf(ks, mean_k)
    return p / p.sum()


def sf_pmf(tau, N, kmin=2):
    ks = np.arange(kmin, 10 * N)
    w = ks.astype(float) ** (-tau)
    p = np.zeros(10 * N)
    p[ks] = w / w.sum()
    return p


# ----------------------------------------------------------------------
# well-mixed Gillespie（W1/W2 质量作用约定，供有限 N 相图）
# ----------------------------------------------------------------------
def gillespie(N, alpha, mu, eps_count, rng):
    ic, s1, s2 = N - eps_count, eps_count, 0
    x1 = x2 = 0
    while s1 + s2 > 0:
        S = s1 + s2
        a1 = 1.0 * ic * S / N
        a2 = alpha * 1.0 * s1 * (S - 1) / N
        a3 = mu * s1 * (S - 1) / N
        a4 = mu * s2 * (S - 1) / N
        a = a1 + a2 + a3 + a4
        if a <= 0:
            break
        u = rng.random() * a
        if u < a1:
            ic -= 1; s1 += 1
        elif u < a1 + a2:
            s1 -= 1; s2 += 1
        elif u < a1 + a2 + a3:
            s1 -= 1; x1 += 1
        else:
            s2 -= 1; x2 += 1
    return ic, x1, x2


# ----------------------------------------------------------------------
# 自检：网络引擎的不变量与已知极限
# ----------------------------------------------------------------------
def selftest():
    rng = np.random.default_rng(7)
    adj, z, n = build_graph("poisson", 2000, mean_k=8, rng=rng)
    sim = NetSim(adj, 1.0, 0.5, 2, rng)
    assert len(sim.bk["SI"]) == sum(
        1 for v in range(n) if sim.state[v] != I
        for w in adj[v] if sim.state[w] == I), "SI 计数不符"
    assert len(sim.bk["SS"]) == 2 * sum(
        1 for v in range(n) if sim.state[v] in POOL
        for w in adj[v] if w > v and sim.state[w] in POOL), "SS 计数不符"
    r = sim.run()
    assert sim.cs1 + sim.cs2 + sim.ci + (sim.x1 + sim.x2) == n, "守恒破坏"
    assert sim.cs1 + sim.cs2 == 0 or sim.t > 1e6, "未吸收"
    print(f"selftest OK: N={n}, <k>={z:.2f}, t={r['t']:.1f}, "
          f"i_f={r['i_f']}, x1={r['x1']}, x2={r['x2']}")


if __name__ == "__main__":
    selftest()
