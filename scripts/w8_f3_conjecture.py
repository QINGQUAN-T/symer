# -*- coding: utf-8 -*-
"""W8 f≥3 阶账本猜想的数值验证（审稿意见 Moderate #7）。

f=3 阶链: I -> S1 -> S2 -> S3 -> R; 速率(k1=1):
  招募 i S; 武装 α1 s1 S (S1->S2); 武装 α2 s2 S (S2->S3); 出传 μ(s1+s2+s3) S。
终端账本 x_m = ∫ μ s_m/i di（沿 i→0 积分）, Σx_m = 1。

猜想 B（嵌套伯努利链; f=2 时精确退化回定理 2）:
  x_m = [∏_{j<m} α_j/(μ+α_j)] · μ/(μ+α_m),  m=1..f (α_f≡0: 末级恒出传)。
  等价地: x1 = μ/(μ+α1); x2/(1−x1) = μ/(μ+α2); x3/(1−x1−x2) = α2/(μ+α2)。
  （注: 直觉推广 A, 即 (x1,x2,x3) ∝ (μ²,μα1,α1α2), 已被数值否定并保留在
   JSON 中作对照。）

输出: results/w8_f3_conjecture.json
"""
import json, os, sys
import numpy as np
from scipy.integrate import solve_ivp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EPS, I0 = 1e-4, 1.0 - 1e-4


def ledger(mu, a1, a2):
    """积分 f=3 链, 返回归一化终端账本 (x1,x2,x3)/Σ。"""
    def rhs(t, y):
        i, s1, s2, s3, X1, X2, X3 = y
        S = s1 + s2 + s3
        return [-i * S,
                i * S - (mu + a1) * s1 * S,
                a1 * s1 * S - (mu + a2) * s2 * S,
                a2 * s2 * S - mu * s3 * S,
                mu * s1 * S, mu * s2 * S, mu * s3 * S]
    y0 = [I0, EPS, 0, 0, 0, 0, 0]
    # 两段积分: 幂律尾衰减慢, 用账本增量判据确认饱和
    X_old, X, y = None, None, None
    t_end = 1e6
    for t_end in (1e6, 1e8, 1e10):
        sol = solve_ivp(rhs, (0, t_end), y0, method='LSODA', rtol=1e-11,
                        atol=1e-15)
        y = sol.y[:, -1]
        X = np.array(y[4:7])
        if X_old is not None and np.max(np.abs(X - X_old)) / X.sum() < 1e-9:
            break
        X_old = X
    i_f, S_f = y[0], y[1] + y[2] + y[3]
    assert S_f < 1e-4 and i_f < 1e-5, f"未饱和: i_f={i_f}, S_f={S_f}"
    return X / X.sum(), i_f


def conj_a(mu, a1, a2):
    v = np.array([mu**2, mu * a1, a1 * a2], float)
    return v / v.sum()


def conj_b(mu, a1, a2):
    """嵌套伯努利链: x_m = [prod_{j<m} a_j/(mu+a_j)] * mu/(mu+a_m)。"""
    x1 = mu / (mu + a1)
    x2 = (a1 / (mu + a1)) * (mu / (mu + a2))
    x3 = (a1 / (mu + a1)) * (a2 / (mu + a2))
    return np.array([x1, x2, x3])


def main():
    rng = None
    out = {'meta': 'f=3 stage-chain ledger conjecture check (Moderate #7)',
           'cases': []}
    cases = [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0), (0.5, 1.0, 2.0),
             (2.0, 3.0, 0.5), (1.5, 2.0, 4.0), (0.75, 0.5, 1.5)]
    for mu, a1, a2 in cases:
        x, i_f = ledger(mu, a1, a2)
        ca, cb = conj_a(mu, a1, a2), conj_b(mu, a1, a2)
        dev_a = np.max(np.abs(x - ca))
        dev_b = np.max(np.abs(x - cb))
        out['cases'].append(dict(mu=mu, alpha1=a1, alpha2=a2,
                                 ledger=[float(v) for v in x],
                                 conjecture_A=[float(v) for v in ca],
                                 conjecture_B_nested=[float(v) for v in cb],
                                 max_abs_dev_A=float(dev_a),
                                 max_abs_dev_B=float(dev_b),
                                 i_f=float(i_f)))
        print(f"mu={mu} a1={a1} a2={a2}: ledger={np.round(x,5)} "
              f"devA={dev_a:.4f} devB={dev_b:.2e}")
    json.dump(out, open(f'{ROOT}/results/w8_f3_conjecture.json', 'w'), indent=1)
    print('-> results/w8_f3_conjecture.json')


if __name__ == '__main__':
    main()
