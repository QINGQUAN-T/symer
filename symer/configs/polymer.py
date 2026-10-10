"""聚合物家族配置：AA+BB 逐步缩聚的两步动力学（ChemRxiv 原文模型）。

状态（以 [BB]_0 归一）：
    b0 : 未反应 BB；b1 : 半反应 B*（冻结库存载体）；b2 : 全反应。
燃料池：sigma = [A]/[BB]_0 = 2(r-1) + 2 b0 + b1（消耗型载体）。
通道：
    P1: BB + A -> B*      速率 2 k sigma b0   （对称因子 2：两个等价 B 基团）
    P2: B* + A -> b2      速率 alpha k sigma b1（k2 = alpha k）
守恒律：pB = 1 - b0 - b1/2（B 侧转化率），pA = pB / r，sigma = 2(r - pB)。
锚点：Carothers Xn=(1+r)/(1+r-2pB)；alpha=1 时 1-pB = sqrt(b0)；
     轨迹式 (6)：b1 = 2(b0^{alpha/2} - b0)/(2 - alpha)。
"""
import sympy as sp

from ..spec import Channel, ReactionSystem


def build_polymer(alpha=None, r=None, k=None):
    """构建 AA+BB 两步缩聚符号系统。

    参数符号缺省自动创建；alpha=2 为共振分支（对数尾），由测试单独验证。
    """
    alpha = alpha if alpha is not None else sp.Symbol("alpha", positive=True)
    r = r if r is not None else sp.Symbol("r", positive=True)
    k = k if k is not None else sp.Symbol("k", positive=True)

    b0, b1, b2 = (sp.Symbol(n, nonnegative=True) for n in ("b0", "b1", "b2"))
    sigma = sp.Symbol("sigma", nonnegative=True)  # [A]/[BB]_0，消耗型燃料

    channels = [
        Channel(
            name="P1_nucleation",
            stoich={b0: -1, b1: +1, sigma: -1},
            rate=2 * k * sigma * b0,
            note="BB + A -> B*，对称因子 2（两个等价 B 基团）",
        ),
        Channel(
            name="P2_growth",
            stoich={b1: -1, b2: +1, sigma: -1},
            rate=alpha * k * sigma * b1,
            note="B* + A -> 链增长，k2 = alpha k",
        ),
    ]

    pB = 1 - b0 - b1 / 2  # B 侧转化率（守恒律，配置化声明）
    system = ReactionSystem(
        name="AA+BB two-step polycondensation",
        states={"b0": b0, "b1": b1, "b2": b2, "sigma": sigma},
        params={"alpha": alpha, "r": r, "k": k},
        channels=channels,
        conserved=[
            (sigma - (2 * (r - 1) + 2 * b0 + b1), 0),   # 燃料账本 sigma = 2(r-1)+2b0+b1
            (pB - (b1 / 2 + b2), 0),                     # B 基团守恒 S7.3
        ],
        carrier=sigma,  # 消耗型载体：轨迹方程中必须完全消去
    )
    return system


def freezing_line(b0, b1, r):
    """冻结线（终端条件配置化）：sigma = 0  <=>  2 b0 + b1 = 2(1 - r)。"""
    return sp.Eq(2 * b0 + b1, 2 * (1 - r))


def xn_all(pB, r):
    """Carothers 数均：Xn = (1+r)/(1+r-2pB)。"""
    return (1 + r) / (1 + r - 2 * pB)


def xw_all(pB, b1, r):
    """SI (S7.7) 双向游走重均：Xw = 1 + 2pB(3+pA-qB)/((1+r)w)。

    pA = pB/r, qB = b1/(2pB), w = 1 - rho, rho = (pB - b1/2)/r。
    """
    pA = pB / r
    qB = b1 / (2 * pB)
    rho = (pB - b1 / 2) / r
    w = 1 - rho
    return 1 + 2 * pB * (3 + pA - qB) / ((1 + r) * w)
