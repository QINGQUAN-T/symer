"""谣言家族配置：两阶段传播者链式模型（W1 锁定版 v1.0-lock）。

状态（人口归一化，N=1）：
    i  : 无知者（单调递减的"燃料"侧，冻结条件 i = 0）
    s1 : 新传播者（被招募一次，未再接触）
    s2 : 被再激活的传播者（再接触后武装/预热）
    x  : 出传者 R（吸收态）
通道（k1=1 无量纲化，k2 = mu，再接触响应增强因子 alpha）：
    REC : S_j x I  -> S_j + S1   速率 k1 * i * (s1+s2)   招募（催化：传播者存活）
    ARM : S_j x S1 -> S_j + S2   速率 alpha * i 换成 alpha * s1 * (s1+s2)
         门控/武装：新传播者被任何传播者再接触，以 alpha*k1 升级（预热响应）
    STIF1/2: S_j x S_k -> S_j + R 速率 k2 * s_k * (s1+s2)
         出传：被再接触的新/熟传播者以 k2 出传（听到第二次而厌倦，DK 机制）
载体：S = s1 + s2（催化型传播者池，在轨迹方程中完全消去）。
无量纲群：mu = k2/k1（出传/招募比），gamma = alpha + mu（职业层引擎指数）。

W1 裁决（详见备忘录）：
  * 被接触驱动（被动）：招募与武装的速率都正比于伙伴池密度 -> 载体消去；
  * S2 吸收通道与 S1 相同（对称的"听到第二次"出传），保证饱和后全 R；
  * f>=3 阶段列入展望。
"""
import sympy as sp

from ..spec import Channel, ReactionSystem


def build_rumor(alpha=None, mu=None, k1=None):
    alpha = alpha if alpha is not None else sp.Symbol("alpha", nonnegative=True)
    mu = mu if mu is not None else sp.Symbol("mu", nonnegative=True)
    k1 = k1 if k1 is not None else sp.Integer(1)  # 无量纲化
    k2 = mu * k1

    i, s1, s2, x = (sp.Symbol(n, nonnegative=True) for n in ("i", "s1", "s2", "x"))
    S = s1 + s2  # 催化型传播者池（载体）

    channels = [
        Channel(
            name="REC",
            stoich={i: -1, s1: +1},
            rate=k1 * i * S,
            note="招募：任何传播者 x 无知者 -> 新传播者（消耗一个无知者）",
        ),
        Channel(
            name="ARM",
            stoich={s1: -1, s2: +1},
            rate=alpha * k1 * s1 * S,
            note="门控/武装：S1 被再接触 -> S2，速率 alpha*k1（预热响应）",
        ),
        Channel(
            name="STIF1",
            stoich={s1: -1, x: +1},
            rate=k2 * s1 * S,
            note="出传：新传播者被再接触 -> R（听到第二次而厌倦）",
        ),
        Channel(
            name="STIF2",
            stoich={s2: -1, x: +1},
            rate=k2 * s2 * S,
            note="出传：熟传播者被再接触 -> R",
        ),
    ]

    system = ReactionSystem(
        name="two-stage rumor cascade (v1.0-lock)",
        states={"i": i, "s1": s1, "s2": s2, "x": x},
        params={"alpha": alpha, "mu": mu, "k1": k1},
        channels=channels,
        conserved=[(i + s1 + s2 + x, 1)],  # 总人数守恒（归一化）
        carrier=S,  # 催化型载体：轨迹方程中必须完全消去
    )
    return system


def pooled_slope(mu):
    """合并层（alpha 盲）引擎：dS/di = -1 + mu * S / i。"""
    S, i = sp.Symbol("S", nonnegative=True), sp.Symbol("i", nonnegative=True)
    return -1 + mu * S / i


def career_slope(alpha, mu):
    """职业层引擎：ds1/di = -1 + gamma * s1 / i，gamma = alpha + mu。"""
    s1, i = sp.Symbol("s1", nonnegative=True), sp.Symbol("i", nonnegative=True)
    gamma = alpha + mu
    return -1 + gamma * s1 / i, gamma


def frozen_inventory_law(alpha, mu):
    """冻结库存律（W1 手推主结果）：(x1, x2) = (mu, alpha)/(alpha+mu)。

    x1 = 出传时仍为新传播者的账号数（转一次即沉默）
    x2 = 出传时已被再激活的账号数（多波转发后沉默）
    返回 alpha 反演式：alpha = mu * x2 / x1。
    """
    return sp.Eq(alpha, mu * sp.Symbol("x2") / sp.Symbol("x1"))
