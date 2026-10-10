"""符号反应系统声明层：状态、通道、守恒律的配置化。

约定
----
* ``states``: 状态符号表（浓度或归一化份数）。
* ``params``: 无量纲参数符号表（如聚合的 alpha, r；谣言的 mu, alpha）。
* ``channels``: 反应通道。每个通道给出化学计量变化 ``stoich`` 与总速率
  表达式 ``rate``（速率表达式须已含伙伴池因子与对称因子）。
* ``conserved``: 守恒律列表 ``(expr, value)``，表示 expr == value；
  用于消元与冻结线（终端条件）的配置化。
* ``carrier``: 预期在轨迹方程中完全消去的公共因子——聚合家族是燃料池
  [A]（消耗型），谣言家族是传播者池 S（催化型）。这是"催化家族接口"
  的核心：两类家族的区分只体现在 carrier 是否被通道消耗，引擎本身不变。
"""
from dataclasses import dataclass, field

import sympy as sp


@dataclass
class Channel:
    """一条反应通道。

    stoich : 状态符号 -> 化学计量变化量（整数）。
    rate   : 该通道的总速率表达式（平均场意义下，含所有密度因子）。
    note   : 人类可读说明。
    """

    name: str
    stoich: dict
    rate: sp.Expr
    note: str = ""


@dataclass
class ReactionSystem:
    name: str
    states: dict            # 名称 -> sympy 符号
    params: dict            # 名称 -> sympy 符号
    channels: list          # Channel 列表
    conserved: list         # [(表达式, 值), ...]
    carrier: sp.Expr        # 预期消去的公共载体因子
    time: sp.Symbol = field(default_factory=lambda: sp.Symbol("t", positive=True))

    # ------------------------------------------------------------------
    def rhs(self):
        """返回 {状态符号: d(state)/dt 的符号表达式}。"""
        out = {s: sp.Integer(0) for s in self.states.values()}
        for ch in self.channels:
            for s, d in ch.stoich.items():
                out[s] = out[s] + d * ch.rate
        return {s: sp.expand(e) for s, e in out.items()}

    def substitute(self, expr, subs_map):
        return sp.simplify(expr.subs(subs_map))

    # ------------------------------------------------------------------
    def describe(self):
        lines = [f"ReactionSystem<{self.name}>"]
        for ch in self.channels:
            stoich = " ".join(
                f"{'+' if d > 0 else ''}{d}{s}" for s, d in ch.stoich.items()
            )
            lines.append(f"  [{ch.name}] {stoich}   rate = {ch.rate}")
        for expr, val in self.conserved:
            lines.append(f"  conserved: {sp.simplify(expr - val)} == 0")
        return "\n".join(lines)
