"""symer -- 符号核 v0：链式两步配对过程系统的符号引擎。

设计目标（对应项目计划书 W1）：
  1. 转移图 + 守恒律配置化：用 Channel/ReactionSystem 声明状态、通道、守恒律；
  2. 催化家族接口：partner_consumed 标志区分"消耗型伙伴"（聚合的 A 端基，
     每次反应被消耗）与"催化型伙伴"（谣言的传播者池，催化不消耗）——
     两个家族共享同一数学引擎；
  3. 轨迹方程推导：自动组装 ODE 右端、求 dY/dX、验证公共载体因子
     （聚合 [A] / 谣言传播者池 S）完全消去、符号求解闭式解、残差归零校验。
"""
from .spec import Channel, ReactionSystem
from .engine import (
    derive_trajectory,
    solve_trajectory,
    verify_solution,
    carrier_cancellation_check,
    TrajectoryResult,
)

__all__ = [
    "Channel",
    "ReactionSystem",
    "TrajectoryResult",
    "derive_trajectory",
    "solve_trajectory",
    "verify_solution",
    "carrier_cancellation_check",
]
