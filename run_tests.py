"""一键回归：聚合物三锚点 + 谣言 Gate-1（符号 / 数值 / SSA 三段）。"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
stages = [
    ("聚合物回归锚点", ["tests/test_polymer_anchors.py"], None),
    ("谣言 Gate-1 符号段", ["tests/test_rumor_gate1.py"], "sym"),
    ("谣言 Gate-1 数值段", ["tests/test_rumor_gate1.py"], "num"),
    ("谣言 Gate-1 SSA 段", ["tests/test_rumor_gate1.py"], "ssa"),
]

if __name__ == "__main__":
    failures = []
    for name, args, stage in stages:
        env = {"PYTHONPATH": str(HERE), "PATH": "/usr/bin:/bin:/usr/local/bin",
               "HOME": str(Path.home()), "RUMOR_STAGE": stage or "all"}
        print(f"\n===== {name} =====", flush=True)
        r = subprocess.run([sys.executable, *args], cwd=HERE, env=env)
        if r.returncode != 0:
            failures.append(name)
    print("\n" + "=" * 46)
    if failures:
        print("FAILED:", ", ".join(failures))
        sys.exit(1)
    print("全部回归通过 —— 闸门 1 正式关闭")
