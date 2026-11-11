"""Gate verdict + reasons for the 12 pr8 L8 runs, and the same for live bo4 chain L8 runs (comparison)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.bo2d import loop as LP  # noqa: E402

C = LP.C
for pre in ("pr8_l8_", "bo4_l8_"):
    n = f = 0
    for d in sorted((C.ROOT / "runs").glob(pre + "*")):
        if not (d / "params.json").exists() or "failed" in d.name:
            continue
        try:
            v, info = C.warm_run_gate(d)
        except Exception as e:  # noqa: BLE001
            v, info = "ERR", repr(e)[:80]
        n += 1; f += v == "FAILED-GATE"
        if pre == "pr8_l8_" or v == "FAILED-GATE":
            print(d.name, v, {k: info[k] for k in list(info)[:8]} if isinstance(info, dict) else info)
    print(pre, "total", n, "FAILED-GATE", f)
