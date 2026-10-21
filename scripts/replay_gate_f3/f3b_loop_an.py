"""Periodic-flow replay at L8: compare looped replays (K = 1, 2, 4 recorded periods, 22 cycles) with the coupled run (diary 2026-10-21).
PASS if kLa_1T_10, kLa_1T_25, dtmix_0.50 are all within 5 % of the coupled run for some K; MARGINAL within 10 %; else FAIL.
Usage: uv run python scripts/replay_gate_f3/f3b_loop_an.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import replay_gate as rg

Q = ("kLa_1T_10", "kLa_1T_25", "dtmix_0.50")
ref = rg.load_results(rg.run_dir("f3b_coupled"))
print("coupled", {k: round(ref[k], 3) for k in Q})
out = {}
best = None
for k in (1, 2, 4):
    try:
        r = rg.load_results(rg.run_dir(f"f3b_loop{k}"))
    except Exception as e:  # noqa: BLE001
        print(k, "missing", repr(e)[:80])
        continue
    d = {q: 100 * (r[q] / ref[q] - 1) for q in Q}
    out[k] = d
    worst = max(abs(v) for v in d.values())
    best = worst if best is None else min(best, worst)
    print(f"K = {k}: " + "  ".join(f"{q} {r[q]:.3f} ({d[q]:+.1f} %)" for q in Q))
verdict = "PASS" if best is not None and best <= 5 else "MARGINAL" if best is not None and best <= 10 else "FAIL"
print("verdict:", verdict, "(smallest worst-case deviation %.1f %%)" % best if best is not None else "")
Path(__file__).resolve().parents[2].joinpath("experiments/multifidelity/f3b_loop.json").write_text(json.dumps(dict(deviation_percent=out, verdict=verdict), indent=1))
