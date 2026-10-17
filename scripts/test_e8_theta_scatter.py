"""Test E8 (pre-registered, diary 2026-10-16): does the L8 scatter of E0 (theta 7) also hold at theta 4 and 9?

Per theta: 21 runs on the 0.625-rpm grid 20-32.5 (t6_l8 at the 2.5-rpm points, e8_l8 elsewhere). E0 analyse() quadratic trend in log:
residual sd s and lag-1 autocorrelation rho1. The +-0.1 rpm pairs of E0 do not exist at theta 4/9, so the NOISE-LIKE reading uses rho1
only (deviation from the registered rule, disclosed).
SAME-SCATTER if s(theta 4) and s(theta 9) are each within a factor 1.5 of s(theta 7) for both primaries (kLa_1T_25 0.189, dtmix_0.95 0.141).
Usage: uv run python scripts/test_e8_theta_scatter.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.test_e0_scatter import RHO_BAR, analyse  # noqa: E402

GRID = [20 + 0.625 * i for i in range(21)]
S7 = {"kLa_1T_25": 0.189, "dtmix_0.95": 0.141, "kLa_1T_10": 0.095, "dtmix_0.50": 0.104}


def run_id(r, th):
    if min(abs(r - g) for g in (20, 22.5, 25, 27.5, 30, 32.5)) < 1e-9:
        return f"t6_l8_rpm{r:g}_th{th:g}"
    return f"e8_l8_rpm{r:g}_th{th:g}"


def val(run, q):
    f = ROOT / "runs" / run / "results.json"
    return float(json.loads(f.read_text()).get(q, np.nan)) if f.exists() else np.nan


def main():
    out, same = {}, {}
    for th in (4.0, 9.0):
        for q, s7 in S7.items():
            y = np.log([val(run_id(r, th), q) for r in GRID])
            ok = np.isfinite(y)
            if ok.sum() < 15:
                print(f"theta {th:g} {q}: only {ok.sum()} finite runs, skipped")
                continue
            res = analyse(np.array(GRID)[ok], y[ok], [np.nan], 2)
            reading = "NOISE-LIKE (rho1 only)" if abs(res["rho1"]) <= RHO_BAR else "STRUCTURED" if res["rho1"] > RHO_BAR else "ANTI-CORRELATED"
            ratio = res["s"] / s7
            out[f"{th:g}/{q}"] = dict(s=res["s"], rho1=res["rho1"], n=int(ok.sum()), ratio_to_theta7=ratio, reading=reading)
            same[(th, q)] = 1 / 1.5 <= ratio <= 1.5
            tag = "PRIMARY" if q in ("kLa_1T_25", "dtmix_0.95") else "reported"
            print(f"theta {th:g} {q:11s} [{tag}] n {ok.sum()}  s {res['s']:.3f} (theta 7: {s7:.3f}, ratio {ratio:.2f})  rho1 {res['rho1']:+.2f}  {reading}")
    prim = [same[(th, q)] for th in (4.0, 9.0) for q in ("kLa_1T_25", "dtmix_0.95") if (th, q) in same]
    verdict = "SAME-SCATTER" if prim and all(prim) else "THETA-DEPENDENT"
    out["verdict"] = verdict
    print("E8 verdict:", verdict, f"({sum(prim)}/{len(prim)} primary ratios within 1.5x)")
    (ROOT / "experiments/multifidelity/test_e8_theta_scatter.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
