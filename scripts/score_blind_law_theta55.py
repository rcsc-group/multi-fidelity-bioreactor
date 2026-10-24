"""Score the pre-registered blind test of the kLa-EDR law (diary 2026-10-21): frozen predictions in experiments/multifidelity/physics_law_frozen_theta55.json
against the measured kLa_1T_10 of the five theta-5.5 L10 runs (t10_l10_rpm{r}_th5.5).
PASS: law RMS ln error <= 0.20 and <= 0.8 x that of the (rpm, theta) regression; MARGINAL: <= 0.20 only; else FAIL.
Also scores the W9 replicate test (diary 2026-10-23): sd of ln dtmix_0.50 over the 3 warm runs and the gap to the cold run.
Usage: uv run python scripts/score_blind_law_theta55.py"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
fz = json.loads((ROOT / "experiments/multifidelity/physics_law_frozen_theta55.json").read_text())["predictions"]
el, er = [], []
for r, p in fz.items():
    d = json.loads((ROOT / f"runs/t10_l10_rpm{r}_th5.5/results.json").read_text())
    k = d["kLa_1T_10"]
    el.append(np.log(p["law"] / k)); er.append(np.log(p["rpm_theta"] / k))
    print(f"{r:>5} rpm: measured {k:.2f}  law {p['law']:.2f} ({100*(p['law']/k-1):+.0f} %)  rpm+theta {p['rpm_theta']:.2f} ({100*(p['rpm_theta']/k-1):+.0f} %)")
rl, rr = float(np.sqrt(np.mean(np.square(el)))), float(np.sqrt(np.mean(np.square(er))))
verdict = "PASS" if rl <= 0.20 and rl <= 0.8 * rr else "MARGINAL" if rl <= 0.20 else "FAIL"
print(f"RMS ln error: law {rl:.3f}, rpm+theta {rr:.3f}, ratio {rl/rr:.2f} -> {verdict}")

ld = lambda n: np.log(json.loads((ROOT / f"runs/{n}/results.json").read_text())["dtmix_0.50"])  # noqa: E731
warm = np.array([ld("w9n8_l9_rpm32.5_th7"), ld("w9_l9_rpm32.5_th7"), ld("w9n18_l9_rpm32.5_th7")]); cold = ld("fig9_l9_rpm32.5")
sd = float(warm.std(ddof=1)); gap = float(cold - warm.mean())
v = "SYSTEMATIC offset" if sd <= 0.15 and abs(gap) >= 0.3 else "INTRINSIC sensitivity" if sd >= 0.30 else "MIXED"
print(f"W9 replicate: ln dtmix_0.50 warm {np.round(warm, 3).tolist()} sd {sd:.3f}; cold {cold:.3f}; cold - warm mean {gap:+.3f} -> {v}")
