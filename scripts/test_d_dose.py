"""Test D (pre-registered, diary 2026-10-16): dose-response of the rpm jump. Base rpm 22.5 and 30, eps in 1e-4 ... 3e-2 (plus +-0.1 from E0 pairs).

Delta(eps) = |log q(base + eps) - log q(base)| for dtmix_0.50, dtmix_0.95 (primary), kLa_1T_10 (reported).
CHAOS-LIKE: Delta(1e-3) >= 0.5 x median Delta(|eps| = 0.1) at both bases. SMOOTH: log-log slope 0.8-1.2 over 1e-4..0.1 and Delta(0.1) <= 0.05.
THRESHOLD: Delta <= 0.02 for all eps <= e* and >= 0.08 for all eps >= 3 e* for some e* in the tested range, at that base. Else MIXED.
Usage: uv run python scripts/test_d_dose.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

EPS = (1e-4, 1e-3, 3e-3, 1e-2, 3e-2)
QS = ("dtmix_0.50", "dtmix_0.95", "kLa_1T_10")


def val(run, q):
    f = ROOT / "runs" / run / "results.json"
    return float(json.loads(f.read_text()).get(q, np.nan)) if f.exists() else np.nan


def main():
    ids = json.loads((ROOT / "experiments/multifidelity/test_t5_roughness.json").read_text())["run_ids"]
    out = {}
    for base in (22.5, 30.0):
        b = ids[f"{base:g}"]
        pair = {22.5: ("t5_l8_rpm22.4_th7", "t5_l8_rpm22.6_th7"), 30.0: ("t5_l8_rpm29.9_th7", "t5_l8_rpm30.1_th7")}[base]
        for q in QS:
            q0 = val(b, q)
            d = {e: abs(np.log(val(f"dose_l8_rpm{base + e:g}_th7", q)) - np.log(q0)) for e in EPS}
            d01 = [abs(np.log(val(p, q)) - np.log(q0)) for p in pair]
            out[f"{base:g}/{q}"] = dict(base_value=q0, delta={f"{e:g}": float(v) for e, v in d.items()}, delta_0p1=d01)
            row = "  ".join(f"{e:g}:{d[e]:.3f}" for e in EPS)
            print(f"base {base:g} {q:11s} base {q0:7.2f} | Delta(eps) {row} | +-0.1: {d01[0]:.3f}, {d01[1]:.3f}")
    verdicts = {}
    for q in QS[:2]:
        chaos, thr = [], []
        for base in (22.5, 30.0):
            r = out[f"{base:g}/{q}"]
            ref = float(np.median(r["delta_0p1"]))
            dl = [r["delta"][f"{e:g}"] for e in EPS]
            chaos.append(np.isfinite(ref) and r["delta"]["0.001"] >= 0.5 * ref)
            ok_thr = False
            for i, e in enumerate(EPS[:-1]):
                low = all(x <= 0.02 for x in dl[: i + 1])
                high = [x for ee, x in zip(EPS, dl) if ee >= 3 * e]
                if low and high and all(x >= 0.08 for x in high):
                    ok_thr = True
            thr.append(ok_thr)
        smooth = all(out[f"{b:g}/{q}"]["delta"]["0.03"] <= 0.05 for b in (22.5, 30.0))
        verdicts[q] = "CHAOS-LIKE" if all(chaos) else "THRESHOLD" if all(thr) else "SMOOTH?" if smooth else "MIXED"
        print(f"{q}: chaos per base {chaos}, threshold per base {thr} -> {verdicts[q]} (smooth screen: Delta(0.03) <= 0.05 at both: {smooth})")
    out["verdict"] = verdicts
    (ROOT / "experiments/multifidelity/test_d_dose.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
