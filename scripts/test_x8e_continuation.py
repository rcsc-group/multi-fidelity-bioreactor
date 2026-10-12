"""Verdict of the pre-registered x8 continuation (diary 2026-10-11): TRANSIENT / BIAS / UNSETTLED.

Runs x8e_l10_rpm<r>_th7 continue x8_l10_rpm<r>_th7 at L10 for 6 more cycles (cycles 7-12 after the L8 -> L10 restart).
Per cycle: tau_mean, tau_95, EDR (scripts/test_warmstart_x8.cyc_vals, same units and columns as the warm-start test).
Rule (tau_mean): TRANSIENT if, at both rpm, |mean(last 3 cycles)/xlevel - 1| <= 5% AND every |cycle-to-cycle change| in
the last 3 cycles < 0.5%. BIAS if at 17.5 rpm mean(last 3)/xlevel - 1 > 5% with all |changes| < 0.5%. Else UNSETTLED.
Usage: uv run python scripts/test_x8e_continuation.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.test_warmstart_x8 import QS, XL, cyc_vals, load, window_info  # noqa: E402

RPMS = ["17.5", "25"]
TOL, SETTLE = 0.05, 0.005


def per_rpm(rpm):
    xl = window_info(XL.format(rpm))
    x8 = load(f"x8_l10_rpm{rpm}_th7")
    xe = load(f"x8e_l10_rpm{rpm}_th7")
    if xl is None or x8 is None or xe is None or not (ROOT / f"runs/x8e_l10_rpm{rpm}_th7/results.json").exists():
        return None
    c0 = x8["p"]["t_checkpoint"] / x8["T"]
    c1 = xe["p"]["t_checkpoint"] / xe["T"]
    pre = cyc_vals(x8, c0, 6)
    post = [c for c in cyc_vals(xe, c1, 7) if c["n_steps"] >= 25]  # whole cycles only
    seq = pre + post
    ref = xl["q"]
    last3 = post[-3:]
    tm = np.array([c["tau_mean_t"] for c in last3])
    dev = {q: float(np.mean([c[q] for c in last3]) / ref[q] - 1) for q in QS}
    changes = [float(b / a - 1) for a, b in zip(tm[:-1], tm[1:])]
    changes = [float(tm[0] / post[-4]["tau_mean_t"] - 1)] + changes if len(post) >= 4 else changes
    settled = all(abs(c) < SETTLE for c in changes)
    return dict(xlevel=ref, n_post_cycles=len(post), dev_last3=dev, last_changes_tau_mean=changes, settled=settled,
                tau_mean_by_cycle=[c["tau_mean_t"] for c in seq], edr_by_cycle=[c["edr_mean_t"] for c in seq],
                tau95_by_cycle=[c["tau_95_t"] for c in seq])


def main():
    res = {r: per_rpm(r) for r in RPMS}
    for r, e in res.items():
        if e is None:
            print(f"{r} rpm: pending")
            continue
        print(f"{r} rpm: {e['n_post_cycles']} continuation cycles; xlevel tau_mean {e['xlevel']['tau_mean_t']:.5g}")
        print("   tau_mean by cycle 1..: " + ", ".join(f"{v:.4g}" for v in e["tau_mean_by_cycle"]))
        print("   EDR by cycle 1..:      " + ", ".join(f"{v:.4g}" for v in e["edr_by_cycle"]))
        print("   dev of last 3 vs xlevel: " + ", ".join(f"{q} {v:+.2%}" for q, v in e["dev_last3"].items()))
        print("   last cycle-to-cycle changes in tau_mean: " + ", ".join(f"{c:+.2%}" for c in e["last_changes_tau_mean"])
              + f"  settled={e['settled']}")
    if any(e is None for e in res.values()):
        verdict = "PENDING"
    elif all(abs(e["dev_last3"]["tau_mean_t"]) <= TOL and e["settled"] for e in res.values()):
        verdict = "TRANSIENT"
    elif res["17.5"]["dev_last3"]["tau_mean_t"] > TOL and res["17.5"]["settled"]:
        verdict = "BIAS"
    else:
        verdict = "UNSETTLED"
    print("VERDICT:", verdict)
    out = ROOT / "experiments/multifidelity/test_x8e_continuation.json"
    out.write_text(json.dumps(dict(verdict=verdict, per_rpm=res), indent=1))
    print("wrote", out)


if __name__ == "__main__":
    main()
