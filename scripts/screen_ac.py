"""Zero-compute screens for the mfbml fallback routes (pre-registered, diary 2026-10-13).

Screen A: does the L8 -> L10 correction of tau_mean / tau_95 / EDR at theta 7 track the L8 mean circulation better
than rpm? Screen C: do the test K D-scan curves of L8 and L9 collapse when plotted against D_eff (test V), with no
free parameter? Rules are in the diary; this script prints the numbers and the readings.
Usage: uv run python scripts/screen_ac.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.diag_x8_slow_mode import normf_cycles  # noqa: E402
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.test_warmstart_x8 import load  # noqa: E402

RPMS = [17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
QS = ("tau_mean_t", "tau_95_t", "edr_mean_t")


def screen_a():
    rows = []
    for r in RPMS:
        l8, l10 = f"kmix_l8_rpm{r:g}", f"l10_fig13a_xlevel_rpm{r:g}"
        q8, q10, s8 = qois(l8), qois(l10), load(l8)
        if not (q8 and q10 and s8):
            print(f"A: {r} missing")
            continue
        circ = normf_cycles(l8, s8, np.floor(s8["cyc"].max()) - 3, 3)
        rows.append(dict(rpm=r, circ=float(np.mean([c["Omega_liq_avg"] for c in circ])),
                         **{q: q10[q] / q8[q] for q in QS}))
    rpm = np.array([x["rpm"] for x in rows])
    circ = np.array([x["circ"] for x in rows])
    print("Screen A (theta 7): rpm, L8 circulation, L10/L8 ratios")
    for x in rows:
        print(f"  {x['rpm']:5.1f}  circ {x['circ']:6.3f}  " + "  ".join(f"{q} {x[q]:.3f}" for q in QS))
    out = {}
    for q in QS:
        y = np.array([x[q] for x in rows])
        rho = spearmanr(y, circ).correlation
        res_c = np.std(y - np.polyval(np.polyfit(circ, y, 1), circ))
        res_r = np.std(y - np.polyval(np.polyfit(rpm, y, 1), rpm))
        out[q] = dict(spearman_circ=float(rho), rms_res_circ=float(res_c), rms_res_rpm=float(res_r))
        print(f"  {q}: Spearman(ratio, circ) {rho:+.2f}; RMS residual linear-in-circ {res_c:.4f} vs "
              f"linear-in-rpm {res_r:.4f} (ratio {res_c / res_r:.2f})")
    a = out["tau_mean_t"]
    verdict = "WORTH BUILDING" if abs(a["spearman_circ"]) >= 0.8 and a["rms_res_circ"] <= 0.7 * a["rms_res_rpm"] \
        else "NOT SUPPORTED at theta 7"
    print("  Screen A reading:", verdict)
    return dict(rows=rows, stats=out, verdict=verdict)


def screen_c():
    runs = json.loads((ROOT / "experiments/multifidelity/test_k_summary.json").read_text())["runs"]
    r8, r9 = 669.0, 239.0
    s9, s8p = 1e3, 1e3 + r9 - r8
    out, big = {}, 0
    print(f"Screen C: L9 at s = {s9:g} vs L8 at s' = {s8p:.0f} (equal D_eff = {s9 + r9:.0f} D_phys)")
    for q in ("dtmix_0.50", "dtmix_0.75", "kLa_1T_25"):
        lo, hi = runs["L8"]["1e2"][q], runs["L8"]["1e3"][q]
        w = (np.log(s8p) - np.log(1e2)) / (np.log(1e3) - np.log(1e2))
        q8 = float(np.exp((1 - w) * np.log(lo) + w * np.log(hi)))
        q9 = runs["L9"]["1e3"][q]
        d = float(np.log(q9 / q8))
        big += abs(d) > 0.30
        out[q] = dict(L9=q9, L8_interp=q8, dlog=d)
        print(f"  {q}: L9 {q9:.3f} vs L8(interp) {q8:.3f}  dlog {d:+.2f}")
    if all(abs(v["dlog"]) <= 0.10 for v in out.values()):
        verdict = "COLLAPSE"
    elif big >= 2:
        verdict = "NO COLLAPSE"
    else:
        verdict = "INCONCLUSIVE"
    print("  Screen C reading:", verdict)
    return dict(points=out, verdict=verdict)


def main():
    res = dict(A=screen_a(), C=screen_c())
    (ROOT / "experiments/multifidelity/screen_ac.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
