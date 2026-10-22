"""Zero-compute check of the premise behind the repeated-period low fidelity (meeting 2026-10-22, card 04):
is the flow of the LAST whole cycles periodic? Per run: relative cycle-to-cycle spread (sd / mean, 3 cycles of the QoI window) of the
EDR and the signed spatial-mean shear (tau_amp), and the largest cycle-to-cycle drift of the EDR (max |c_k - c_k-1| / mean).
Levels L8 (kmix_l8 / t6_l8), L9 (fig9_l9 / x8_l9), L10 (hf runs). Output: median and max over runs per level.
Usage: uv run python scripts/test_periodicity_levels.py
"""
import sys
import numpy as np

sys.path[:0] = [".", "scripts"]
from scripts.fig_h3_rpm import cycle_sd  # noqa: E402
from scripts.test_h1_mfbml import hf_points, hf_run, lf_run, has_results, RPMS, LF_THETAS  # noqa: E402
from scripts.test_h2_mfbml import l9_points, l9_run  # noqa: E402

LEVELS = {"L8": [lf_run(r, t) for r in RPMS for t in LF_THETAS],
          "L9": [l9_run(r, t) for r, t in l9_points()],
          "L10": [hf_run(r, t) for r, t in hf_points()]}
for lv, runs in LEVELS.items():
    out = {"edr_max_t": [], "tau_amp_t": [], "edr_mean_t": []}
    for run in runs:
        if not has_results(run):
            continue
        try:
            for q in out:
                out[q].append(cycle_sd(run, q))
        except Exception:
            pass
    n = len(out["edr_mean_t"])
    print(lv, "runs", n, " | ".join(f"{q}: median {100*np.median(v):.1f}% max {100*np.max(v):.1f}%" for q, v in out.items() if v))
