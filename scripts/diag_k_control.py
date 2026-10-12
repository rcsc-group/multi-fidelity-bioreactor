"""Why does the test K positive control (s = 1e5) fail at L9 for dtmix_0.95 only? (diagnostic, not pre-registered)

Prints chi(t) near the 0.95 crossing for k_l<L>_s1e5 at L7, L8, L9 (solver-consistent chi via
scripts/diag_chi_slope.chi_of) and the liquid tracer mass drift; writes experiments/multifidelity/diag_k_control.png.
Usage: uv run python scripts/diag_k_control.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts import postprocess as pp  # noqa: E402
from scripts.chi_scale import locate  # noqa: E402
from scripts.diag_chi_slope import chi_of  # noqa: E402

RUNS = {7: "k_l7_s1e5", 8: "k_l8_s1e5", 9: "k_l9_s1e5"}


def mass_drift(run_id):
    d = locate(run_id)["data_dir"]
    t, c_sum = pp._joined_cols(d, "tr_oxy.dat", [pp._COL_C2_LIQ_SUM])
    i0 = int(np.where(c_sum > 1e-10 * np.abs(c_sum).max())[0][0])
    return float(c_sum[-1] / c_sum[i0] - 1.0), float(np.max(np.abs(c_sum[i0:] / c_sum[i0] - 1.0)))


def main():
    fig, ax = plt.subplots(figsize=(7, 4))
    for L, r in RUNS.items():
        try:
            ts, chi, per = chi_of(r)
        except Exception as e:  # noqa: BLE001
            print(r, "unavailable:", e)
            continue
        res = json.loads((ROOT / "runs" / r / "results.json").read_text())
        dm_end, dm_max = mass_drift(r)
        print(f"{r}: period {per:.3f} s, dtmix_0.95 (results) {res.get('dtmix_0.95')}, chi max {chi.max():.4f}, "
              f"tracer mass drift end {dm_end:+.3%} max {dm_max:.3%}")
        for thr in (0.90, 0.93, 0.94, 0.95):
            i = int(np.argmax(chi >= thr))
            print(f"   first chi >= {thr}: " + (f"{ts[i]:.2f} s" if chi[i] >= thr else "never"))
        m = ts < 15
        print("   chi at t = 2,3,4,5,6,8,10,12 s: " + ", ".join(
            f"{np.interp(x, ts, chi):.4f}" for x in (2, 3, 4, 5, 6, 8, 10, 12)))
        ax.plot(ts[m], chi[m], ls=":", color=fs.LEVEL_COLORS[L] if hasattr(fs, "LEVEL_COLORS") else None,
                label=f"L{L}")
    ax.axhline(0.95, color="0.5", lw=0.6)
    ax.set_xlabel("t after release (s)")
    ax.set_ylabel(r"$\chi$")
    ax.set_title(r"$D \times 10^5$")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    out = ROOT / "experiments/multifidelity/diag_k_control.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
