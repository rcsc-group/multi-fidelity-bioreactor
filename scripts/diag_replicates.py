"""Run-to-run spread (aleatoric noise s(x, h)) from release-shift replicates.

For each level and rpm: runs released after 80, 82 and 85 cycles (same binary
and protocol; L6/L8 release-80 runs are the kmix ladder, L7 release-80 is the
cold rep_l7_*_rel80 because the L7 ladder was warm-started). Reports the
coefficient of variation (sd/mean, n = 3) per QoI, and plots it against level.
Tests assumption A6 (the noise may change with h, in either direction).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs   # noqa: E402

RPMS = (17.5, 25, 32.5)
QOIS = [("dtmix_0.50", r"$\Delta t_{0.50}$"), ("dtmix_0.95", r"$\Delta t_{0.95}$"),
        ("kLa_1T_25", r"$k_La$ ($C^*=25\%$)")]


def runs(level, rpm):
    base = {6: f"kmix_l6_rpm{rpm:g}", 7: f"rep_l7_rpm{rpm:g}_rel80", 8: f"kmix_l8_rpm{rpm:g}",
            9: f"fig9_l9_rpm{rpm:g}"}[level]
    return [base, f"rep_l{level}_rpm{rpm:g}_rel82", f"rep_l{level}_rpm{rpm:g}_rel85"]


def cv(level, rpm, key):
    vals = []
    for r in runs(level, rpm):
        f = ROOT / "runs" / r / "results.json"
        if not f.exists():
            return np.nan
        vals.append(json.load(open(f)).get(key, np.nan))
    v = np.array(vals, float)
    return v.std(ddof=1) / v.mean() if np.all(np.isfinite(v)) else np.nan


def main():
    levels = (6, 7, 8, 9)
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, len(QOIS), figsize=(12, 3.6), sharey=True)
    for ax, (key, lab) in zip(axes, QOIS):
        for rpm, m in zip(RPMS, ("o", "s", "^")):
            c = np.array([cv(l, rpm, key) for l in levels])
            print(f"{key:11s} {rpm:5g} rpm  CV% by level 6/7/8/9: {np.round(100 * c, 2).tolist()}")
            ok = np.isfinite(c) & (c > 0)
            ax.semilogy(np.array(levels)[ok], 100 * c[ok], m + "-", ms=5, lw=1, label=f"{rpm:g} rpm")
        ax.set_title(lab, loc="left", fontsize=10)
        ax.set_xlabel("grid level")
        ax.set_xticks(levels)
        ax.grid(**fs.GRID_KW)
    axes[0].set_ylabel("run-to-run CV (%), n = 3")
    axes[-1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    out = ROOT / "experiments/multifidelity/replicate_cv.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
