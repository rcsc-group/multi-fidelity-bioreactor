"""Spin-up test (assumption A11): does the release cycle change the QoIs?

L8, 32.5 rpm, theta 7 deg, cold start, tracer and oxygen released after
33, 50, 80 (ladder run kmix_l8_rpm32.5), 82 and 85 cycles. 82/85 vs 80 is the
replicate spread; 33/50 vs 80 is the effect of a short spin-up. The only L10
point (fig9_l10b_seg2) was released after 33 cycles.
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

RUNS = {33: "rep_l8_rpm32.5_rel33", 50: "rep_l8_rpm32.5_rel50", 80: "kmix_l8_rpm32.5",
        82: "rep_l8_rpm32.5_rel82", 85: "rep_l8_rpm32.5_rel85"}
KEYS = [("dtmix_0.50", r"$\Delta t_{0.50}$"), ("dtmix_0.75", r"$\Delta t_{0.75}$"),
        ("dtmix_0.95", r"$\Delta t_{0.95}$"), ("kLa_1T_25", r"$k_La$ ($C^*=25\%$)")]


def main():
    rel = np.array(sorted(RUNS))
    val = np.array([[json.load(open(ROOT / "runs" / RUNS[r] / "results.json"))[k] for k, _ in KEYS]
                    for r in rel])
    ref = val[rel >= 80].mean(0)
    rep_cv = val[rel >= 80].std(0, ddof=1) / ref
    print("release:", rel.tolist())
    for j, (k, _) in enumerate(KEYS):
        print(f"{k:11s} values {np.round(val[:, j], 2).tolist()}  rel. to mean(80,82,85): "
              f"{np.round(100 * (val[:, j] / ref[j] - 1), 1).tolist()} %  replicate CV {100 * rep_cv[j]:.2f}%")
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for j, (k, lab) in enumerate(KEYS):
        ax.plot(rel, 100 * (val[:, j] / ref[j] - 1), "o-", ms=5, lw=1.0, label=lab)
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlabel("cycles before tracer and oxygen release")
    ax.set_ylabel("change against release at 80-85 cycles (%)")
    ax.set_xticks(rel)
    ax.grid(**fs.GRID_KW)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    out = ROOT / "experiments/multifidelity/spinup_test_l8_32p5.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
