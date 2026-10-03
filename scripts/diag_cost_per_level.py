"""Measured cost of the mixing runs: core-seconds per simulated second, L6-L9.

Source: logstats.dat (rank-0 wall clock against simulation time) of the
kmix_l6/l7/l8 and fig9_l9 rpm sweeps. Cost rate = slope of wall clock against
physical time over the post-release part, times the rank count of the level.
Also the cost to reach chi = 0.95 at each level: rate x (spin-up + dtmix_0.95).
Answers: the cost prior c ~ 2^(gamma*level); does gamma depend on rpm?
"""
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs                       # noqa: E402
from scripts.postprocess import _t_scales                # noqa: E402
from scripts.pilot_truncated_chi import RPMS             # noqa: E402

FAMILIES = {6: ("kmix_l6_rpm{:g}", 4), 7: ("kmix_l7_rpm{:g}", 8),
            8: ("kmix_l8_rpm{:g}", 16), 9: ("fig9_l9_rpm{:g}", 32)}
PAT = re.compile(r"t: ([0-9.eE+-]+) .*Wall clock time \(s\): ([0-9.eE+-]+)")


def rate(run, ranks):
    d = ROOT / "runs" / run
    tw = np.array([[float(a), float(b)] for a, b in PAT.findall((d / "logstats.dat").read_text())])
    T_bio, _ = _t_scales(json.load(open(d / "params.json")))
    t, w = tw[:, 0] * T_bio, tw[:, 1]
    keep = np.r_[True, np.diff(w) > 0]           # drop segment restarts
    dt, dw = np.diff(t[keep]), np.diff(w[keep])
    ok = (dt > 0) & (dw > 0)
    return ranks * dw[ok].sum() / dt[ok].sum()   # core-s per physical s


def main():
    rows = []
    for lv, (tmpl, ranks) in FAMILIES.items():
        for r in RPMS:
            run = tmpl.format(r)
            try:
                k = rate(run, ranks)
            except (FileNotFoundError, ValueError, IndexError):
                continue
            res = json.load(open(ROOT / "runs" / run / "results.json"))
            rows.append((lv, r, k, res.get("dtmix_0.95", np.nan)))
    a = np.array(rows, dtype=float)
    print(f"{'rpm':>5s}" + "".join(f"  L{lv} core-s/s" for lv in FAMILIES))
    for r in RPMS:
        print(f"{r:5g}" + "".join(
            f"{a[(a[:, 0] == lv) & (a[:, 1] == r), 2][0]:14.0f}" if ((a[:, 0] == lv) & (a[:, 1] == r)).any()
            else f"{'-':>14s}" for lv in FAMILIES))
    g = []
    for lv in (7, 8, 9):
        ratio = [a[(a[:, 0] == lv) & (a[:, 1] == r), 2][0] / a[(a[:, 0] == lv - 1) & (a[:, 1] == r), 2][0]
                 for r in RPMS if ((a[:, 0] == lv) & (a[:, 1] == r)).any() and ((a[:, 0] == lv - 1) & (a[:, 1] == r)).any()]
        g.append(np.log2(ratio))
        print(f"L{lv-1}->L{lv}: cost-rate ratio median {np.median(ratio):.1f} "
              f"(range {min(ratio):.1f}-{max(ratio):.1f}), gamma = {np.median(np.log2(ratio)):.2f}")

    plt.rcParams.update(fs.rcparams())
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
    cmap = plt.get_cmap("viridis")
    for k, r in enumerate(RPMS):
        m = a[:, 1] == r
        col = cmap(k / (len(RPMS) - 1))
        a1.semilogy(a[m, 0], a[m, 2], "o-", color=col, ms=3, lw=0.8, label=f"{r:g} rpm")
        a2.semilogy(a[m, 0], a[m, 2] * a[m, 3] / 3600, "o-", color=col, ms=3, lw=0.8)
    a1.set_xlabel("grid level"); a1.set_ylabel("cost rate (core-s per simulated s)")
    a2.set_xlabel("grid level"); a2.set_ylabel(r"core-h to observe $\Delta t_{0.95}$ after release")
    for ax in (a1, a2):
        ax.set_xticks([6, 7, 8, 9]); ax.grid(**fs.GRID_KW)
    a2.legend(*a1.get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    fig.tight_layout()
    out = ROOT / "experiments/multifidelity/cost_per_level.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
