"""Liquid velocity across a restart seam: does a continuation continue the flow?

    uv run python scripts/figures/fig_restart_seam.py --parent RUN --old RUN [RUN ...] --fixed RUN [RUN ...] --out PNG

Top panel: extensions written by the old rule (no *_prev keys: the solver ramps the forcing again from zero
after the restart). Bottom panel: extensions written by submit_continue.build_params (*_prev = the current
forcing). Every panel also shows the parent run before the seam. Line = the raw normf.dat series; hollow
markers = the mean of each period. Pitfall "restart-ramp" in docs_site/developer/solver-pitfalls.md.
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.figstyle import MK, level_colour, rcparams, series_kw, GUIDE  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402

BEFORE, AFTER = 4, 8   # periods shown before and after the seam


def velocity(rd: Path):
    lines = (rd / "normf.dat").read_text().splitlines()
    head = lines[0].split()
    d = np.array([[float(x) for x in ln.split()] for ln in lines[1:] if ln.strip()])
    return d[:, head.index("t")], np.hypot(d[:, head.index("ux_liq_rms")], d[:, head.index("uy_liq_rms")])


def draw(ax, rd: Path, t_seam: float, T: float, colour: str, alpha: float, label: str | None):
    t, v = velocity(rd)
    x = (t - t_seam) / T
    keep = (x >= -BEFORE) & (x <= AFTER)
    ax.plot(x[keep], v[keep], color=colour, lw=0.6, alpha=0.35 * alpha)
    k = np.floor(x[keep] + 1e-9)
    xs = sorted(set(k))
    means = [v[keep][k == i].mean() for i in xs]
    ax.plot(np.array(xs) + 0.5, means, alpha=alpha, label=label,
            **series_kw(colour, MK["velocity"], stat="mean", ours=True))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", type=Path, required=True)
    ap.add_argument("--old", type=Path, nargs="+", required=True)
    ap.add_argument("--fixed", type=Path, nargs="*", default=[])
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    p = json.loads((a.old[0] / "params.json").read_text())
    T, t_seam = _t_scales(p)[1], p["t_checkpoint"]
    colour = level_colour(p["fidelity"])
    plt.rcParams.update(rcparams())
    rows = [("old", a.old)] + ([("fixed", a.fixed)] if a.fixed else [])
    fig, axes = plt.subplots(len(rows), 1, figsize=(6.4, 2.4 * len(rows)), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, (rule, runs) in zip(axes, rows):
        draw(ax, a.parent, t_seam, T, colour, 0.35, "parent")
        for i, rd in enumerate(runs):
            draw(ax, rd, t_seam, T, colour, 1.0 if i == 0 else 0.6, rd.name)
        ax.axvline(0, color=GUIDE, lw=0.8)
        ax.set_ylabel("liquid velocity rms")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=8)
    axes[-1].set_xlabel("periods since the restart")
    fig.savefig(a.out, dpi=150, bbox_inches="tight")
    print("saved", a.out)


if __name__ == "__main__":
    main()
