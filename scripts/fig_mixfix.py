"""Mixing-time bias check (diary 2026-10-29): free surface and dtmix, OLD vs NEW production, 3 conditions at L8.
Usage: uv run python scripts/fig_mixfix.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts import postprocess as pp  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts.mixfix_check import CONDS, OUT, rid  # noqa: E402

TH = ("0.50", "0.75", "0.95")


def main() -> None:
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    cols = ["tab:orange", "tab:blue", "tab:green"]
    for (rpm, th), c in zip(CONDS, cols):
        d = {}
        for kind, ls, lw in (("old", "-", 3.0), ("new", "--", 1.3)):
            run = rg.RUNS / rid(kind, rpm, th)
            p = json.loads((run / "params.json").read_text())
            v = np.loadtxt(run / "vol_frac_interf.dat", skiprows=1)
            dmax = v[:, 4]
            print(f"{rpm} {th} {kind}: posY_max mean {dmax[len(dmax)//2:].mean():.6f} max {dmax.max():.6f}")
            d[kind] = pp._compute_mixing_metrics(run, p, use_bulk=(kind == "new"))
        r = [d["new"][f"dtmix_{k}"] / d["old"][f"dtmix_{k}"] for k in TH]
        ax[0].plot([float(k) for k in TH], r, "o-", color=c, label=f"{rpm:g} rpm, {th:g}°")
        for k, m in zip(TH, ("o", "s", "^")):
            ax[1].plot(d["old"][f"dtmix_{k}"], d["new"][f"dtmix_{k}"], m, color=c, ms=8)
    lim = [0, 115]
    ax[1].plot(lim, lim, color="0.6", lw=0.8)
    ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].set_xlabel("dtmix, old (s)"); ax[1].set_ylabel("dtmix, new (s)")
    ax[0].set_xlabel(r"mixing threshold $\chi$"); ax[0].set_ylabel("dtmix new / old")
    for a in ax:
        a.grid(**fs.GRID_KW)
    fig.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"mixfix.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "mixfix.png")


if __name__ == "__main__":
    main()
