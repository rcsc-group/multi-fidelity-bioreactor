"""Figure: hydrodynamics-informed laws. (A) kLa(C* = 0.10) against EDR at L9 and L10; (B) 50 % mixing time against EDR at L8 and L9; (C) leave-one-out error for kLa at L10 (13 runs) by method.
Data: scripts/physics_law_kla_edr.py (table), physics_law_more.py. Usage: uv run python scripts/fig_physics_law.py
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
from scripts.physics_law_kla_edr import table  # noqa: E402


def dtmix_table(level):
    out = []
    pats = ("t6_l8_*", "t5_l8_*", "e8_l8_*") if level == 8 else ("fig9_l9_rpm*",)
    from scripts.hydro_dataset_v2 import qois
    for pat in pats:
        for d in sorted((ROOT / "runs").glob(pat)):
            f = d / "results.json"
            if "_ext" in d.name or not f.exists():
                continue
            v = json.loads(f.read_text()).get("dtmix_0.50")
            try:
                h = qois(d.name)
            except Exception:
                continue
            if v and v > 0 and h:
                out.append((h["edr_mean_t"], v))
    return np.array(out)


def main():
    T = table()
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8))
    ax = axes[0]
    for L in (9, 10):
        a = T[L]
        mk = fs.MK["cstar_10"]
        kw = fs.series_kw(fs.level_colour(L), mk, stat="mean", ours=True)
        kw.pop("ls")
        ax.plot(a[:, 3], a[:, 2], ls="none", **kw)
        s, c = np.polyfit(np.log(a[:, 3]), np.log(a[:, 2]), 1)
        xx = np.geomspace(a[:, 3].min(), a[:, 3].max(), 20)
        ax.plot(xx, np.exp(c) * xx ** s, color=fs.level_colour(L), ls=":", lw=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("dissipation rate EDR (W/m$^3$)")
    ax.set_ylabel(r"$k_La$, $C^*\!=\!0.10$ (1/h)")
    ax.grid(**fs.GRID_KW)
    ax = axes[1]
    for L in (8, 9):
        a = dtmix_table(L)
        kw = fs.series_kw(fs.level_colour(L), fs.MK["chi_0.50"], stat="mean", ours=True)
        kw.pop("ls")
        ax.plot(a[:, 0], a[:, 1], ls="none", **kw)
        s, c = np.polyfit(np.log(a[:, 0]), np.log(a[:, 1]), 1)
        xx = np.geomspace(a[:, 0].min(), a[:, 0].max(), 20)
        ax.plot(xx, np.exp(c) * xx ** s, color=fs.level_colour(L), ls=":", lw=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("dissipation rate EDR (W/m$^3$)")
    ax.set_ylabel(r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$ (s)")
    ax.grid(**fs.GRID_KW)
    ax = axes[2]
    names = ["constant", "rpm and angle", "L10 runs only", "chain from L9", "L9 times ratio", "EDR law"]
    vals = [0.569, 0.207, 0.225, 0.188, 0.167, 0.158]
    colors = ["0.6", "0.6", fs.level_colour(10), fs.level_colour(9), fs.level_colour(9), fs.level_colour(10)]
    ax.barh(range(len(names)), vals, color=colors, alpha=0.85, height=0.6)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("error on unseen L10 runs (RMS of ln ratio)")
    ax.grid(**fs.GRID_KW, axis="x")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / f"experiments/multifidelity/fig_physics_law.{ext}", dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
