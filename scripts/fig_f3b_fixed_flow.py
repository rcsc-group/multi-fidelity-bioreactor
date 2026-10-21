"""Fixed-flow test: the SAME recorded L8 flow (32.5 rpm, 7 deg, 22 cycles after release), oxygen and tracer solved on three meshes (L8, L9, L10).
Values from the replay runs (diary 2026-10-19, F3b): kLa_1T_10 46.93 / 38.31 / 50.21 1/h, kLa_1T_25 28.64 / 38.66 / 51.56 1/h, dtmix_0.50 11.79 / 19.03 / 24.74 s.
Plotted relative to the L8 mesh. Usage: uv run python scripts/fig_f3b_fixed_flow.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

DATA = {"kLa_1T_10": (fs.MK["cstar_10"], r"$k_La$, $C^*\!=\!0.10$", [46.93, 38.31, 50.21]),
        "kLa_1T_25": (fs.MK["cstar_25"], r"$k_La$, $C^*\!=\!0.25$", [28.64, 38.66, 51.56]),
        "dtmix_0.50": (fs.MK["chi_0.50"], r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$", [11.79, 19.03, 24.74])}


def main():
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for k, (mk, lab, v) in DATA.items():
        ax.plot([8, 9, 10], [100 * x / v[0] for x in v], ls=":", lw=1.0, color=fs.level_colour(8), marker=mk, mfc="w", mec=fs.level_colour(8), ms=7, mew=1.2, label=lab)
    ax.axhline(100, color=fs.GUIDE, lw=1.0)
    ax.set_xticks([8, 9, 10])
    ax.set_xticklabels(["L8", "L9", "L10"])
    ax.set_xlabel("mesh of the oxygen and tracer equations")
    ax.set_ylabel("value relative to the L8 mesh (%)")
    ax.grid(**fs.GRID_KW, axis="y")
    lg = ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / f"experiments/multifidelity/fig_f3b_fixed_flow.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=[lg])


if __name__ == "__main__":
    main()
