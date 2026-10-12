"""Figure (test V): ratio of the tracer variance actually destroyed to what physical diffusion destroys on the
resolved field, D_eff/D_phys, during mixing (chi 0.5 -> 0.95), L8 and L9 at 32.5 rpm. A physically controlled decay
sits on the line D_eff/D_phys = 1. Data: experiments/multifidelity/test_v_variance_budget.json.
Usage: uv run python scripts/fig_variance_budget.py
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

RUNS = {8: "scal_l8_rpm32.5_th7", 9: "scal_l9_rpm32.5_th7"}


def main():
    data = json.loads((ROOT / "experiments/multifidelity/test_v_variance_budget.json").read_text())["runs"]
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for L, run in RUNS.items():
        chi, r = np.array(data[run]["chi"]), np.array(data[run]["r"], dtype=float)
        m = (chi >= 0.5) & (chi <= 0.95) & np.isfinite(r)
        ax.plot(chi[m], r[m], ls=":", lw=1.2, color=fs.level_colour(L), label=f"L{L}")
    ax.axhline(1.0, color=fs.NEUTRAL, ls="--", lw=1.0, label="physical diffusion")
    ax.set_yscale("log")
    ax.set_ylim(0.5, 3e3)
    ax.set_xlim(0.5, 0.95)
    ax.set_xlabel(r"$\chi$")
    ax.set_ylabel(r"$D_{\mathrm{eff}}/D_{\mathrm{phys}}$")
    ax.grid(**fs.GRID_KW)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = ROOT / f"experiments/multifidelity/fig_variance_budget.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print("wrote", out)


if __name__ == "__main__":
    main()
