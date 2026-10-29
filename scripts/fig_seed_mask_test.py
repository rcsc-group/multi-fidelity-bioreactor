"""Two tests (diary 2026-10-29, pre-registered). Left: bulk chi of the DNS, ours with our seed rule (f > 0.5),
ours with Kim's (f == 1), Kim's. Right: tracer mass of the flux-from-u looped replay, ours without and with
Kim's solid-face zeroing, Kim's. All from prereg_diag.dat (same code both solvers).
Usage: uv run python scripts/fig_seed_mask_test.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

T = 0.607329
OURS = Path("/oscar/scratch/eaguerov/replay_gate/runs")
KIM = Path("/oscar/scratch/eaguerov/kim_pipeline/prereg_out/pre_A_N256")
OUT = ROOT / "docs/meetings/prereg"
C_OURS, C_KIM = fs.level_colour(8), fs.KIM

# (dir, label, colour, width, style): thick first
LEFT = [(KIM / "off_m1", "Kim", C_KIM, 3.5, "-"),
        (OURS / "pr28_A", "ours, seed f > 0.5", C_OURS, 2.0, "-"),
        (Path("/oscar/scratch/eaguerov/replay_gate/mpi_runs/pr28_Aseed"), "ours, seed f = 1", C_OURS, 1.2, "--")]
RIGHT = [(KIM / "his_loop", "Kim", C_KIM, 3.5, "-"),
         (OURS / "pr28_repWS", "ours", C_OURS, 2.0, "-"),
         (OURS / "pr28_repWSm", "ours, solid faces zeroed", C_OURS, 1.2, "--")]


def load(d: Path) -> dict[str, np.ndarray]:
    a = np.genfromtxt(d / "prereg_diag.dat", names=True)
    a = a[int(np.argmax(a["C1"] > 0)):]
    a = a[a["t"] <= a["t"][0] + 3 * T + 1e-9]
    mean = a["B1"] / a["B0"]
    var = a["B2"] / a["B0"] - mean ** 2
    return {"x": (a["t"] - a["t"][0]) / T, "chi": 1 - var / var[0], "m": 100 * (a["C1"] / a["C1"][0] - 1)}


def main() -> None:
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
    for k, (runs, key, lab) in enumerate(((LEFT, "chi", r"DNS mixing $\chi$"),
                                          (RIGHT, "m", "flux from u: tracer mass change (%)"))):
        for d, name, c, lw, ls in runs:
            r = load(d)
            ax[k].plot(r["x"], r[key], color=c, lw=lw, ls=ls, alpha=0.45 if lw > 3 else 1.0, label=name)
        ax[k].set_ylabel(lab)
        ax[k].set_xlabel("time from release (cycles)")
        ax[k].grid(**fs.GRID_KW)
        ax[k].legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.2), ncol=1)
    ax[0].set_xlim(0, 3)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"seed_mask_test.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "seed_mask_test.png")


if __name__ == "__main__":
    main()
