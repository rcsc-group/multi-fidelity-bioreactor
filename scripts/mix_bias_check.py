"""Size of the mixing-time bias in the database (diary 2026-10-29). On our L8 DNS (20 rpm, theta 7), time to reach
chi = 0.10, 0.15, 0.20 measured three ways:
  old   : the production pipeline (postprocess.py on tr_oxy.dat: c2*f sums, sampled every t_out), seed f > 0.5
  sums  : corrected sums (bulk liquid, phi1), every step, same run (seed f > 0.5)
  seed  : corrected sums, Kim's seed f == 1 (pr28_Aseed)
Kim's DNS (corrected sums, his seed) is the reference.  Usage: uv run python scripts/mix_bias_check.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import prereg_offline as po  # noqa: E402

T = 0.607329
LEVELS = (0.10, 0.15, 0.20)
RUNS = Path("/oscar/scratch/eaguerov/replay_gate")


def from_diag(p: Path) -> tuple[np.ndarray, np.ndarray]:
    a = np.genfromtxt(p, names=True)
    a = a[int(np.argmax(a["C1"] > 0)):]
    mean = a["B1"] / a["B0"]
    var = a["B2"] / a["B0"] - mean ** 2
    return (a["t"] - a["t"][0]) / T, 1 - var / var[0]


def t_reach(t: np.ndarray, chi: np.ndarray, lev: float) -> float:
    k = int(np.argmax(chi >= lev))
    return float(t[k]) if chi[k] >= lev else np.nan


def main() -> None:
    t_old, _, chi_old, _ = po.series("pr28_A")           # production pipeline on the same run
    cases = {"old": (t_old, chi_old),
             "sums": from_diag(RUNS / "runs/pr28_A/prereg_diag.dat"),
             "seed": from_diag(RUNS / "mpi_runs/pr28_Aseed/prereg_diag.dat"),
             "Kim": from_diag(Path("/oscar/scratch/eaguerov/kim_pipeline/prereg_out/pre_A_N256/off_m1/prereg_diag.dat"))}
    print("time (cycles) to reach chi =", LEVELS)
    for name, (t, chi) in cases.items():
        print(f"{name:5s}", "  ".join(f"{t_reach(t, chi, lev):6.3f}" for lev in LEVELS))


if __name__ == "__main__":
    main()
