"""Test R12 (pre-registered, diary 2026-10-16): does releasing the tracer and oxygen after 12 cycles instead of 80
change the L8 scalar QoIs? Same binary, level, geometry and theta as the T5 runs; only n_mix_cycles (12) and t_end
((12 + 150) cycles) differ. Compared against the T5 run at the same rpm (run ids in test_t5_roughness.json).
Usage: uv run python scripts/submit_r12.py [--dry-run]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm  # noqa: E402
from scripts.submit_l8_t5_t6 import NTASKS, POST, SPINUP, WALL, params  # noqa: E402

SPIN_R12 = 12
RPMS = (20, 21.25, 22.5, 23.75, 25, 26.25, 27.5, 28.75, 30, 31.25, 31.875, 32.5)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for rpm in RPMS:
        p = params(rpm, 7.0, f"r12_l8_rpm{rpm:g}_th7")
        p["t_end"] = round(p["t_end"] / (SPINUP + POST) * (SPIN_R12 + POST), 4)
        p["n_mix_cycles"] = SPIN_R12
        print(f"{p['run_id']}: n_mix_cycles {p['n_mix_cycles']}, t_end {p['t_end']}")
        if not a.dry_run:
            job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=WALL,
                               template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, ntasks=NTASKS, mem="4G")
            print("  job", job)


if __name__ == "__main__":
    main()
