"""Perturbation-size scan at L8 (pre-registered, diary 2026-10-10, test P): does the jump in dtmix shrink with the size of
the rpm perturbation (smooth but steep structure) or not (sensitive dependence: effectively noise)?

Same protocol and binary as kmix_l8 / T5 (scripts/submit_l8_t5_t6.params); only rpm changes.
Base rpm 22.5 and 30 (T5 measured +-0.1 rpm jumps of 0.20-0.35 in log there); offsets +1e-2, +1e-3, +1e-6 rpm.
Run ids: p8_rpm<base>_d<offset as 1e-k>_th7.

Usage: uv run python scripts/submit_l8_perturb.py [--dry-run]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm  # noqa: E402
from scripts.submit_l8_t5_t6 import NTASKS, WALL, params  # noqa: E402

BASES = (22.5, 30.0)
OFFSETS = {"1e-2": 1e-2, "1e-3": 1e-3, "1e-6": 1e-6}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for base in BASES:
        for tag, d in OFFSETS.items():
            rid = f"p8_rpm{base:g}_d{tag}_th7"
            p = params(base + d, 7.0, rid)
            print(f"  {rid}: rpm {base + d!r}, omega_b {p['omega_b']!r}")
            if a.dry_run:
                continue
            job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=WALL,
                               template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, ntasks=NTASKS, mem="4G")
            print(f"    job={job}")


if __name__ == "__main__":
    main()
