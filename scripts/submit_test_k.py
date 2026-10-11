"""Test K (pre-registered, diary 2026-10-10): diffusivity scan. Is dtmix (and kLa) set by the grid or by physics?

32.5 rpm, theta 7, kmix protocol (80-cycle spin-up, release, 150 cycles), lean binary with runtime D multipliers
(commit fab19eb). tracer_D_scale = oxy_D_scale = s.
  L6, L7, L8: s in {1, 1e1, 1e2, 1e3, 1e4, 1e5}; L9: s in {1, 1e3, 1e5}.
Run ids: k_l<L>_s<s as 1eK>.

Usage: uv run python scripts/submit_test_k.py [--dry-run] [--only RUN ...]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm  # noqa: E402
from scripts.submit_l8_t5_t6 import params  # noqa: E402

BINARY = "/oscar/data/dharri15/eaguerov/bin/BioReactor-mpi-lean-dscale-fab19eb"
RPM = 32.5
SCALES = {6: (0, 1, 2, 3, 4, 5), 7: (0, 1, 2, 3, 4, 5), 8: (0, 1, 2, 3, 4, 5), 9: (0, 3, 5)}
RESOURCES = {6: (4, "04:00:00"), 7: (8, "06:00:00"), 8: (16, "08:00:00"), 9: (32, "30:00:00")}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", nargs="+")
    a = ap.parse_args()
    for level, exps in SCALES.items():
        ntasks, wall = RESOURCES[level]
        for k in exps:
            rid = f"k_l{level}_s1e{k}"
            if a.only and rid not in a.only:
                continue
            p = params(RPM, 7.0, rid)
            p.update(fidelity=level, tracer_D_scale=10.0 ** k, oxy_D_scale=10.0 ** k, _binary=BINARY)
            print(f"  {rid}: s={10.0 ** k:g}, {ntasks} ranks, {wall}")
            if a.dry_run:
                continue
            job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=wall,
                               template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, ntasks=ntasks, mem="4G")
            print(f"    job={job}")


if __name__ == "__main__":
    main()
