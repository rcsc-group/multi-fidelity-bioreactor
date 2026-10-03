"""Replicate and spin-up-length runs for the noise model of the MF problem definition.

Same protocol and binary as the existing ladder (kmix_l6/l8 cold starts,
fig9_l9 cold starts): cold start, spin-up of `release` cycles, tracer and
oxygen released at that cycle, then POST cycles. Only the release cycle
changes, so the spread across release cycles measures how sensitive the QoI
is to the flow state at release (the noise s(x, h) of the model).

Uses:
  replicates     --release 82 85   (with the existing release-80 runs: 3 per point)
  spin-up test   --release 33 50   at L8, 32.5 rpm: does the spin-up length change
                 dtmix? The L10 point fig9_l10b_seg2 released after 33 cycles.

Usage: uv run python scripts/submit_replicates.py --level 8 --rpms 17.5 25 32.5 \
           --release 82 85 [--dry-run]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm   # noqa: E402

GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
THETA = 7.0
LEAN = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-f1c11e0"   # kmix_l6/l7/l8 binary
PROD = "/oscar/scratch/eaguerov/BioReactor-mpi-prod-4b3a435"   # fig9_l9 binary
# level -> (ranks, walltime, binary, post-release cycles)
LEVELS = {6: (4, "02:00:00", LEAN, 150), 7: (8, "04:00:00", LEAN, 150),
          8: (16, "08:00:00", LEAN, 150), 9: (32, "30:00:00", PROD, 160)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, required=True, choices=sorted(LEVELS))
    ap.add_argument("--rpms", type=float, nargs="+", required=True)
    ap.add_argument("--release", type=int, nargs="+", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    ntasks, wall, binary, post = LEVELS[a.level]
    if not Path(binary).exists():
        sys.exit(f"binary missing: {binary}")
    zeros = [0.0, 0.0, 0.0]
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    for rpm in a.rpms:
        omega_b = rpm * 2 * math.pi / 60.0
        T_per = 2 * math.pi / omega_b
        V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA)))
        T_nd = T_per / (L / (V / (H * 0.5) / T_per))
        for rel in a.release:
            run_id = f"rep_l{a.level}_rpm{rpm:g}_rel{rel}"
            if (ROOT / "runs" / run_id).exists():
                print(f"  {run_id}: exists, skipped")
                continue
            p = {"run_id": run_id, "fidelity": a.level, "geometry": GEOMETRY,
                 "fill_level": 0.5, "n_harmonics": 1,
                 "theta_max": [THETA, 0.0, 0.0], "phi_angular": zeros,
                 "omega_b": omega_b, "omega_h": 0.0, "amplitude_h": zeros,
                 "phi_horizontal": zeros, "n_mix_cycles": rel,
                 "t_end": round((rel + post) * T_nd, 4), "_binary": binary}
            print(f"  {run_id}: release at cycle {rel}, {post} cycles after, t_end={p['t_end']}")
            if a.dry_run:
                continue
            job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=wall,
                               template=ROOT / "config" / "slurm_mpi_template.sh",
                               cpus=1, ntasks=ntasks, mem="4G")
            print(f"    job={job}")


if __name__ == "__main__":
    main()
