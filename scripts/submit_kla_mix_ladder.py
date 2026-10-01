"""L6 / L8 rpm sweeps with tracer + oxygen, for an E&H grid extrapolation at every rpm.

With L7 (submit_l7_kla_mix.py) and L9 (fig9_l9_rpm*) this gives L6-L9 at all ten
rpm, so Kim's kLa and dtmix can be extrapolated to h -> 0 across the sweep
(grid_ladder_32p5.png did it at one rpm). Cold starts with Kim's protocol, the
same as the L9 sweep: 80-cycle spin-up, release at cycle 80, then POST cycles.
Binary f1c11e0, the same as L7/L9.

Usage:  uv run python scripts/submit_kla_mix_ladder.py --level 6|8 [--dry-run]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm   # noqa: E402

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-f1c11e0"
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
SPINUP, POST, THETA = 80, 150, 7.0
LEVELS = {6: (4, "02:00:00"), 8: (16, "06:00:00")}   # ranks, walltime


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, required=True, choices=sorted(LEVELS))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    ntasks, wall = LEVELS[a.level]
    zeros = [0.0, 0.0, 0.0]
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    for rpm in RPMS:
        omega_b = rpm * 2 * math.pi / 60.0
        T_per = 2 * math.pi / omega_b
        V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA)))
        T_nd = T_per / (L / (V / (H * 0.5) / T_per))
        p = {"run_id": f"kmix_l{a.level}_rpm{rpm:g}", "fidelity": a.level,
             "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
             "theta_max": [THETA, 0.0, 0.0], "phi_angular": zeros,
             "omega_b": omega_b, "omega_h": 0.0, "amplitude_h": zeros,
             "phi_horizontal": zeros, "n_mix_cycles": SPINUP,
             "t_end": round((SPINUP + POST) * T_nd, 4), "_binary": BINARY}
        print(f"  {p['run_id']}: {SPINUP}+{POST} cycles, t_end={p['t_end']}")
        if a.dry_run:
            continue
        job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=wall,
                           template=ROOT / "config" / "slurm_mpi_template.sh",
                           cpus=1, ntasks=ntasks, mem="4G")
        print(f"    job={job}")


if __name__ == "__main__":
    main()
