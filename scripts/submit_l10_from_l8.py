"""L10 hydrodynamics warm-started from a settled L8 run (two-level cross-level restart).

Same mechanics as scripts/submit_l10_fig13a_xlevel.py (which restarts L10 from L9): restart from the source
run's checkpoint.dump, refine() onto the L10 grid, run N_CYCLES, release at the end; tau/EDR are read over the
last 3 cycles before release (hydro_dataset_v2 window). Only the source level differs (L8, not L9).
Sources: kmix_l8_rpm<r> (theta 7) or t6_l8_rpm<r>_th<theta> / t5_l8_rpm<r>_th7.

Validation first (pre-registered, diary 2026-10-10): theta 7, rpm 17.5/25/32.5, compared with
l10_fig13a_xlevel_rpm<r> (L10 warm-started from L9).

Usage: uv run python scripts/submit_l10_from_l8.py --points RPM:THETA [RPM:THETA ...] [--dry-run]
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm  # noqa: E402
from scripts.submit_l8_t5_t6 import params as l8_params  # noqa: E402

BINARY = "/oscar/data/dharri15/eaguerov/bin/BioReactor-mpi-xlevel-chain"
N_CYCLES, NTASKS, MARGIN, SPINUP_POST = 6, 32, 1.5, 230
_RATE_32 = {17.5: 168.1, 32.5: 101.0, 37.5: 42.0}  # min per L10 cycle at 32 ranks, scripts/cost_model.py


def source_run(rpm: float, theta: float) -> str:
    for name in ([f"kmix_l8_rpm{rpm:g}"] if theta == 7.0 else []) + \
                [f"t6_l8_rpm{rpm:g}_th{theta:g}", f"t5_l8_rpm{rpm:g}_th{theta:g}"]:
        if (ROOT / "runs" / name / "checkpoint.dump").exists():
            return name
    raise FileNotFoundError(f"no settled L8 run with a checkpoint at rpm {rpm:g}, theta {theta:g}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", nargs="+", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for pt in a.points:
        rpm, theta = map(float, pt.split(":"))
        src = source_run(rpm, theta)
        sd = ROOT / "runs" / src
        t_dump = float(np.loadtxt(sd / "shear_stress.dat", skiprows=1)[-1, 1])
        p = l8_params(rpm, theta, f"x8_l10_rpm{rpm:g}_th{theta:g}")
        t_nd = p["t_end"] / SPINUP_POST  # nondimensional period (depends on theta through the fill volume)
        zeros = [0.0, 0.0, 0.0]
        p.update(fidelity=10, t_end=round(N_CYCLES * t_nd, 6), n_mix_cycles=N_CYCLES, t_checkpoint=t_dump,
                 omega_b_prev=p["omega_b"], theta_max_prev=[theta, 0.0, 0.0], phi_angular_prev=zeros,
                 amplitude_h_prev=zeros, phi_horizontal_prev=zeros, omega_h_prev=0.0, _binary=BINARY)
        xs = sorted(_RATE_32)
        rate = float(np.interp(min(max(rpm, xs[0]), xs[-1]), xs, [_RATE_32[x] for x in xs]))
        hours = max(4, math.ceil(N_CYCLES * rate * MARGIN / 60.0))
        print(f"  {p['run_id']}: source {src} (t {t_dump:g}), t_end {p['t_end']}, walltime {hours} h")
        if a.dry_run:
            continue
        job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=f"{hours:02d}:00:00",
                           template=ROOT / "config" / "slurm_mpi_template.sh",
                           checkpoint=str((sd / "checkpoint.dump").resolve()), cpus=1, ntasks=NTASKS, mem="4G")
        print(f"    job={job}")


if __name__ == "__main__":
    main()
