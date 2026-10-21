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
N_CYCLES, NTASKS, MARGIN, SPINUP_POST = 12, 32, 1.2, 230
# BINARY predates f1c11e0 (deadline checkpoint): a job killed at walltime leaves NO dump.
_RATE_32 = {17.5: 194.0, 25: 162.0, 32.5: 51.0, 37.5: 42.0}  # min per L10 cycle, 32 ranks; 17.5-32.5 measured on x8_l10 (2026-10-12)


def source_run(rpm: float, theta: float) -> str:
    for name in ([f"kmix_l8_rpm{rpm:g}"] if theta == 7.0 else []) + \
                [f"t6_l8_rpm{rpm:g}_th{theta:g}", f"t5_l8_rpm{rpm:g}_th{theta:g}"]:
        if (ROOT / "runs" / name / "checkpoint.dump").exists():
            return name
    raise FileNotFoundError(f"no settled L8 run with a checkpoint at rpm {rpm:g}, theta {theta:g}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--points", nargs="+", required=True)
    ap.add_argument("--cycles", type=int, default=N_CYCLES, help="L10 cycles after the restart (12 after 2026-10-12)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--exclude", default=None, help="nodes to avoid (sbatch --exclude)")
    ap.add_argument("--level", type=int, default=10, help="target level (9 or 10); 9 for the H2 middle level")
    ap.add_argument("--rate", type=float, default=None, help="measured min per cycle (overrides the theta-7 table)")
    ap.add_argument("--post", type=int, default=0, help="extra cycles AFTER the tracer/oxygen release (release stays at --cycles); one job does warm start, release and the scalar run")
    ap.add_argument("--tag", default="x8", help="run id prefix")
    ap.add_argument("--account", default=None, help="sbatch account override (mbessa-condo only for a named experiment the user authorised)")
    a = ap.parse_args()
    for pt in a.points:
        rpm, theta = map(float, pt.split(":"))
        src = source_run(rpm, theta)
        sd = ROOT / "runs" / src
        t_dump = float(np.loadtxt(sd / "shear_stress.dat", skiprows=1)[-1, 1])
        p = l8_params(rpm, theta, f"{a.tag}_l{a.level}_rpm{rpm:g}_th{theta:g}")
        t_nd = p["t_end"] / SPINUP_POST  # nondimensional period (depends on theta through the fill volume)
        zeros = [0.0, 0.0, 0.0]
        p.update(fidelity=a.level, t_end=round((a.cycles + a.post) * t_nd, 6), n_mix_cycles=a.cycles, t_checkpoint=t_dump,
                 omega_b_prev=p["omega_b"], theta_max_prev=[theta, 0.0, 0.0], phi_angular_prev=zeros,
                 amplitude_h_prev=zeros, phi_horizontal_prev=zeros, omega_h_prev=0.0, _binary=BINARY)
        xs = sorted(_RATE_32)
        rate = float(np.interp(min(max(rpm, xs[0]), xs[-1]), xs, [_RATE_32[x] for x in xs]))
        rate = rate if a.level == 10 else rate / 6.0  # L9 ~1/8 of L10 per cycle; /6 kept conservative
        rate = a.rate if a.rate else rate  # theta 4 runs ~2x slower than the theta-7 table (x8_l10_rpm32.5_th4: 97)
        hours = min(48, max(4, math.ceil((a.cycles + a.post) * rate * MARGIN / 60.0)))  # 48 h cap; no dump at kill
        print(f"  {p['run_id']}: source {src} (t {t_dump:g}), t_end {p['t_end']}, walltime {hours} h")
        if a.dry_run:
            continue
        job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=f"{hours:02d}:00:00",
                           template=ROOT / "config" / "slurm_mpi_template.sh",
                           checkpoint=str((sd / "checkpoint.dump").resolve()), cpus=1, ntasks=NTASKS, mem="4G",
                           exclude=a.exclude, account=a.account)
        print(f"    job={job}")


if __name__ == "__main__":
    main()
