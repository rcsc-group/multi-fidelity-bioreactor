"""Velocity fields at one common instant, n_L = 2^6 to 2^10, for Fig. A.16(c).

Panel (c) is the discrete L2 norm of the velocity error against a finer
reference, so it needs the per-cell velocity at ONE instant at every
resolution. Nothing on disk has that: normf.dat carries rms scalars,
frames_tau carries (f, tau, eps), and the Basilisk checkpoint stores a tree
the analysis scripts cannot read. The solver now writes uv_field.bin at the
checkpoint instant (BioReactor.c, write_uv_field), so these runs exist only to
reach a settled state and stop.

Each level cold-starts 30 cycles. Warm-starting a coarse level from a fine
dump is a cross-level restore, which this project has a documented history of
getting silently wrong, and it would in any case give each level a different
starting state -- the one thing a convergence comparison cannot tolerate.

All five stop at the same t/T_p, which the analysis checks before differencing
anything: an L2 norm taken across levels at different phases measures the
phase, not the error.

Cost is dominated by the finest point -- L10 at 32 ranks is 58.2 min/cycle, so
~29 h of the ~32 h total. The other four are hours or minutes.

Usage:
    uv run python scripts/submit_figA16_c.py --dry-run
    uv run python scripts/submit_figA16_c.py --levels 6 7 8
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.autoextend import walltime_safety     # noqa: E402
from scripts.cost_model import min_per_cycle       # noqa: E402
from scripts.simulate import submit_slurm          # noqa: E402

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-uvfield-25dbf66"
RPM = 32.5
N_CYCLES = 30
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
THETA = [7.0, 0.0, 0.0]
# Ranks per level: the coarse grids do not have the cells to keep 32 ranks busy,
# and over-decomposing them costs more in halo exchange than it saves.
NTASKS = {6: 8, 7: 8, 8: 8, 9: 32, 10: 32}


def t_scales() -> tuple[float, float]:
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA[0])))
    return T_per, L / (V / (H * 0.5) / T_per)


def submit(level: int, cycles: float, dry: bool) -> None:
    T_per, T_bio = t_scales()
    ntasks = NTASKS[level]
    mpc, conf = min_per_cycle(level=level, ntasks=ntasks, rpm=RPM)
    hours = cycles * mpc / 60.0 * walltime_safety(conf)
    if hours > 47:
        raise SystemExit(
            f"L{level} needs {hours:.1f} h for {cycles:g} cycles, over the cap "
            f"-- this point has to be segmented, not submitted as one job.")
    params = {
        "run_id": f"figA16c_l{level}", "fidelity": level,
        "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
        "theta_max": THETA, "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": RPM * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": round(cycles * T_per / T_bio, 4),
        "n_mix_cycles": int(cycles) + 1,   # past t_end: no tracer, no oxygen
        "_binary": BINARY,
    }
    walltime = f"{int(hours) + 1:02d}:00:00"
    print(f"  L{level}  {ntasks:2d} ranks  {cycles:g} cyc  "
          f"t_end={params['t_end']:.4f}  walltime={walltime}  "
          f"({mpc:.2f} min/cyc, {conf})")
    if dry:
        print(f"       DRY RUN -- not submitted ({params['run_id']})")
        return
    job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                       walltime=walltime,
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=ntasks, mem="4G")
    rec = ROOT / "logs" / f"submitted_{job}.txt"
    rec.parent.mkdir(exist_ok=True)
    rec.write_text(f"{params['run_id']}\n")
    print(f"       submitted {params['run_id']}  job={job}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", type=int, nargs="*", default=[6, 7, 8, 9, 10])
    ap.add_argument("--cycles", type=float, default=N_CYCLES)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    print(f"Fig A.16(c) field set: {a.cycles:g} cycles at {RPM:g} rpm, "
          f"{THETA[0]:g} deg")
    for level in a.levels:
        submit(level, a.cycles, a.dry_run)


if __name__ == "__main__":
    main()
