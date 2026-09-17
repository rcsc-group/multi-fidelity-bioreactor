"""One run for Kim's Figs. 2-7: the reference condition, fully instrumented.

Figs. 2, 3, 4, 5, 6 and 7 are all the SAME condition -- 7 deg, 32.5 rpm -- seen
six ways: the rocking profile, the instantaneous vortical structure, the steady
streaming, the degree of mixing for four initial tracer configurations, those
tracer fields, and oxygen transfer. Six runs would differ only in noise, so
this submits one and records everything the six need:

  * the signed spatially averaged velocities (normf.dat, always on)   -> Fig 2
  * field snapshots over the release cycles                           -> Figs 3, 6, 7
  * the steady-streaming field (streaming_field.bin, always on)       -> Fig 4
  * all four tracer configurations (EXTRA_TRACERS=1 binary)           -> Figs 5, 6

The snapshot window spans the tracer/oxygen release through the mixing that
follows, because Fig. 3 wants four phases of the release cycle while Figs. 6
and 7 want instants defined by thresholds -- chi = 0.5, C* = 0.10/0.25/0.50 --
that nothing knows in advance. Recording the window and picking the frame in
analysis keeps every threshold out of the solver.

Level. L9 rather than Kim's L10: the figures here are structural (where the
vortices sit, how a tracer folds, what the oxygen front looks like) rather than
the converged scalars Figs. 9-13 report, and L10 would cost 13x for the same
picture.

Cadence. Five frames per period, not thirteen. Each snapshot is ten planes --
~10 MB at L9 -- so the window is a gigabyte-scale decision, and five is
enough: the interval is T_p/(N + 0.618), so frames drift across phase instead
of pinning to N points, and 150 frames over 31 cycles cover phase densely.
Fig. 3's four quarter-cycle phases are matched on PHASE rather than on
absolute time, which a settled periodic flow makes equivalent and which is
what lets the cadence be this low.

Usage:
    uv run python scripts/submit_figset.py --dry-run
    uv run python scripts/submit_figset.py --level 10
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.autoextend import walltime_safety    # noqa: E402
from scripts.cost_model import min_per_cycle      # noqa: E402
from scripts.simulate import submit_slurm         # noqa: E402

# EXTRA_TRACERS=1 build: Figs 5 and 6 need all four initial configurations.
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-figset-3645d13"
RPM = 32.5
SPINUP = 80                 # Kim's release instant, t/T_p = 80
SNAP_START, SNAP_END = 79.0, 110.0
FRAMES_PER_PERIOD = 5
MIX_CYCLES_AFTER = 40       # enough for chi = 0.5 on every configuration
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
THETA = [7.0, 0.0, 0.0]
NTASKS = {6: 8, 7: 8, 8: 8, 9: 32, 10: 32}


def t_scales() -> tuple[float, float]:
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA[0])))
    return T_per, L / (V / (H * 0.5) / T_per)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, default=9)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    ntasks = NTASKS[a.level]
    T_per, T_bio = t_scales()
    cycles = SPINUP + MIX_CYCLES_AFTER
    mpc, conf = min_per_cycle(level=a.level, ntasks=ntasks, rpm=RPM)
    hours = cycles * mpc / 60.0 * walltime_safety(conf)
    if hours > 47:
        raise SystemExit(
            f"L{a.level} needs {hours:.1f} h for {cycles} cycles, over the cap "
            f"-- segment it or drop a level.")
    params = {
        "run_id": f"figset_l{a.level}", "fidelity": a.level,
        "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
        "theta_max": THETA, "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": RPM * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": round(cycles * T_per / T_bio, 4),
        "n_mix_cycles": SPINUP,
        "snap_start_cycle": SNAP_START,
        "snap_end_cycle": SNAP_END,
        "frames_per_period": FRAMES_PER_PERIOD,
        "_binary": BINARY,
    }
    walltime = f"{int(hours) + 1:02d}:00:00"
    n_frames = int((SNAP_END - SNAP_START) * FRAMES_PER_PERIOD)
    print(f"Fig 2-7 set: L{a.level}, {ntasks} ranks, {RPM:g} rpm, "
          f"{THETA[0]:g} deg, {cycles} cycles")
    print(f"  t_end={params['t_end']}, walltime={walltime} "
          f"({mpc:.2f} min/cyc, {conf})")
    print(f"  snapshots over t/T_p {SNAP_START:g}-{SNAP_END:g} "
          f"-> ~{n_frames} frames, ~{n_frames * 10 * (2 ** a.level) ** 2 * 4 / 1e9:.1f} GB")
    if a.dry_run:
        print("  DRY RUN -- not submitted")
        return
    job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                       walltime=walltime,
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=ntasks, mem="4G")
    rec = ROOT / "logs" / f"submitted_{job}.txt"
    rec.parent.mkdir(exist_ok=True)
    rec.write_text(f"{params['run_id']}\n")
    print(f"  submitted {params['run_id']}  job={job}")


if __name__ == "__main__":
    main()
