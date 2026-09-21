"""Surface-elevation sweep for Kim's Figs. 14 and 15.

Fig. 14 is the maximum surface elevation at the left end of the bag against
rpm and against angle; Fig. 15 is the frequency spectrum of that elevation at
three of each. Both read `posY_left`, a column the solver only started writing
on 2026-09-16, so no existing run carries it and the Fig. 9 and 10 sweeps
cannot supply it retroactively.

They do not need to. Surface elevation is a settled-state property of the
forcing: no tracer, no oxygen, no 80-cycle protocol, just enough cycles to
shed the soft start and resolve a spectrum. Sixteen points at 20 cycles each
is ~23 h in total at L9 -- against ~1200 cycles for one Fig. 9 point.

`n_mix_cycles` is set past `t_end` on purpose: the tracer and oxygen events
then never fire, so these runs pay for neither.

Usage:
    uv run python scripts/submit_fig14.py --dry-run
    uv run python scripts/submit_fig14.py --rpms 32.5 --angles
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

# Lean build (one tracer): these runs release nothing, so the extra tracer
# configurations would be advected for nothing.
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-f1c11e0"
LEVEL, NTASKS = 9, 32
N_CYCLES = 20                # 5 discarded as soft start, 15 for the spectrum
RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
ANGLES = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
REF_RPM, REF_ANGLE = 32.5, 7.0
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}


def t_scales(rpm: float, theta_deg: float) -> tuple[float, float]:
    omega_b = rpm * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta_deg)))
    return T_per, L / (V / (H * 0.5) / T_per)


def submit(run_id: str, rpm: float, theta: float, cycles: float,
           dry: bool) -> None:
    T_per, T_bio = t_scales(rpm, theta)
    mpc, conf = min_per_cycle(level=LEVEL, ntasks=NTASKS, rpm=rpm)
    hours = cycles * mpc / 60.0 * walltime_safety(conf)
    params = {
        "run_id": run_id, "fidelity": LEVEL,
        "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
        "theta_max": [theta, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": rpm * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": round(cycles * T_per / T_bio, 4),
        "n_mix_cycles": int(cycles) + 10,     # past t_end: never released
        "_binary": BINARY,
    }
    walltime = f"{min(int(hours) + 1, 47):02d}:00:00"
    print(f"  {run_id:<22} {rpm:>5.1f} rpm  {theta:g} deg  "
          f"t_end={params['t_end']:8.4f}  walltime={walltime}")
    if dry:
        print("       DRY RUN -- not submitted")
        return
    job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                       walltime=walltime,
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=NTASKS, mem="4G")
    rec = ROOT / "logs" / f"submitted_{job}.txt"
    rec.parent.mkdir(exist_ok=True)
    rec.write_text(f"{run_id}\n")
    print(f"       submitted {run_id}  job={job}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rpms", type=float, nargs="*", default=RPMS)
    ap.add_argument("--angles", type=float, nargs="*", default=ANGLES)
    ap.add_argument("--cycles", type=float, default=N_CYCLES)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    print(f"Fig 14/15 surface-elevation sweep: L{LEVEL}, {NTASKS} ranks, "
          f"{a.cycles:g} cycles")
    for rpm in a.rpms:
        submit(f"fig14_l{LEVEL}_rpm{rpm:g}", rpm, REF_ANGLE, a.cycles, a.dry_run)
    for theta in a.angles:
        # theta = 7 at 32.5 rpm is the same condition as the rpm sweep's own
        # 32.5 point; submitting it twice would pay for one number twice.
        if abs(theta - REF_ANGLE) < 1e-9 and REF_RPM in a.rpms:
            print(f"  fig14_l{LEVEL}_th{theta:g}: same condition as "
                  f"rpm{REF_RPM:g} -- skipped")
            continue
        submit(f"fig14_l{LEVEL}_th{theta:g}", REF_RPM, theta, a.cycles,
               a.dry_run)


if __name__ == "__main__":
    main()
