"""The 2D/3D pair for Kim's Fig. A.18.

Kim compares a 2D and a 3D simulation of the same condition to justify running
the study in 2D. Two runs, identical in every parameter, differing only in the
grid: quadtree against octree.

The 3D binary needs no source changes -- the embedded bag is a function of x
and y, so it extrudes into a slab under -grid=octree, and the liquid volume
comes out identical (0.286 at both, measured 2026-09-16). What it does need is
CELLS: one level in 3D is 8x the cells of the same level in 2D, on top of the
timestep shrinking, so the level here is deliberately well below the 2D sweeps'.
The figure is a statement about whether out-of-plane motion exists at all, not
a converged scalar, and it is read as a ratio between components of the SAME
run -- which is why a coarse pair still answers it and a mismatched pair would
not.

Both runs are therefore submitted at the SAME level. Comparing a fine 2D run
against a coarse 3D one would confound dimensionality with resolution, which
is the one confusion this figure cannot survive.

Spanwise extent. The domain is a cube, so the 3D bag is as deep as it is wide.
That is not Kim's cellbag, which is deeper than it is tall but not as deep as
it is long. A spanwise-periodic slab of the right depth is the physically
faithful choice and is a modelling decision, not a bug -- recorded here
because the figure's conclusion depends on it and the cube is the convenient
default, not the argued one.

Usage:
    uv run python scripts/submit_figA18.py --dry-run
    uv run python scripts/submit_figA18.py --level 6
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm         # noqa: E402

BINARY_2D = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-af1afce"
BINARY_3D = "/oscar/scratch/eaguerov/BioReactor-mpi-3d-af1afce"
RPM, THETA = 32.5, [7.0, 0.0, 0.0]
N_CYCLES = 32               # 30 to settle, 2 for the window the figure shows
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
# Measured 2026-09-16, serial, fidelity 5 octree: 180 s for one cycle. One
# level is 8x cells and ~2x steps in 3D, so ~16x per level; 32 ranks buy back
# ~20x. These are estimates, not sacct measurements, and the walltime carries
# a matching cushion.
MIN_PER_CYCLE_3D = {5: 0.15, 6: 2.4, 7: 38.0}
MIN_PER_CYCLE_2D = {5: 0.02, 6: 0.05, 7: 0.17}
NTASKS = 32


def t_scales() -> tuple[float, float]:
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA[0])))
    return T_per, L / (V / (H * 0.5) / T_per)


def submit(tag: str, binary: str, level: int, mpc: float, cycles: float,
           dry: bool) -> None:
    T_per, T_bio = t_scales()
    hours = cycles * mpc / 60.0 * 3.0        # estimate, not a measurement
    if hours > 47:
        raise SystemExit(f"{tag} L{level}: {hours:.1f} h, over the cap")
    params = {
        "run_id": f"figA18_{tag}_l{level}", "fidelity": level,
        "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
        "theta_max": THETA, "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": RPM * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": round(cycles * T_per / T_bio, 4),
        "n_mix_cycles": int(cycles) + 10,    # past t_end: nothing released
        "_binary": binary,
    }
    walltime = f"{max(int(hours) + 1, 2):02d}:00:00"
    print(f"  {params['run_id']:<20} L{level}  {cycles:g} cyc  "
          f"t_end={params['t_end']:.4f}  walltime={walltime} "
          f"(~{mpc:g} min/cyc, ESTIMATED)")
    if dry:
        print("       DRY RUN -- not submitted")
        return
    job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                       walltime=walltime,
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=NTASKS, mem="8G")
    rec = ROOT / "logs" / f"submitted_{job}.txt"
    rec.parent.mkdir(exist_ok=True)
    rec.write_text(f"{params['run_id']}\n")
    print(f"       submitted {params['run_id']}  job={job}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, default=6, choices=sorted(MIN_PER_CYCLE_3D))
    ap.add_argument("--cycles", type=float, default=N_CYCLES)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    print(f"Fig A.18 pair: L{a.level}, {NTASKS} ranks, {RPM:g} rpm, "
          f"{THETA[0]:g} deg, {a.cycles:g} cycles")
    submit("2d", BINARY_2D, a.level, MIN_PER_CYCLE_2D[a.level], a.cycles,
           a.dry_run)
    submit("3d", BINARY_3D, a.level, MIN_PER_CYCLE_3D[a.level], a.cycles,
           a.dry_run)


if __name__ == "__main__":
    main()
