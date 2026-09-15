"""L10 Fig13a sweep, cross-level warm-started from the converged L9 sweep.

This is the production use of the 2026-09-14 cross-level fix (diary.md):
each L10 point restarts from the SAME-rpm L9 checkpoint and refine()s onto
the L10 grid, instead of cold-starting L10 and paying the full settling
time at the most expensive resolution. The nine L9 sources are exactly the
runs already trusted for the L9 series in replicated_Fig13.png.

Why per-rpm cross-level rather than chain.py's rpm->rpm chaining: each
point stays independent (no condition-chain contamination), and the nine
jobs are embarrassingly parallel.

Sizing: 32 ranks, pinned to ONE node (slurm_mpi_template.sh sets
--nodes=1; multi-node MPI has an unfixed output-corruption bug, diary.md
2026-09-10 (2), re-confirmed 2026-09-15). Per-cycle cost from
scripts/cost_model.py where measured, otherwise interpolated across the three measured 32-rank anchors.
Walltime carries a 1.5x margin on top, because BioReactor only writes its
checkpoint at the PLANNED END -- a segment killed at walltime loses
everything, so over-reserving is free and under-reserving is total loss.

Usage:
    uv run python scripts/submit_l10_fig13a_xlevel.py [--dry-run]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.simulate import submit_slurm

root = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
DRY = "--dry-run" in sys.argv

# The exact L9 runs backing the L9 series in replicated_Fig13.png.
L9_SOURCES = {
    17.5: "l9_sweep_rpm17.5", 20.0: "20a14369", 22.5: "255e0b87",
    25.0: "bce29aa5", 27.5: "9ec56180", 30.0: "60d09a80",
    32.5: "a34fc4d4", 35.0: "a281a16f", 37.5: "0f0ad3ea",
}

N_CYCLES = 6         # warm-started from a converged L9 state: ~2 to settle + ~4 QSS
NTASKS = 32
MARGIN = 1.5

# 32 ranks, not 64: checked 2026-09-15, ZERO batch nodes had >=64 free cores
# (597 idle but fragmented), so a 64-rank --nodes=1 job would pend exactly
# like the 64h job 6326987 did. 32-rank fits common fragments and, unlike
# 64, is measured at all three anchor rpms.
_RATE_32 = {17.5: 168.1, 32.5: 101.0, 37.5: 42.0}   # scripts/cost_model.py


def rate32(rpm: float) -> float:
    xs = sorted(_RATE_32)
    return float(np.interp(min(max(rpm, xs[0]), xs[-1]), xs, [_RATE_32[x] for x in xs]))


manifest = {}
print(f"{'rpm':>6} {'L9 source':<18} {'min/cyc':>8} {'walltime':>9}  job")
for rpm, src in sorted(L9_SOURCES.items()):
    sd = root / "runs" / src
    t_dump = float(np.loadtxt(sd / "shear_stress.dat", skiprows=1)[-1, 1])
    omega = rpm * 2 * math.pi / 60.0
    r = rate32(rpm)
    hours = max(4, math.ceil(N_CYCLES * r * MARGIN / 60.0))
    walltime = f"{hours:02d}:00:00"
    run_id = f"l10_fig13a_xlevel_rpm{rpm:g}"
    params = {
        "run_id": run_id, "fidelity": 10,
        "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5,
        "n_harmonics": 1, "theta_max": [7.0, 0.0, 0.0],
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": omega, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": N_CYCLES * 0.607329, "n_mix_cycles": N_CYCLES,
        "t_checkpoint": t_dump,
        "omega_b_prev": omega, "theta_max_prev": [7.0, 0.0, 0.0],
        "phi_angular_prev": [0.0, 0.0, 0.0], "amplitude_h_prev": [0.0, 0.0, 0.0],
        "phi_horizontal_prev": [0.0, 0.0, 0.0], "omega_h_prev": 0.0,
        "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-xlevel-chain",
    }
    if DRY:
        job = "DRY"
    else:
        job = submit_slurm(params, project_root=root, runs_root=root / "runs",
                           walltime=walltime,
                           template=root / "config" / "slurm_mpi_template.sh",
                           checkpoint=str((sd / "checkpoint.dump").resolve()),
                           cpus=1, ntasks=NTASKS, mem="4G")
    manifest[run_id] = {"rpm": rpm, "l9_source": src, "min_per_cycle_est": round(r, 1),
                        "walltime": walltime, "n_cycles": N_CYCLES, "job": job}
    print(f"{rpm:6.1f} {src:<18} {r:8.1f} {walltime:>9}  {job}")

if not DRY:
    out = root / "experiments" / "l10_fig13a_xlevel_manifest.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2))
    print(f"\nmanifest -> {out}")
