"""L10 at the three Fig 13(b) training angles, cross-level from the L9 states.

Why. The 13(b) MF borrows its HF grid uncertainty, U(L9)/q, from the rpm sweep
at 7 deg and assumes it at every angle (plot_mf_fig13_unc.py). With L10 at
2, 4 and 7 deg, U is measured at those angles from L8/L9/L10 (mf_l8_th* is the
L8 side). L10 at 7 deg also checks the 7-deg-from-rpm-sweep value directly.

Seeds. Each restarts from fig10_l9_th{th} (80 cold-start cycles, lean binary)
and refines onto L10 at the SAME condition, with every *_prev equal to the
current value, so the restart ramp has nothing to interpolate. The
seeds are cold starts rather than chained states, which is why N_CYCLES is 8,
not the 6 used by l10_fig13a_xlevel: 2 cycles to settle after the refine, then 6 settled.

Cost. 58.2 min/cycle measured at (L10, 32 ranks, 32.5 rpm, 7 deg). Smaller
angles move less fluid and take larger steps, so this is an upper bound. The
graceful walltime checkpoint (f1c11e0) covers an underestimate.

Usage:
    uv run python scripts/submit_mf_l10_angle.py --dry-run
    uv run python scripts/submit_mf_l10_angle.py
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm               # noqa: E402
from scripts.dump_fields import fields as _dump_fields  # noqa: E402

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-f1c11e0"
RPM, NTASKS, N_CYCLES, MIN_PER_CYCLE, SAFETY = 32.5, 32, 8, 58.2, 1.3
ANGLES = [2.0, 4.0, 7.0]
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}


def t_scales(theta: float) -> tuple[float, float]:
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    return T_per, L / (V / (H * 0.5) / T_per)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    hours = N_CYCLES * MIN_PER_CYCLE / 60.0 * SAFETY
    walltime = f"{int(hours) + 1:02d}:00:00"
    omega_b = RPM * 2 * math.pi / 60.0
    zeros = [0.0, 0.0, 0.0]
    for th in ANGLES:
        ckpt = ROOT / "runs" / f"fig10_l9_th{th:g}" / "checkpoint.dump"
        t_ck = _dump_fields(ckpt)[0]["t"]
        if t_ck <= 0:
            raise SystemExit(f"{ckpt} reports t={t_ck}; cannot arm a restart")
        T_per, T_bio = t_scales(th)
        theta = [th, 0.0, 0.0]
        params = {
            "run_id": f"mf_l10_th{th:g}", "fidelity": 10, "t_checkpoint": t_ck,
            "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
            "theta_max": theta, "phi_angular": zeros,
            "omega_b": omega_b, "omega_h": 0.0,
            "amplitude_h": zeros, "phi_horizontal": zeros,
            "omega_b_prev": omega_b, "theta_max_prev": theta,
            "phi_angular_prev": zeros, "amplitude_h_prev": zeros,
            "phi_horizontal_prev": zeros, "omega_h_prev": 0.0,
            "t_end": round(N_CYCLES * T_per / T_bio, 4),
            "n_mix_cycles": N_CYCLES + 10,   # no tracer/oxygen release
            "_binary": BINARY,
        }
        print(f"  {params['run_id']}: seed t={t_ck:.4f}, {N_CYCLES} cycles, "
              f"t_end={params['t_end']}, walltime {walltime}")
        if a.dry_run:
            continue
        job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                           walltime=walltime,
                           template=ROOT / "config" / "slurm_mpi_template.sh",
                           cpus=1, ntasks=NTASKS, mem="4G", checkpoint=str(ckpt))
        print(f"    submitted job={job}")


if __name__ == "__main__":
    main()
