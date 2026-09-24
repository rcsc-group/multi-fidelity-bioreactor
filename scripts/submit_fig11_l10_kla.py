"""L10 kLa points for Kim Fig. 11 at the two cheapest rpms, 35 and 37.5.

Why these two. L10 per-cycle cost falls ~4x from 17.5 to 37.5 rpm and the
oxygen crossings come sooner at high rpm, so these are the only points that
fit a week. With the 32.5 rpm L10 anchor (fig9_l10_seg1, lean binary, carries
oxygen) they put three L10 points beside the L9 curve -- enough to show
whether the L9 kLa gap to Kim closes on his own mesh (n_L = 2^10).

Seeds. Each warm-starts from the L10 Fig 13 state at its OWN condition
(l10_fig13a_xlevel_rpm*), whose params already set every *_prev field equal to
the current one. That matters: with *_prev unset a continuation ramps the
forcing up from ZERO over N_RAMP_CYCLES, and fig8_hist_l10 showed that
transient still carrying a 2x/1x harmonic of 0.23 eleven cycles later. With
*_prev == current the ramp has nothing to interpolate. Measured on the seeds
themselves: A2/A1 of the signed wall shear <= 0.026 from the first cycle,
A1 flat to 1% over six cycles.

Protocol deviation, stated: oxygen and tracer release N_SETTLE cycles after
the restart rather than at Kim's absolute cycle 80. The flow has ~250 cycles
of history (a chained rpm ladder); the 18% dtmix shift measured for a cycle-25
release was a cold start releasing into an undeveloped flow, which this is
not.

Usage:
    uv run python scripts/submit_fig11_l10_kla.py --dry-run
    uv run python scripts/submit_fig11_l10_kla.py
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm            # noqa: E402
from scripts.dump_fields import fields as _dump_fields  # noqa: E402

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-f1c11e0"
NTASKS = 32
N_SETTLE = 5
WALLTIME_CAP_H = 47
SAFETY = 1.15
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
THETA = [7.0, 0.0, 0.0]

# post-release cycles: the longer of C*=0.50 on our L9 kLa and 1.3x our L9 dt95
# (37.5) / 1.3x C*=0.25 on Kim's kLa (35), rounded up.
# min/cycle: 37.5 is MEASURED at (10,32) but on a video-writing run, so it is
# an upper bound; 35 is log-interpolated between 32.5 (58.2, no video) and
# 37.5 -- unmeasured. The graceful walltime checkpoint covers an underestimate.
POINTS = {
    37.5: {"post_cycles": 42, "min_per_cycle": 42.0},
    35.0: {"post_cycles": 45, "min_per_cycle": 50.0},
}


def t_scales(rpm: float) -> tuple[float, float]:
    omega_b = rpm * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA[0])))
    return T_per, L / (V / (H * 0.5) / T_per)


def submit(rpm: float, dry: bool) -> None:
    pt = POINTS[rpm]
    seed = ROOT / "runs" / f"l10_fig13a_xlevel_rpm{rpm:g}"
    ckpt = seed / "checkpoint.dump"
    if not ckpt.exists():
        raise SystemExit(f"checkpoint not found: {ckpt}")
    # restart is gated on t_checkpoint > 0, not on a staged dump
    t_ck = _dump_fields(ckpt)[0]["t"]
    if t_ck <= 0:
        raise SystemExit(f"{ckpt} reports t={t_ck}; cannot arm a restart")

    T_per, T_bio = t_scales(rpm)
    cycles = N_SETTLE + pt["post_cycles"]
    hours = cycles * pt["min_per_cycle"] / 60.0 * SAFETY
    walltime = f"{min(int(hours) + 1, WALLTIME_CAP_H):02d}:00:00"
    omega_b = rpm * 2 * math.pi / 60.0
    zeros = [0.0, 0.0, 0.0]
    params = {
        "run_id": f"fig11_l10_rpm{rpm:g}", "fidelity": 10,
        "t_checkpoint": t_ck,
        "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
        "theta_max": THETA, "phi_angular": zeros,
        "omega_b": omega_b, "omega_h": 0.0,
        "amplitude_h": zeros, "phi_horizontal": zeros,
        "omega_b_prev": omega_b, "theta_max_prev": THETA,
        "phi_angular_prev": zeros, "amplitude_h_prev": zeros,
        "phi_horizontal_prev": zeros, "omega_h_prev": 0.0,
        "t_end": round(cycles * T_per / T_bio, 4),
        "n_mix_cycles": N_SETTLE,
        "_binary": BINARY,
    }
    print(f"  {params['run_id']}: seed t={t_ck:.4f}, {cycles} cycles "
          f"({N_SETTLE} settle + {pt['post_cycles']} after release), "
          f"t_end={params['t_end']}, ~{hours:.0f} h -> walltime {walltime}")
    if hours > WALLTIME_CAP_H:
        print("    exceeds the cap: graceful checkpoint + a continuation segment will finish it")
    if dry:
        print("    DRY RUN")
        return
    job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                       walltime=walltime,
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=NTASKS, mem="4G", checkpoint=str(ckpt))
    print(f"    submitted job={job}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rpm", type=float, nargs="*", default=sorted(POINTS))
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for rpm in a.rpm:
        submit(rpm, a.dry_run)


if __name__ == "__main__":
    main()
