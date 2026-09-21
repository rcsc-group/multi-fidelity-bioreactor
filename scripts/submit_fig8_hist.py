"""Re-record Kim Fig. 8(b)/(c) at L10 with enough distinct phases to find the peak.

Panels (b) and (c) are the distributions of shear stress and EDR across the
liquid AT THE INSTANT their spatial means peak. Choosing that instant is the
whole measurement, and the existing L10 source (runs/57f68830) cannot do it:
its frame interval is exactly T_p/5, so however long it runs it samples five
phases and nothing else. The peak it offers is the largest of five, roughly a
tenth below the true one, and the histogram is drawn one phase away from where
Kim drew his.

That cadence bug is fixed in the solver -- frames land at T_p/(N + 0.618),
the golden-ratio offset being the interval hardest to approximate by a
rational, so successive frames spread over phase instead of pinning to N
points (BioReactor.c, event movies_output_tau). This run is the first L10
recording that uses it: 13 frames per period over 10 periods is ~130 distinct
phases rather than 5.

Warm-started from that same converged state with restart_continue=1, so
nothing is re-injected and no spin-up is repaid -- the panels need the settled
flow field, not a tracer protocol.

Cost: L10 at 32 ranks is 58.2 min/cycle without video and video costs a
further ~43 (cost_model.py), so ~17 h for the ten cycles. It needs 32 CPUs
free, which the Fig 9 sweep currently does not leave.

Usage:
    uv run python scripts/submit_fig8_hist.py --dry-run
    uv run python scripts/submit_fig8_hist.py
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.dump_fields import fields as dump_fields   # noqa: E402
from scripts.simulate import submit_slurm               # noqa: E402

# VIDEOS=1 build of 25dbf66. The older video binaries on scratch all predate
# the off-period cadence (ff57efe), the bag-mask fix (d566037) and the
# Kim-exact tau/EDR columns (4b3a435) -- recording with one of those would
# reproduce the very sampling and masking defects this run exists to remove.
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-video-f1c11e0"
SEED_RUN = "57f68830"          # converged L10 at 7 deg, 32.5 rpm
# [2026-09-21] The first re-record (fig8_hist_l10) never settled. Measured on
# its own frames, the SECOND harmonic of <tau'_w> -- which half-period
# antisymmetry forbids in a settled signal -- decays 0.880 -> 0.440 -> 0.228
# over cycles 0-4, 4-7, 7-11 after the restart, i.e. halving every ~3.5
# cycles and still 0.228 at the end. Its waveform is a relaxing transient,
# not L10's shape, which is why the curve read as broken rather than
# sinusoidal and very likely why its EDR peaked 1.7x Kim's.
#
# A continuation re-arms the N_RAMP_CYCLES = 3 soft start, so each restart
# buys a fresh transient regardless of how settled its seed is. The fix is
# therefore length, not a better seed: run long enough that the new
# transient decays below the level that matters before the recording window.
# At a 3.5-cycle half-life, 0.88 -> under 0.05 takes ~15 cycles, so 20 gives
# margin and the analysis uses the tail.
RPM, LEVEL, NTASKS = 32.5, 10, 32
N_CYCLES = 20            # ~15 to settle, tail to record
FRAMES_PER_PERIOD = 13
MIN_PER_CYCLE_VIDEO = 101.2    # 58.2 measured + ~43 video (cost_model.py)
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
THETA = [7.0, 0.0, 0.0]


def t_scales() -> tuple[float, float]:
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA[0])))
    return T_per, L / (V / (H * 0.5) / T_per)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cycles", type=float, default=N_CYCLES)
    ap.add_argument("--seed", default=SEED_RUN,
                    help="run whose checkpoint.dump to continue from")
    ap.add_argument("--run-id", default="fig8_hist_l10")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    ckpt = ROOT / "runs" / a.seed / "checkpoint.dump"
    if not ckpt.exists():
        raise SystemExit(f"checkpoint not found: {ckpt}")
    # Gated on params.t_checkpoint > 0, not on the dump being staged: without
    # it the binary silently cold-starts (see submit_fig9_l10_anchor.py).
    t_ck = dump_fields(ckpt)[0]["t"]
    if t_ck <= 0:
        raise SystemExit(f"checkpoint reports t={t_ck}; cannot arm a restart")

    T_per, T_bio = t_scales()
    params = {
        "run_id": a.run_id, "fidelity": LEVEL,
        "t_checkpoint": t_ck,
        "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
        "theta_max": THETA, "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": RPM * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": round(a.cycles * T_per / T_bio, 4),
        "n_mix_cycles": 1,
        "restart_continue": 1,
        "frames_per_period": FRAMES_PER_PERIOD,
        "_binary": BINARY,
        "_parent_run": a.seed,
    }
    hours = a.cycles * MIN_PER_CYCLE_VIDEO / 60.0 * 1.4
    # Clamped to the 47 h Exploratory cap. Measured on the first re-record:
    # 10 cycles took 18.7 h, i.e. 112 min/cycle, so 20 cycles is ~37 h and the
    # clamp is not binding. If the estimate is wrong the run now checkpoints
    # before the kill (BIOREACTOR_MAXRUNTIME_S) instead of losing itself.
    walltime = f"{min(int(hours) + 1, 47):02d}:00:00"
    print(f"Fig 8(b,c) L10 re-record: {a.cycles:g} cycles, "
          f"{FRAMES_PER_PERIOD} frames/period "
          f"(~{a.cycles * FRAMES_PER_PERIOD:.0f} distinct phases)")
    # t_end is a RELATIVE duration -- BioReactor.c:542 forms
    # t_end_abs = t_checkpoint + t_end -- so print both. Printing only the
    # relative value next to the restart instant reads like the run ends
    # before it starts.
    print(f"  from {a.seed} at t={t_ck:.4f}, +{params['t_end']} "
          f"-> ends at t={t_ck + params['t_end']:.4f}, walltime={walltime}")
    if a.dry_run:
        print("  DRY RUN -- not submitted")
        return
    job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                       walltime=walltime,
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=NTASKS, mem="4G", checkpoint=str(ckpt))
    rec = ROOT / "logs" / f"submitted_{job}.txt"
    rec.parent.mkdir(exist_ok=True)
    rec.write_text(f"{params['run_id']}\n")
    print(f"  submitted {params['run_id']}  job={job}")


if __name__ == "__main__":
    main()
