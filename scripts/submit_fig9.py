"""Submit the Kim Fig. 9 (mixing time vs rocking frequency) replica.

Kim's protocol (Main.tex, Fig. 9 + Fig. 2 caption): the passive tracer is
introduced at t/T_p = 80, and the mixing time is measured from there until
the degree of mixing chi reaches 0.50 / 0.75 / 0.95. So each run must cover
80 spin-up cycles PLUS dtmix_0.95/T_per more -- and dtmix_0.95 grows sharply
as rpm falls (33.6 s at 37.5 rpm, 661.3 s at 15 rpm), which is what makes
this sweep expensive rather than the resolution.

Reference values come from experiments/kimetal2024/csv_raw/
mixing_kla_vs_frequency.csv, which is a bit-exact export of Kim's OWN
published dataset_kimetal2024.xlsx (sheet "Fig. 9,11") -- verified
2026-09-15. Nothing here is digitized.

STAGED ON PURPOSE (diary.md 2026-09-15 (6)). The chi -> dtmix pipeline had
never produced a number anyone compared against Kim: every Fig 13 sweep ran
with t_end < t_mix, so the tracer was never released in any run on disk
(runs/a34fc4d4/tr_oxy.dat is all zeros). The full L7 sweep is ~500
wall-clock hours, so it does not get launched on an unvalidated pipeline.

    --stage validate   ONE L6 point at 32.5 rpm, full Kim protocol (~9 h).
                       Cheap, and still directly comparable to Kim.
    --stage sweep      the 10-point L7 sweep. Only after validate lands and
                       sigma2_max ~ 0.25 and dtmix is sane.

What to check when `validate` finishes
--------------------------------------
  sigma2_max ~ 0.25   the injection invariant (see
                      tests/verification/test_mixing_metrics.py). postprocess
                      returns NaN dtmix if this is off by >20%, because
                      dtmix itself is provably blind to a wrong sigma2_max.
  dtmix_0.50/0.75/0.95 finite, ordered, and within ~2x of Kim at L6.

Usage:
    uv run python scripts/submit_fig9.py --stage validate [--dry-run]
"""
import argparse
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm          # noqa: E402
from scripts.cost_model import min_per_cycle       # noqa: E402

# Built from commit 2736992 at 2026-09-15 12:13 (make build-mpi). Named after
# the commit so a stale-binary mixup is visible rather than inferred -- see
# feedback_binary_deployment.
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-fig9-2736992"
# Same source, built with -DEXTRA_TRACERS=0 (drops the three never-initialised
# tracer variants c/c1/c3 from stracers).
LEAN_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-lean"
# Production binary for the Fig 9 sweep, built from 4b3a435 with
# -DEXTRA_TRACERS=0: lean (1.69x), Kim-exact tau/EDR columns, and the
# steady-streaming vorticity accumulator.
PROD_BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-prod-4b3a435"

KIM_CSV = ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_frequency.csv"
SPINUP_CYCLES = 80     # Kim's tracer release instant, t/T_p = 80
MARGIN = 1.30          # our dtmix may exceed Kim's; don't end the run blind
WALLTIME_SAFETY = 1.6  # cost model is measured pre-injection; tracer gradients cost more

GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
THETA_MAX = [7.0, 0.0, 0.0]


def t_scales(rpm: float) -> tuple[float, float]:
    """(T_per [s], T_bio [s]) -- same definitions as postprocess._t_scales."""
    omega_b = rpm * 2 * math.pi / 60.0
    L = GEOMETRY["a"]
    H = 2 * GEOMETRY["b"]
    th = math.radians(THETA_MAX[0])
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(th))
    U = V / (H * 0.5) / T_per
    return T_per, L / U


def plan(rpm: float, level: int, ntasks: int, allow_missing_cost: bool = False) -> dict:
    """Cycles, non-dimensional t_end and a measured-cost walltime for one point."""
    kim = pd.read_csv(KIM_CSV).set_index("RPM")
    T_per, T_bio = t_scales(rpm)
    cyc_mix = MARGIN * float(kim.loc[rpm, "dtmix_strict_0.95"]) / T_per
    cyc_tot = SPINUP_CYCLES + cyc_mix
    t_end = cyc_tot * T_per / T_bio
    try:
        mpc, conf = min_per_cycle(level=level, ntasks=ntasks, rpm=rpm)
    except ValueError:
        if not allow_missing_cost:
            raise
        mpc, conf = float("nan"), "no cost data; explicit walltime"
    hours = cyc_tot * mpc / 60.0 * WALLTIME_SAFETY
    return {"rpm": rpm, "cyc_tot": cyc_tot, "t_end": t_end,
            "hours": hours, "min_per_cycle": mpc, "confidence": conf,
            "kim_dtmix_0.95": float(kim.loc[rpm, "dtmix_strict_0.95"]),
            "kim_dtmix_0.50": float(kim.loc[rpm, "dtmix_strict_0.5"])}


def submit(rpm: float, level: int, ntasks: int, prefix: str, dry: bool,
           walltime_override: str | None = None,
           binary: str | None = None) -> None:
    p = plan(rpm, level, ntasks, allow_missing_cost=walltime_override is not None)
    if walltime_override is None and p["hours"] > 48:
        raise SystemExit(
            f"L{level} {rpm} rpm needs {p['hours']:.1f} h, over the 48 h "
            f"Exploratory cap -- this point must be chained (scripts/chain.py), "
            f"not submitted as one job.")
    walltime = walltime_override or f"{int(p['hours']) + 1:02d}:00:00"
    run_id = f"{prefix}_l{level}_rpm{rpm:g}"
    params = {
        "run_id": run_id, "fidelity": level,
        "geometry": GEOMETRY, "fill_level": 0.5,
        "n_harmonics": 1, "theta_max": THETA_MAX,
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": rpm * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": round(p["t_end"], 4),
        "n_mix_cycles": SPINUP_CYCLES,
        "_binary": binary or BINARY,
    }
    print(f"  L{level} {rpm:g} rpm  {p['cyc_tot']:6.1f} cyc  t_end={p['t_end']:8.3f}  "
          f"walltime={walltime}  ({p['min_per_cycle']:.2f} min/cyc, {p['confidence']})")
    print(f"       Kim: dtmix_0.50={p['kim_dtmix_0.50']:.1f}s  "
          f"dtmix_0.95={p['kim_dtmix_0.95']:.1f}s")
    if dry:
        print(f"       DRY RUN -- not submitted ({run_id})")
        return
    job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                       walltime=walltime,
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=ntasks, mem="4G")
    print(f"       submitted {run_id}  job={job}")


# The Fig 8 late-time probe left a CONVERGED L10 state at exactly Fig 9's
# reference condition (7 deg, 32.5 rpm), at absolute cycle 47.0. Reusing it
# skips 47 cycles of spin-up already paid for.
L10_SEED_RUN = "57f68830"
L10_SEED_CYCLE = 47.0


def submit_l10_warm(dry: bool) -> None:
    """The single most useful L10 run available: Fig 9 at Kim's own mesh.

    Fig 8 already has L10 and agrees with Kim (0.93-0.97x). Fig 9 has NO L10
    data at all, and Kim's production mesh IS n_L=2^10 -- so this is the only
    run that can test whether our mixing times converge to his, rather than
    comparing a coarse grid against his fine one.

    Warm-started from runs/57f68830's checkpoint (cycle 47.0, same condition,
    depth 10) with n_mix_cycles = 80 - 47 = 33, so the tracer is released at
    ABSOLUTE cycle 80 -- Kim's exact protocol. The saving is therefore free of
    protocol risk: it does not depend on whether the separate spin-up-25 test
    (job 6417004) shows Kim's 80 cycles to be unnecessary.

    restart_continue stays 0 (the default): this IS a new tracer experiment,
    and 57f68830 never released a tracer anyway (its t_end < t_mix), so its
    stracers are zero already.
    """
    rpm = 32.5
    kim = pd.read_csv(KIM_CSV).set_index("RPM")
    T_per, T_bio = t_scales(rpm)
    n_mix = SPINUP_CYCLES - L10_SEED_CYCLE          # 33.0
    cyc_mix = MARGIN * float(kim.loc[rpm, "dtmix_strict_0.95"]) / T_per
    cyc_seg = n_mix + cyc_mix
    t_end = cyc_seg * T_per / T_bio                 # relative to the checkpoint
    seed_ck = ROOT / "runs" / L10_SEED_RUN / "checkpoint.dump"
    if not seed_ck.exists():
        raise SystemExit(f"seed checkpoint missing: {seed_ck}")

    params = {
        "run_id": f"fig9_l10_rpm{rpm:g}", "fidelity": 10,
        "geometry": GEOMETRY, "fill_level": 0.5,
        "n_harmonics": 1, "theta_max": THETA_MAX,
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": rpm * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": round(t_end, 4),
        "n_mix_cycles": int(n_mix),
        "restart_continue": 0,
        "_binary": LEAN_BINARY,
    }
    print(f"L10 {rpm:g} rpm, warm from {L10_SEED_RUN} @ cycle {L10_SEED_CYCLE}")
    print(f"  tracer released at absolute cycle "
          f"{L10_SEED_CYCLE + n_mix:.0f} (Kim's protocol: 80)")
    print(f"  segment = {cyc_seg:.1f} cycles (cold would be "
          f"{SPINUP_CYCLES + cyc_mix:.1f}), t_end={t_end:.3f}")
    print(f"  Kim @32.5rpm: dtmix 0.50/0.75/0.95 = "
          f"{kim.loc[rpm,'dtmix_strict_0.5']:.1f} / "
          f"{kim.loc[rpm,'dtmix_strict_0.75']:.1f} / "
          f"{kim.loc[rpm,'dtmix_strict_0.95']:.1f} s")
    if dry:
        print("  DRY RUN -- not submitted")
        return
    job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                       walltime="48:00:00",
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=32, mem="4G", checkpoint=str(seed_ck))
    print(f"  submitted {params['run_id']}  job={job}  (32 ranks, 48 h)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage",
                    choices=["validate", "ladder", "lean", "l10", "streamcheck", "kimcheck",
                             "sweep"],
                    required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    if a.stage == "validate":
        print("Stage 'validate': one L6 point, full Kim protocol.")
        submit(32.5, level=6, ntasks=8, prefix="fig9_validate", dry=a.dry_run)
    elif a.stage == "ladder":
        # Kim's production mesh is n_L = 2^10 (Main.tex sec. 4), i.e. our
        # L10 -- verified by his own stated 0.24 mm cell size and 1.05e6
        # cell count, which our L10 reproduces exactly. Coarser levels
        # sitting below his curve are therefore a resolution difference, not
        # a discrepancy. This ladder measures the convergence RATE toward his
        # mesh so the sweep level is chosen from data: L6 12.9x, L7 6.5x on
        # dtmix_0.95 (gain 1.98/level, extrapolating to parity at ~L10).
        # 32.5 rpm, the condition with the most reference data.
        #
        # Walltimes are scaled from the MEASURED L6 run (209.7 cycles in
        # 8m36s on 8 ranks, no video) at 4x per level, not from cost_model,
        # whose (6,8) entry is 63x pessimistic here because it was measured
        # on a video-writing binary. Generous margins since 4x/level is
        # itself an assumption.
        # L9 runs on the LEAN binary. On the baseline it needs 24-37 h
        # (6.5-8x/level from L7's 34m17s) against a 30 h cap, and a timeout
        # loses everything because results are only copied back on completion.
        # Lean is 1.69x faster -> 14-22 h. L7 exists on both binaries with
        # identical dtmix, so the ladder stays internally consistent.
        for level, walltime, lean in [(7, "02:00:00", False),
                                      (8, "08:00:00", False),
                                      (9, "30:00:00", True)]:
            submit(32.5, level=level, ntasks=8, prefix="fig9_ladder",
                   dry=a.dry_run, walltime_override=walltime,
                   binary=LEAN_BINARY if lean else None)
    elif a.stage == "lean":
        # A/B against fig9_ladder_l7 (34m17s): identical config, identical
        # walltime cap, only EXTRA_TRACERS=0. c/c1/c3 are never initialised
        # and never read, but each costs 2 VOF-advected fields plus a
        # multigrid diffusion solve every timestep.
        submit(32.5, level=7, ntasks=8, prefix="fig9_lean", dry=a.dry_run,
               walltime_override="02:00:00", binary=LEAN_BINARY)
    elif a.stage == "l10":
        submit_l10_warm(dry=a.dry_run)
    elif a.stage == "kimcheck":
        # Short L8 on the Kim-exact binary, to size the mask difference BEFORE
        # committing to reruns of the Fig 8/13 sweeps. Kim masks on cs==1 AND
        # f==1 and uses SIGNED tau; we used cs>0, f>0.5, |tau|.
        rpm, T_per, T_bio = 32.5, *t_scales(32.5)
        params = {
            "run_id": "kimcheck_l8_rpm32.5", "fidelity": 8,
            "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
            "theta_max": THETA_MAX, "phi_angular": [0.0, 0.0, 0.0],
            "omega_b": rpm * 2 * math.pi / 60.0, "omega_h": 0.0,
            "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
            "t_end": round(15 * T_per / T_bio, 4), "n_mix_cycles": 80,
            "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-kimexact",
        }
        print(f"kimcheck: L8 15 cycles, t_end={params['t_end']}")
        if not a.dry_run:
            job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                               walltime="02:00:00",
                               template=ROOT / "config" / "slurm_mpi_template.sh",
                               cpus=1, ntasks=8, mem="4G")
            print(f"  submitted job={job}")
    elif a.stage == "streamcheck":
        # Cheap L6 run, just long enough for steady streaming to establish
        # (settled by ~cycle 23 at L6), to check the new streaming.dat output.
        # Target: Kim reports <|xi_bar|> = 1.3295 1/s at 32.5 rpm, i.e. 4.04 in
        # non-dim units (x T_bio = 3.039 s). The OLD, wrong quantity sat at
        # ~6.2 non-dim. A correct streaming vorticity must be far below that.
        rpm, T_per, T_bio = 32.5, *t_scales(32.5)
        params = {
            "run_id": "streamcheck_l6_rpm32.5", "fidelity": 6,
            "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
            "theta_max": THETA_MAX, "phi_angular": [0.0, 0.0, 0.0],
            "omega_b": rpm * 2 * math.pi / 60.0, "omega_h": 0.0,
            "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
            "t_end": round(40 * T_per / T_bio, 4), "n_mix_cycles": 80,
            "_binary": "/oscar/scratch/eaguerov/BioReactor-mpi-streaming",
        }
        print(f"streamcheck: L6 40 cycles, t_end={params['t_end']}")
        if a.dry_run:
            print("  DRY RUN")
        else:
            job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                               walltime="01:00:00",
                               template=ROOT / "config" / "slurm_mpi_template.sh",
                               cpus=1, ntasks=8, mem="4G")
            print(f"  submitted job={job}")
    else:
        # L9, 32 ranks, all 10 of Kim's rpm points, run in PARALLEL.
        #
        # Not L10: measured L10 cost is 58.2 min/cycle at 32 ranks, so the same
        # sweep there is ~2200 h and the longest point alone exceeds the 48 h
        # cap. L9 is ~11x cheaper (5.2 min/cycle) -- every point fits in one
        # job, longest 26 h. The L10 anchor at 32.5 rpm is run separately, in
        # segments, to tie the sweep to Kim's own mesh.
        print("Stage 'sweep': L9, 32 ranks, all 10 rpm points (parallel).")
        for rpm in [37.5, 35, 32.5, 30, 27.5, 25, 22.5, 20, 17.5, 15]:
            submit(float(rpm), level=9, ntasks=32, prefix="fig9",
                   dry=a.dry_run, binary=PROD_BINARY,
                   walltime_override=_l9_walltime(float(rpm)))


def _l9_walltime(rpm: float) -> str:
    """Per-point walltime from the MEASURED L10 cost divided by the measured
    per-level factor, with 1.6x margin. Never a guess -- see cost_model."""
    from scripts.cost_model import min_per_cycle, PER_LEVEL_WORK_FACTOR
    mc10, _ = min_per_cycle(10, 32, 32.5)
    mc9 = mc10 / PER_LEVEL_WORK_FACTOR
    kim = pd.read_csv(KIM_CSV).set_index("RPM")
    T_per, _ = t_scales(rpm)
    cyc = SPINUP_CYCLES + MARGIN * float(kim.loc[rpm, "dtmix_strict_0.95"]) / T_per
    hours = cyc * mc9 / 60.0 * 1.6
    if hours > 47:
        raise SystemExit(f"L9 {rpm} rpm needs {hours:.0f} h > 48 h cap")
    return f"{int(hours) + 1:02d}:00:00"


if __name__ == "__main__":
    main()
