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
           walltime_override: str | None = None) -> None:
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
        "_binary": BINARY,
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["validate", "ladder", "sweep"], required=True)
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
        for level, walltime in [(7, "02:00:00"), (8, "08:00:00"), (9, "30:00:00")]:
            submit(32.5, level=level, ntasks=8, prefix="fig9_ladder",
                   dry=a.dry_run, walltime_override=walltime)
    else:
        print("Stage 'sweep': L7, all 10 of Kim's rpm points.")
        for rpm in [37.5, 35, 32.5, 30, 27.5, 25, 22.5, 20, 17.5, 15]:
            submit(float(rpm), level=7, ntasks=8, prefix="fig9", dry=a.dry_run)


if __name__ == "__main__":
    main()
