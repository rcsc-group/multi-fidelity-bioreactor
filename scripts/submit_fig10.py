"""Kim Fig. 10/12: mixing time and kLa vs rocking ANGLE, at 32.5 rpm.

One sweep serves both figures. Kim's own numbers for them live in a single
block of dataset_kimetal2024.xlsx labelled "Fig. 10,12" -- mixing time on one
axis, kLa on the other, same six runs -- exactly as "Fig. 9,11" serves the rpm
pair. postprocess already emits both KPIs from one run, so Fig 12 costs
nothing beyond Fig 10.

Cheaper than Fig 9 in one respect: Kim publishes dtmix_strict_0.5 for the
angle series and nothing tighter, so the target is chi = 0.50 rather than
0.95. That matters because the tail is where our grid disagrees with his
(scripts/autoextend.py) -- at chi = 0.50 the 37.5 rpm point came in at 0.71x
his time, well inside any sane margin.

EACH ANGLE COLD-STARTS. Warm-starting theta=6 from theta=7's converged state
looks like it would save five 80-cycle spin-ups, and does not: the measured
theta transitions at this rpm take LONGER to settle than a spin-up costs --
121 cycles at L6 and 91 at L8 for 7 deg -> 4 deg (experiments/
settling_model_fit.json) -- and every measured 4 deg -> 2 deg transition is a
flagged anomaly that never settled to its reference at any level. Chaining
here would buy a longer, less faithful run. The saving is real for a repeated
CONDITION, not for a changed one.

The low-angle points do not fit one job. theta=2 needs ~337 cycles against a
48 h cap, so t_end is truncated to what the cap holds and the point is
finished by its own continuation segments -- the same mechanism that finishes
Fig 9's short points, not a hand-stitched chain.

Usage:
    uv run python scripts/submit_fig10.py --dry-run
    uv run python scripts/submit_fig10.py --angles 7 6
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.autoextend import walltime_safety          # noqa: E402
from scripts.cost_model import min_per_cycle            # noqa: E402
from scripts.simulate import submit_slurm               # noqa: E402

KIM_CSV = ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_angle.csv"
# Lean build (one tracer) at 3645d13. Newer than prod-4b3a435 in one way that
# matters here: it writes posY_left/posY_right, so the angle sweep also feeds
# Figs. 14(b) and 15(d-f) instead of needing its own re-run.
BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-lean-3645d13"

RPM = 32.5
LEVEL, NTASKS = 9, 32
ANGLES = [7.0, 6.0, 5.0, 4.0, 3.0, 2.0]
SPINUP_CYCLES = 80      # Kim's tracer release instant, t/T_p = 80
MARGIN = 1.30
WALLTIME_CAP_H = 47     # against the 48 h Exploratory cap
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}


def t_scales(theta_deg: float) -> tuple[float, float]:
    """(T_per [s], T_bio [s]) -- same definitions as postprocess._t_scales.

    T_bio depends on theta through the swept volume, so it must be recomputed
    per angle; reusing the 7 deg value would mis-size every other point's
    t_end by the ratio of their tangents.
    """
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta_deg)))
    return T_per, L / (V / (H * 0.5) / T_per)


def plan(theta: float) -> dict:
    kim = pd.read_csv(KIM_CSV).set_index("Angle_deg")
    T_per, T_bio = t_scales(theta)
    dtmix_50 = float(kim.loc[theta, "dtmix_strict_0.5"])
    want = SPINUP_CYCLES + MARGIN * dtmix_50 / T_per
    mpc, conf = min_per_cycle(level=LEVEL, ntasks=NTASKS, rpm=RPM)
    safety = walltime_safety(conf)
    affordable = WALLTIME_CAP_H * 60.0 / (mpc * safety)
    cycles = min(want, affordable)
    return {"theta": theta, "cycles": cycles, "want": want,
            "truncated": cycles < want - 1e-6,
            "t_end": cycles * T_per / T_bio,
            "hours": cycles * mpc / 60.0 * safety,
            "kim_dtmix_0.50": dtmix_50, "min_per_cycle": mpc,
            "confidence": conf}


def submit(theta: float, dry: bool) -> None:
    p = plan(theta)
    T_per, T_bio = t_scales(theta)
    run_id = f"fig10_l{LEVEL}_th{theta:g}"
    params = {
        "run_id": run_id, "fidelity": LEVEL,
        "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
        "theta_max": [theta, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": RPM * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": round(p["t_end"], 4),
        "n_mix_cycles": SPINUP_CYCLES,
        "_binary": BINARY,
    }
    walltime = f"{min(int(p['hours']) + 1, WALLTIME_CAP_H):02d}:00:00"
    note = (f"  TRUNCATED at the walltime cap ({p['want']:.0f} cyc wanted); "
            f"autoextend will finish it" if p["truncated"] else "")
    print(f"  th={theta:g} deg  {p['cycles']:6.1f} cyc  t_end={p['t_end']:8.3f}  "
          f"walltime={walltime}  ({p['min_per_cycle']:.2f} min/cyc, {p['confidence']})")
    print(f"       Kim: dtmix_0.50={p['kim_dtmix_0.50']:.1f} s{note}")
    if dry:
        print(f"       DRY RUN -- not submitted ({run_id})")
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
    ap.add_argument("--angles", type=float, nargs="*", default=ANGLES)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    print(f"Fig 10/12 angle sweep: L{LEVEL}, {NTASKS} ranks, {RPM:g} rpm")
    for theta in a.angles:
        submit(theta, a.dry_run)


if __name__ == "__main__":
    main()
