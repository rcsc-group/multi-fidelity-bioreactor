"""The L10 anchor for Kim Fig. 9: one point at 32.5 rpm on Kim's own mesh.

Why a separate script. The Fig 9 SWEEP runs at L9 (10 rpm points, parallel,
longest 26 h). L10 is not viable for a sweep -- measured 58.2 min/cycle at 32
ranks means ~2200 h for ten points -- but one L10 point is essential, because
Kim's production mesh IS n_L=2^10 (Main.tex sec. 4) and without it every
comparison is our coarse grid against his fine one. This point ties the L9
sweep to his mesh.

Why segments. 162.7 cycles x 58.2 min = 158 h against a 48 h cap. Segments of
~41 cycles (~40 h) chain through it. This is the first production use of
`restart_continue=1` (010fa71): without it every segment after the first
would re-inject the tracer at t_checkpoint + n_mix_cycles and restart chi
from zero mid-experiment.

Protocol fidelity. Segment 1 warm-starts from runs/57f68830 (a converged L10
state at this exact condition, absolute cycle 47.0) with n_mix_cycles = 33,
so the tracer releases at ABSOLUTE cycle 80 -- Kim's protocol exactly. That
saving reuses 47 cycles already paid for; it does NOT shorten the protocol,
which matters because job 6417004 measured an 18% shift in dtmix_0.95 when
release moved from cycle 80 to 25 (diary.md 2026-09-15 (11)).

Usage:
    uv run python scripts/submit_fig9_l10_anchor.py --segment 1
    uv run python scripts/submit_fig9_l10_anchor.py --segment 2 --from-run fig9_l10_seg1
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.simulate import submit_slurm            # noqa: E402
from scripts.dump_fields import fields as _dump_fields  # noqa: E402
from scripts.cost_model import min_per_cycle         # noqa: E402

BINARY = "/oscar/scratch/eaguerov/BioReactor-mpi-prod-4b3a435"
SEED_RUN, SEED_CYCLE = "57f68830", 47.0
RPM, SPINUP, MARGIN = 32.5, 80, 1.30
NTASKS = 32
SEG_HOURS = 40.0            # per segment, against the 48 h Exploratory cap
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
THETA = [7.0, 0.0, 0.0]


def t_scales():
    omega_b = RPM * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(THETA[0])))
    return T_per, L / (V / (H * 0.5) / T_per)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--segment", type=int, required=True)
    ap.add_argument("--from-run", help="previous segment's run_id (segments >1)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    T_per, T_bio = t_scales()
    mc, conf = min_per_cycle(10, NTASKS, RPM)
    cyc_per_seg = SEG_HOURS * 60.0 / mc

    kim = pd.read_csv(ROOT / "experiments/kimetal2024/csv_raw/"
                      "mixing_kla_vs_frequency.csv").set_index("RPM")
    total_cyc = (SPINUP - SEED_CYCLE) + MARGIN * float(
        kim.loc[RPM, "dtmix_strict_0.95"]) / T_per

    if a.segment == 1:
        ckpt = ROOT / "runs" / SEED_RUN / "checkpoint.dump"
        n_mix, cont = SPINUP - SEED_CYCLE, 0      # release at absolute cycle 80
        parent = None                             # the seed is a different experiment
    else:
        if not a.from_run:
            raise SystemExit("--from-run required for segments > 1")
        ckpt = ROOT / "runs" / a.from_run / "checkpoint.dump"
        n_mix, cont = 1, 1                        # continuation: no re-injection
        parent = a.from_run
    if not ckpt.exists():
        raise SystemExit(f"checkpoint not found: {ckpt}")

    seg_cyc = min(cyc_per_seg, total_cyc - (a.segment - 1) * cyc_per_seg)
    if seg_cyc <= 0:
        raise SystemExit(f"segment {a.segment} is past the end ({total_cyc:.1f} cyc total)")

    # The restart branch is gated on params.t_checkpoint > 0 (BioReactor.c:~437),
    # NOT on a dump being staged. submit_slurm's checkpoint= only copies the file
    # into place; without t_checkpoint the binary IGNORES argv[2] and silently
    # runs a COLD start. That is exactly what happened to kimcheck_l10_rpm32.5,
    # which was reported as warm-started from 57f68830 and was not. The C code
    # then overwrites this with the dump's own time, so the value here only has
    # to be the true one for the pre-restore timing to be right.
    t_ck = _dump_fields(ckpt)[0]["t"]
    if t_ck <= 0:
        raise SystemExit(f"checkpoint {ckpt} reports t={t_ck}; cannot arm a restart")

    params = {
        "run_id": f"fig9_l10_seg{a.segment}", "fidelity": 10,
        "t_checkpoint": t_ck,
        "geometry": GEOMETRY, "fill_level": 0.5, "n_harmonics": 1,
        "theta_max": THETA, "phi_angular": [0.0, 0.0, 0.0],
        "omega_b": RPM * 2 * math.pi / 60.0, "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0], "phi_horizontal": [0.0, 0.0, 0.0],
        "t_end": round(seg_cyc * T_per / T_bio, 4),
        "n_mix_cycles": int(n_mix),
        "restart_continue": cont,
        "_binary": BINARY,
    }
    # postprocess walks _parent_run backwards and joins the raw series before
    # computing dtmix or kLa. Without it a segment is scored in isolation, which
    # for a continuation is meaningless: oxygen starts saturated so every kLa
    # threshold is crossed at the first row, and the tracer starts partly mixed
    # so chi never starts at 0. Only continuations set it -- segment 1 warm-starts
    # from a DIFFERENT experiment and must not be joined to it.
    if parent:
        params["_parent_run"] = parent
    print(f"L10 anchor segment {a.segment}: {seg_cyc:.1f} cycles "
          f"(of {total_cyc:.1f} total), t_end={params['t_end']}")
    print(f"  from {ckpt.parent.name}, restart_continue={cont}, "
          f"n_mix_cycles={n_mix}  [{mc:.1f} min/cyc, {conf}]")
    if a.segment == 1:
        print(f"  tracer releases at ABSOLUTE cycle {SEED_CYCLE + n_mix:.0f} "
              f"(Kim's protocol: {SPINUP})")
    if a.dry_run:
        print("  DRY RUN")
        return
    job = submit_slurm(params, project_root=ROOT, runs_root=ROOT / "runs",
                       walltime=f"{int(SEG_HOURS) + 4:02d}:00:00",
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       cpus=1, ntasks=NTASKS, mem="4G", checkpoint=str(ckpt))
    print(f"  submitted {params['run_id']}  job={job}")


if __name__ == "__main__":
    main()
