"""L7 runs carrying tracer + oxygen, for the H2 test (L7 -> L9 on Kim's numbers).

Target quantities are Kim's: kLa at C* = 10/25/50% and dtmix at chi =
0.50/0.75/0.95. H2: an L7 sweep (~64x cheaper per cycle than L9) plus 3 L9
points predicts the L9 values better than an L9-only GP and an L7 x ratio scaling.

Warm starts (same condition): every rpm/angle except 15 rpm restarts from its
settled 80-cycle L7 checkpoint (mf_l7_rpm*, mf_l7_th*), with all *_prev equal to
the current condition, so the restart ramp is a no-op. Tracer and oxygen
release N_SETTLE cycles after the restart (restart_continue = 0: a fresh
mixing/oxygen experiment). 15 rpm has no L7 checkpoint and cold-starts with
Kim's 80-cycle spin-up.
POST cycles after release: L9 needed <= 126 cycles to reach C* = 50%, and L7
transfers and mixes faster (coarser grid). 150 is a margin, at ~0.1-0.2 min/cycle.

Binary: f1c11e0, the one every L7/L9 run in this comparison used. Source is
newer only by the period-rounding fix (133e55e), which moves the end time by at
most one period. Consistency across levels matters more here.

Usage:  uv run python scripts/submit_l7_kla_mix.py [--dry-run]
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
GEOMETRY = {"a": 0.25, "b": 0.03575, "n": 8.0}
N_SETTLE, POST, SPINUP, NTASKS = 1, 150, 80, 8
POINTS = ([(r, 7.0, f"mf_l7_rpm{r:g}") for r in (17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5)]
          + [(15.0, 7.0, None)]
          + [(32.5, th, f"mf_l7_th{th:g}") for th in (2.0, 3.0, 4.0, 5.0, 6.0)])


def scales(rpm, th):
    omega_b = rpm * 2 * math.pi / 60.0
    L, H = GEOMETRY["a"], 2 * GEOMETRY["b"]
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(th)))
    return omega_b, T_per / (L / (V / (H * 0.5) / T_per))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    zeros = [0.0, 0.0, 0.0]
    for rpm, th, seed in POINTS:
        omega_b, T_nd = scales(rpm, th)
        theta = [th, 0.0, 0.0]
        tag = f"rpm{rpm:g}" if th == 7.0 else f"th{th:g}"
        p = {"run_id": f"kmix_l7_{tag}", "fidelity": 7, "geometry": GEOMETRY,
             "fill_level": 0.5, "n_harmonics": 1, "theta_max": theta,
             "phi_angular": zeros, "omega_b": omega_b, "omega_h": 0.0,
             "amplitude_h": zeros, "phi_horizontal": zeros, "_binary": BINARY}
        ckpt = None
        if seed:
            ckpt = ROOT / "runs" / seed / "checkpoint.dump"
            t_ck = _dump_fields(ckpt)[0]["t"]
            if t_ck <= 0:
                raise SystemExit(f"{ckpt}: t={t_ck}, cannot arm a restart")
            p.update(t_checkpoint=t_ck, omega_b_prev=omega_b, theta_max_prev=theta,
                     phi_angular_prev=zeros, amplitude_h_prev=zeros,
                     phi_horizontal_prev=zeros, omega_h_prev=0.0,
                     n_mix_cycles=N_SETTLE,
                     t_end=round((N_SETTLE + POST) * T_nd, 4))
            how = f"warm from {seed} (t={t_ck:.2f})"
        else:
            p.update(n_mix_cycles=SPINUP, t_end=round((SPINUP + POST) * T_nd, 4))
            how = "cold start"
        print(f"  {p['run_id']:18s} {how}, t_end={p['t_end']}")
        if a.dry_run:
            continue
        job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs",
                           walltime="03:00:00",
                           template=ROOT / "config" / "slurm_mpi_template.sh",
                           cpus=1, ntasks=NTASKS, mem="4G",
                           **({"checkpoint": str(ckpt)} if ckpt else {}))
        print(f"    job={job}")


if __name__ == "__main__":
    main()
