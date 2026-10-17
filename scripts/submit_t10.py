"""Test T10 (pre-registered, diary 2026-10-17): L10 kLa_1T_10 (and dtmix_0.50 where reached) at the x8 design points.

Restart from the settled x8_l10 checkpoint (same level), release the tracer and oxygen 1 cycle after the restart
(the v_l10_rpm32.5_th7 protocol, scripts/submit_v_record.py), run N cycles, no snapshots. Binary lean-eec1f15 as v_l10.
node2356 is excluded (it killed 15 jobs 1 s after start, 2026-10-16/17).
Usage: uv run python scripts/submit_t10.py --source x8_l10_rpm17.5_th7 [...] [--cycles 30] [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _t_scales  # noqa: E402
from scripts.simulate import submit_slurm  # noqa: E402

BINARY = "/oscar/data/dharri15/eaguerov/bin/BioReactor-mpi-lean-eec1f15"
N_MIX = 1
WALL_H = 47          # 0.8-1.3 h per L10 cycle at 32 ranks (v_l10: 17 h for 21 cycles at 32.5 rpm); the solver checkpoints before the limit
EXCLUDE = "node2356"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", nargs="+", required=True)
    ap.add_argument("--cycles", type=float, default=30.0, help="cycles after the checkpoint")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    zeros = [0.0, 0.0, 0.0]
    for src in a.source:
        sd = ROOT / "runs" / src
        p = json.loads((sd / "params.json").read_text())
        rpm = p["omega_b"] * 60 / (2 * math.pi)
        theta = p["theta_max"][0]
        t_ck = float(np.loadtxt(sd / "shear_stress.dat", skiprows=1)[-1, 1])
        _, T = _t_scales(p)
        run_id = f"t10_l10_rpm{rpm:g}_th{theta:g}"
        p.update(run_id=run_id, t_checkpoint=t_ck, t_end=round(a.cycles * T, 6), n_mix_cycles=N_MIX,
                 omega_b_prev=p["omega_b"], theta_max_prev=p["theta_max"], phi_angular_prev=zeros,
                 amplitude_h_prev=zeros, phi_horizontal_prev=zeros, omega_h_prev=0.0, _binary=BINARY,
                 frames_per_period=13, snap_start_cycle=9999.0, snap_end_cycle=9999.0)
        p.pop("restart_continue", None)
        print(f"{run_id}: source {src} t_ck {t_ck:g} (cycle {t_ck / T:.2f}), release at cycle {t_ck / T + N_MIX:.2f}, "
              f"+{a.cycles:g} cycles, 32 ranks, walltime {WALL_H} h")
        if a.dry_run:
            continue
        job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=f"{WALL_H}:00:00",
                           template=ROOT / "config" / "slurm_mpi_template.sh",
                           checkpoint=str((sd / "checkpoint.dump").resolve()), cpus=1, ntasks=32, mem="4G",
                           exclude=EXCLUDE)
        print(f"  job={job}")


if __name__ == "__main__":
    main()
