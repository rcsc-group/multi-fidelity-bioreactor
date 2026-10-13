"""Recorded tracer run for test V at L10 (and an L8 smoke test of the same pipeline).

Restart from a settled run's checkpoint (same level), release the tracer after N_MIX cycles, write uniform-grid
snapshots (13/period) from release to the end. Purpose: D_eff/D_phys at L10 for the test V figure (diary 2026-10-12).
Source checkpoints: x8_l10_rpm32.5_th7 (L10, settled: x8 passed at 32.5 rpm) or kmix_l8_rpm32.5 (L8 smoke).
*_prev = current (no forcing ramp effect); restart_continue unset, so the tracer is released at t_ck + N_MIX periods.

Usage: uv run python scripts/submit_v_record.py --source x8_l10_rpm32.5_th7 --level 10 --cycles 21 [--dry-run]
       uv run python scripts/submit_v_record.py --source kmix_l8_rpm32.5 --level 8 --cycles 2 --smoke
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
N_MIX, FRAMES = 1, 13
MIN_PER_CYCLE = {6: 1.0, 7: 1.0, 8: 2.0, 10: 51.0}  # L6/L7 generous guesses (< L8); measured: x8_l10_rpm32.5 21157 s for 6.9 cycles at 32 ranks; L8 generous


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--cycles", type=float, required=True, help="total cycles after the checkpoint")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    sd = ROOT / "runs" / a.source
    p = json.loads((sd / "params.json").read_text())
    t_ck = float(np.loadtxt(sd / "shear_stress.dat", skiprows=1)[-1, 1])
    _, T = _t_scales(p)
    c_ck = t_ck / T
    run_id = f"{'vsmoke' if a.smoke else 'v'}_l{a.level}_rpm32.5_th7"
    zeros = [0.0, 0.0, 0.0]
    p.update(run_id=run_id, fidelity=a.level, t_checkpoint=t_ck, t_end=round(a.cycles * T, 6), n_mix_cycles=N_MIX,
             omega_b_prev=p["omega_b"], theta_max_prev=p["theta_max"], phi_angular_prev=zeros,
             amplitude_h_prev=zeros, phi_horizontal_prev=zeros, omega_h_prev=0.0, _binary=BINARY,
             frames_per_period=FRAMES, snap_start_cycle=round(c_ck + N_MIX - 0.01, 4),
             snap_end_cycle=round(c_ck + a.cycles + 1, 4))
    p.pop("restart_continue", None)
    ntasks = 16 if a.level <= 8 else 32
    hours = max(1, math.ceil(a.cycles * MIN_PER_CYCLE[a.level] * 1.3 / 60.0))
    print(f"{run_id}: source {a.source} t_ck {t_ck:g} (cycle {c_ck:.3f}), release at cycle {c_ck + N_MIX:.3f}, "
          f"t_end +{p['t_end']}, snapshots cycles {p['snap_start_cycle']}-{p['snap_end_cycle']}, "
          f"{ntasks} ranks, walltime {hours} h")
    if a.dry_run:
        return
    job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=f"{hours:02d}:00:00",
                       template=ROOT / "config" / "slurm_mpi_template.sh",
                       checkpoint=str((sd / "checkpoint.dump").resolve()), cpus=1, ntasks=ntasks, mem="4G")
    print(f"  job={job}")


if __name__ == "__main__":
    main()
