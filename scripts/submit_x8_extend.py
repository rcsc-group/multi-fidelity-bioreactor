"""Continue an x8_l10 warm-start run at L10 for N_CYCLES more cycles (same-level chaining, L10 -> L10).

Why (diary 2026-10-11): the pre-registered L8 -> L10 warm-start test FAILED at 17.5 rpm on tau_mean (+7.3% vs
l10_fig13a_xlevel, tol 5%) while tau_mean was still falling at cycle 6. The continuation decides whether the gap is
a slow transient (it closes) or a persistent bias of the L8 start (it does not).

Source: runs/x8_l10_rpm<r>_th7/checkpoint.dump (written at its t_end). Same binary and forcing; *_prev = current.
Usage: uv run python scripts/submit_x8_extend.py --rpm 17.5 25 [--dry-run]
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
from scripts.simulate import submit_slurm  # noqa: E402

N_CYCLES, NTASKS, MARGIN = 6, 32, 1.35
_HOURS_PER_CYCLE = {17.5: 19.36 / 6, 25: 15.2 / 5.6}  # measured from x8_l10 wall clock (logstats.dat)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rpm", nargs="+", type=float, required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for rpm in a.rpm:
        src = ROOT / "runs" / f"x8_l10_rpm{rpm:g}_th7"
        p = json.loads((src / "params.json").read_text())
        t_dump = float(np.loadtxt(src / "shear_stress.dat", skiprows=1)[-1, 1])
        p.update(run_id=f"x8e_l10_rpm{rpm:g}_th7", t_checkpoint=t_dump,
                 t_end=round(p["t_end"] / p["n_mix_cycles"] * N_CYCLES, 6), n_mix_cycles=N_CYCLES,
                 omega_b_prev=p["omega_b"], theta_max_prev=p["theta_max"], phi_angular_prev=p["phi_angular"],
                 amplitude_h_prev=p["amplitude_h"], phi_horizontal_prev=p["phi_horizontal"],
                 omega_h_prev=p["omega_h"])
        hours = max(4, math.ceil(N_CYCLES * _HOURS_PER_CYCLE[rpm] * MARGIN))
        print(f"  {p['run_id']}: source {src.name} (t {t_dump:g}), t_end {p['t_end']}, walltime {hours} h")
        if a.dry_run:
            continue
        job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=f"{hours:02d}:00:00",
                           template=ROOT / "config" / "slurm_mpi_template.sh",
                           checkpoint=str((src / "checkpoint.dump").resolve()), cpus=1, ntasks=NTASKS, mem="4G")
        print(f"    job={job}")


if __name__ == "__main__":
    main()
