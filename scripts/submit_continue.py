"""Continue a finished run at the SAME condition for N more rocking cycles (same-condition continuation).

This is "chaining" in the project's vocabulary, not "checkpointing": the experiment (tracer + oxygen) goes on.
  - params = the source run's params.json (same level, rpm, theta, geometry, binary `_binary`);
  - t_checkpoint = last logged time of the source (shear_stress.dat, as scripts/submit_x8_extend.py);
    the solver re-reads the exact dump time from checkpoint.dump;
  - restart_continue = 1: tracers and oxygen keep their restored values and are NOT re-released;
  - *_prev = current forcing: no ramp from zero (the solver ramps *_prev -> current over N_RAMP_CYCLES);
  - n_mix_cycles = 0 by default: with restart_continue the oxygen event starts at t_mix = t_checkpoint +
    n_mix_cycles*T, so 0 keeps the gas reservoir reset active from the first step, as in an uninterrupted run;
  - t_end = (cycles - 0.005)*T_per_nd after the checkpoint (relative; the solver stops at the next period
    boundary, so an integer number of cycles is reached without overshooting to cycles+1);
  - _parent_run = source run_id, so scripts/postprocess.py joins the segments (_segment_chain/_joined_cols).
checkpoint.phase is copied next to the staged checkpoint (the solver reads it from its working directory).

Usage: uv run python scripts/submit_continue.py --source kmix_l8_rpm25 --tag cont --cycles 120 [--binary B]
       [--ntasks 16] [--walltime 08:00:00] [--n-mix 0] [--dry-run]
New run id: <source>_<tag>.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _t_scales  # noqa: E402
from scripts.simulate import submit_slurm  # noqa: E402

MPI_SCRATCH = Path("/oscar/scratch/eaguerov/mpi_runs")


def build_params(src_params: dict, src_id: str, tag: str, t_ck: float, cycles: float, n_mix: float,
                 binary: str | None) -> dict:
    p = dict(src_params)
    for k in ("_canonical_run_dir", "next_run_id"):
        p.pop(k, None)
    _, T_nd = _t_scales(p)
    p.update(run_id=f"{src_id}_{tag}", t_checkpoint=t_ck, t_end=round((cycles - 0.005) * T_nd, 6),
             n_mix_cycles=n_mix, restart_continue=1, _parent_run=src_id,
             omega_b_prev=p["omega_b"], theta_max_prev=list(p["theta_max"]),
             phi_angular_prev=list(p["phi_angular"]), amplitude_h_prev=list(p["amplitude_h"]),
             phi_horizontal_prev=list(p["phi_horizontal"]), omega_h_prev=p["omega_h"])
    if binary:
        p["_binary"] = binary
    return p


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, help="source run id (runs/<source>/ with checkpoint.dump)")
    ap.add_argument("--tag", required=True, help="new run id = <source>_<tag>")
    ap.add_argument("--cycles", type=float, required=True, help="rocking cycles to run after the checkpoint")
    ap.add_argument("--binary", default=None, help="default: the source's params _binary")
    ap.add_argument("--ntasks", type=int, default=16)
    ap.add_argument("--walltime", default="08:00:00")
    ap.add_argument("--n-mix", type=float, default=0.0, help="n_mix_cycles of the continuation (oxygen event time)")
    ap.add_argument("--exclude", default="node2356")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    src = ROOT / "runs" / a.source
    ck = src / "checkpoint.dump"
    if not ck.exists():
        sys.exit(f"{ck} not found")
    sp = json.loads((src / "params.json").read_text())
    t_ck = float(np.loadtxt(src / "shear_stress.dat", skiprows=1)[-1, 1])
    p = build_params(sp, a.source, a.tag, t_ck, a.cycles, a.n_mix, a.binary)
    _, T_nd = _t_scales(p)
    print(f"  {p['run_id']}: source {a.source}, t_checkpoint {t_ck:.6f} ({t_ck / T_nd:.3f} cycles), "
          f"t_end(rel) {p['t_end']} ({a.cycles:g} cycles), binary {p['_binary']}, {a.ntasks} ranks, {a.walltime}")
    if a.dry_run:
        return
    stage = MPI_SCRATCH / p["run_id"]
    stage.mkdir(parents=True, exist_ok=True)
    if (src / "checkpoint.phase").exists():
        shutil.copy2(src / "checkpoint.phase", stage / "checkpoint.phase")
    job = submit_slurm(p, project_root=ROOT, runs_root=ROOT / "runs", walltime=a.walltime,
                       template=ROOT / "config" / "slurm_mpi_template.sh", checkpoint=str(ck.resolve()),
                       cpus=1, ntasks=a.ntasks, mem="4G", exclude=a.exclude or None)
    print(f"    job={job}")


if __name__ == "__main__":
    main()
