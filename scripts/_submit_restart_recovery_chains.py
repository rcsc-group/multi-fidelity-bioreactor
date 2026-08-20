"""One-off submission helper for the checkpoint-restart recovery test
(diary.md 2026-08-20). Reuses chain.py's build_chain() (correct fresh-vs-
restart params.json convention) and simulate.submit_slurm() directly,
skipping postprocess.validate_params(): Kim et al.'s own published
geometry (b=0.03575) sits just outside config/param_space.yaml's [0.05,
0.15] sweep-exploration bound for geometry.b, which is a deliberate
design-space choice for the optimization problem, not a numerical-
validity guard -- this experiment is a validation anchor against the
established Kim case, not a sweep candidate, so the bound doesn't apply.
Not a change to any production file's validation behavior.

Usage:
    uv run python scripts/_submit_restart_recovery_chains.py config/chain_restart_recovery_test.yaml
    uv run python scripts/_submit_restart_recovery_chains.py config/chain_restart_recovery_baseline.yaml
"""
import sys
from pathlib import Path

import yaml

_PROJECT_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(_PROJECT_ROOT))

import scripts.simulate as simulate
from scripts.chain import build_chain, _DEFAULT_TEMPLATE

cfg = yaml.safe_load(Path(sys.argv[1]).read_text())
chain = build_chain(cfg)
runs_root = _PROJECT_ROOT / "runs"
walltime = cfg.get("walltime", "02:00:00")

job_ids = []
prev_run_id = None
for k, params in enumerate(chain):
    checkpoint = None
    if prev_run_id is not None:
        checkpoint = str((runs_root / prev_run_id / "checkpoint.dump").resolve())
    dependency = f"afterok:{job_ids[k-1]}" if k > 0 else None

    field_val = params.get("theta_max", ["?"])[0]
    print(f"[seg {k}] run={params['run_id']}  theta_max_0={field_val}  "
          f"n_mix={params['n_mix_cycles']}  t_end~{params['t_end']:.1f}  "
          f"theta_max_prev={'omitted (fresh)' if 'theta_max_prev' not in params else params['theta_max_prev']}")

    job_id = simulate.submit_slurm(
        params,
        project_root=_PROJECT_ROOT,
        runs_root=runs_root,
        walltime=walltime,
        template=_DEFAULT_TEMPLATE,
        checkpoint=checkpoint,
        dependency=dependency,
        cpus=4,
    )
    job_ids.append(job_id)
    print(f"  -> job {job_id}")
    prev_run_id = params["run_id"]

print("job_ids:", job_ids)
