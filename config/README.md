# config/ — inputs for runs, sweeps and the optimiser

| File | Used by | What it sets |
|---|---|---|
| `param_space.yaml` | sampling, optimiser | The control parameters (rocking rate, angle, harmonics, ...) and their bounds. |
| `sample_config.yaml` | `scripts/sample.py` | A batch of simulations to generate and submit. |
| `bo_config.yaml` | `scripts/loop.py` | The multi-fidelity optimisation loop. |
| `chain_config.yaml`, `chain_config_smoke.yaml` | `scripts/chain.py` | A sweep where each run restarts from the previous run's checkpoint. `_smoke` is a short check. |
| `chain_restart_recovery_*.yaml` | `scripts/chain.py` | A pair (baseline, test) that checks that a restart gives the same result as a continuous run. |
| `sweep_example.json` | `scripts/sweep.py` | Example of a one-parameter sweep (`_sweep` block). |
| `sweep_fb_*.json` | `scripts/sweep.py` | Angle (`theta`) and fill (`fill`) sweeps at level 7 and 8 behind early results in `experiments/`. |
| `slurm_template.sh` | `make submit` | SLURM script for one serial run. |
| `slurm_mpi_template.sh` | most launchers | SLURM script for one MPI run. |
| `slurm_video_template.sh`, `slurm_health_template.sh` | `scripts/run/submit_video_run.py`, `scripts/run/submit_health_runs.py` | SLURM scripts for the `-video` and `-health` binaries. |

One run's parameters are a JSON file with the keys in `sweep_example.json` (without `_sweep`). The
solver reads it with `src/params_read.h`; the meaning of each key is in
[the parameter reference](../docs_site/reference/params.md) of the documentation.

Files that no script used any more were removed; see `docs/history/archived-files.md`.
