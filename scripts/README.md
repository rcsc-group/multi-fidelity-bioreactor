# scripts/ — Python and shell around the solver

Run every Python script from the repository root, so that `from scripts.X import ...` works:

```bash
uv run python scripts/simulate.py --help
```

Most scripts that run the solver submit SLURM jobs. Do not run them on a login node.

## Start here: the main pipeline

| Step | Script | What it does |
|---|---|---|
| 1. One run | `launch.py`, `simulate.py`, `run_key.py` | Make a run directory with `params.json` and a SLURM script; run it locally or submit it. A submission stops if the same run (same key) is already finished or queued. |
| 2. Many runs | `sample.py`, `sweep.py`, `chain.py` | A batch from `config/param_space.yaml`; a one-parameter sweep; a sweep where each run restarts from the previous checkpoint. |
| 3. Finish a run | `autoextend.py`, `stage_segment.py`, `periods.py` | Extend a run that hit its time limit; stage a restart segment; where a run stops (same rule as the solver). |
| 4. Results | `postprocess.py`, `hydro_dataset_v2.py`, `health_report.py` | Quantities of interest from one run (kLa, shear stress, mixing time, energy dissipation); the hydrodynamic data set; numerical health checks. |
| 5. Surrogate and optimiser | `train_surrogate.py`, `suggest.py`, `loop.py` | Fit the multi-fidelity surrogate; propose the next run (expected improvement); the full optimisation loop (`config/bo_config.yaml`). |
| Figures | `figstyle.py` | The shared style (colour = grid level or Kim et al., line style = Kim solid or ours dotted, marker = quantity). New figure scripts must use it. |

## Layout

The top level of `scripts/` holds only the **pipeline and the modules it imports**: everything that the
Makefile, `config/`, `docs_site/`, `tests/` or the packages below use, and what those import in turn.
Every other script is in a sub-folder, by what it does:

| Folder | Files | What is in it |
|---|---|---|
| `run/` | 33 | Scripts that submit or watch SLURM jobs (`submit_*`, `*_watchdog`). |
| `analyse/` | 130 | Analyses behind diary findings (`analyze_*`, `diag_*`, `check_*`, `probe_*`, `compare_*`, ...). They read run data and print or write a small result. |
| `figures/` | 74 | Figures and videos (`plot_*`, `fig_*`, `render_*`). `render_videos.py` is called by the video SLURM template. |
| `experiments/` | 40 | One-off studies that import each other (pre-registered checks, mf-surrogate studies, their figures and launchers). Nothing in the pipeline imports them. |
| `hypotheses/` | 45 | `test_*.py`: **numerical experiments, not pytest tests** (pytest reads only `tests/`), one per pre-registered hypothesis in the lab notebook (for example `test_t4_retro_bo.py` writes `experiments/multifidelity/test_t4_retro_bo.json`). Many submit jobs. |

Some hypotheses reuse functions of others, and `campaign2d/` imports a few (`from scripts.hypotheses.test_h1_mfbml import ...`).

The diary cites scripts by their path before the moves of 2026-11-13 and 2026-11-14 (`scripts/NAME.py`). The new path of each one is in
`internal/docs/history/moved-scripts.md`; `git log --follow` gives a moved file's full history.

## Packages

| Folder | Contents |
|---|---|
| `campaign2d/` | The 2D optimisation campaign: rules for which runs to submit, warm starts, fit gates (`campaign.py`), and its checks. |
| `bo2d/` | Bayesian optimisation over the 2D campaign data: the loop (`bo.py`, `loop.py`), level-transfer studies (`diag_*`), figures (`fig_*`). |
| `mfbml_local/` | The multi-fidelity surrogate KRR-LR-GPR (Yi, Cheng and Bessa) and its Bayesian variants. |
| `kim_prereg/` | Pre-registered comparison with the original Kim et al. code (shell, build and patch). |
| `replay_gate_f3/` | Drivers and analysis for the flux record/replay checks at fidelity 3. |
| `hooks/` | Git pre-commit hook (`make hooks` installs it). |

Each script's first docstring line says what it makes. Shared modules that other code imports most: `postprocess`, `figstyle`, `simulate`, `chain`,
`replay_gate`, `hydro_dataset_v2`, `loop`, `dump_fields`, `cost_model`,
`autoextend`, `chi_scale`, `fields`, `eca_hoekstra`, `periods`, `lagrangian_mix`.

## Archived scripts

One-off launchers of finished campaigns were removed in 2026-11 and are kept in the git tag
`archive/pre-lean-2026-11`. The list, and how to get one back: `internal/docs/history/archived-files.md`.
