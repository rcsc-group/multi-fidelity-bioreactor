# scripts/ — Python and shell around the solver

Run every Python script from the repository root, so that `from scripts.X import ...` works:

```bash
uv run python scripts/simulate.py --help
```

Most scripts that run the solver submit SLURM jobs. Do not run them on a login node.

## Start here: the main pipeline

| Step | Script | What it does |
|---|---|---|
| 1. One run | `launch.py`, `simulate.py` | Make a run directory with `params.json` and a SLURM script; run it locally or submit it. |
| 2. Many runs | `sample.py`, `sweep.py`, `chain.py` | A batch from `config/param_space.yaml`; a one-parameter sweep; a sweep where each run restarts from the previous checkpoint. |
| 3. Finish a run | `autoextend.py`, `stage_segment.py`, `periods.py` | Extend a run that hit its time limit; stage a restart segment; where a run stops (same rule as the solver). |
| 4. Results | `postprocess.py`, `hydro_dataset_v2.py`, `health_report.py` | Quantities of interest from one run (kLa, shear stress, mixing time, energy dissipation); the hydrodynamic data set; numerical health checks. |
| 5. Surrogate and optimiser | `train_surrogate.py`, `suggest.py`, `loop.py` | Fit the multi-fidelity surrogate; propose the next run (expected improvement); the full optimisation loop (`config/bo_config.yaml`). |
| Figures | `figstyle.py` | The shared style (colour = grid level or Kim et al., line style = Kim solid or ours dotted, marker = quantity). New figure scripts must use it. |

## Packages (sub-folders)

| Folder | Contents |
|---|---|
| `campaign2d/` | The 2D optimisation campaign: rules for which runs to submit, warm starts, fit gates (`campaign.py`), and its checks. |
| `bo2d/` | Bayesian optimisation over the 2D campaign data: the loop (`bo.py`, `loop.py`), level-transfer studies (`diag_*`), figures (`fig_*`). |
| `mfbml_local/` | The multi-fidelity surrogate KRR-LR-GPR (Yi, Cheng and Bessa) and its Bayesian variants. |
| `kim_prereg/` | Pre-registered comparison with the original Kim et al. code (shell, build and patch). |
| `replay_gate_f3/` | Drivers and analysis for the flux record/replay checks at fidelity 3. |
| `hooks/` | Git pre-commit hook (`make hooks` installs it). |

## Naming of the other top-level scripts

The name prefix tells what a script is. Each script's first docstring line says what it makes.

| Prefix | Count | Kind |
|---|---|---|
| `test_*` | 52 | **Numerical experiments, not pytest tests.** One per hypothesis in `diary.md` (for example `test_t4_retro_bo.py` writes `experiments/multifidelity/test_t4_retro_bo.json`). Many submit jobs. |
| `plot_*`, `fig_*`, `render_*` | 83 | Figures and videos from finished runs. |
| `analyze_*`, `diag_*`, `check_*`, `probe_*`, `compare_*` | ~85 | Analyses behind diary findings. They read run data and print or write a small result. |
| `submit_*` | 33 | Launchers for one campaign or test (the finished ones were archived). |
| `replay_*` | 9 | Flux record/replay: the gate (`replay_gate.py`) and low-fidelity replay (`replay_lf.py`). |

Shared modules that other scripts import most: `postprocess`, `figstyle`, `simulate`, `chain`,
`settling_study_grid`, `replay_gate`, `hydro_dataset_v2`, `loop`, `dump_fields`, `cost_model`,
`autoextend`, `chi_scale`, `fields`, `eca_hoekstra`, `periods`, `lagrangian_mix`.

## Archived scripts

One-off launchers of finished campaigns were removed in 2026-11 and are kept in the git tag
`archive/pre-lean-2026-11`. The list, and how to get one back: `docs/history/archived-files.md`.
