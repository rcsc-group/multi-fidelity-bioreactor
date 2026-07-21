# Scripts reference

| Script | How to run | What it does |
|--------|-----------|--------------|
| `simulate.py` | `uv run python simulate.py params.json --slurm --wait` | Submit or run one simulation; wait for `results.json` |
| `launch.py` | `uv run python launch.py params.json` | Set up a run directory and write a SLURM script, but do **not** submit |
| `chain.py` | `uv run python chain.py config.yaml` | Sweep one parameter across values via YAML, checkpoint-restarting between segments — [how-to](../how-to/sweep-one-parameter.md) |
| `sweep.py` | `uv run python sweep.py config.json` | Sweep any combination of parameters from JSON (zip or cartesian) — [how-to](../how-to/sweep-json-multi-param.md) |
| `sample.py` | `uv run python sample.py config.yaml` | Space-filling batch of independent runs (LHS / random / grid / Sobol) — [how-to](../how-to/batch-sampling.md) |
| `loop.py` | `uv run python loop.py config.yaml` | Full multi-fidelity BO loop (DoE → train surrogate → suggest → repeat) — [how-to](../how-to/run-bo-loop.md) |
| `postprocess.py` | `uv run python postprocess.py runs/my_run/` | Extract kLa, mixing time, vorticity, and shear-stress KPIs → `results.json` |
| `collect_results.py` | `uv run python collect_results.py --sweep config/my_sweep.json` | Aggregate all `results.json` files into a single CSV |
| `plot_heatmaps.py` | `uv run python plot_heatmaps.py` | Generate KPI heatmap figures from completed sweep results → `experiments/figures/` |
| `render_videos.py` | `uv run python render_videos.py runs/my_run/` | Convert raw frame dumps to MP4 (called automatically by SLURM jobs) |
| `suggest.py` | `uv run python suggest.py exp_dir param_space.yaml` | Print the next highest-EI parameter point to stdout |
| `train_surrogate.py` | `uv run python train_surrogate.py exp_dir model.pkl` | Train the multi-fidelity surrogate from existing run data |
| `stage_segment.py` | `uv run python stage_segment.py <prev_run_id> <next_run_id>` | Safely stage a chain-restart checkpoint — refuses unless the predecessor has a genuine `results.json` |

See [Project structure](project-structure.md) for the full annotated tree, including the diagnostic one-off scripts not listed here.
