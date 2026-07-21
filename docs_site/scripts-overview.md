# Scripts at a glance

| Script | How to run | What it does |
|--------|-----------|--------------|
| `simulate.py` | `uv run python simulate.py params.json --slurm --wait` | Submit or run one simulation; wait for `results.json` |
| `launch.py` | `uv run python launch.py params.json` | Set up a run directory and write a SLURM script, but do **not** submit |
| `chain.py` | `uv run python chain.py config.yaml` | Sweep **one** parameter across values via a YAML config, using checkpoint restart between segments |
| `sweep.py` | `uv run python sweep.py config.json` | Sweep **any combination** of parameters from a JSON file (zip or cartesian); groups runs by checkpoint compatibility; chains self-submit |
| `sample.py` | `uv run python sample.py config.yaml` | Space-filling batch of independent runs (Latin Hypercube or random) |
| `loop.py` | `uv run python loop.py config.yaml` | Full Bayesian optimisation loop (DoE → train surrogate → suggest → repeat) |
| `postprocess.py` | `uv run python postprocess.py runs/my_run/` | Extract kLa, mixing time, and vorticity from run output files → `results.json` |
| `collect_results.py` | `uv run python collect_results.py --sweep config/my_sweep.json` | Aggregate all `results.json` files into a single CSV |
| `plot_heatmaps.py` | `uv run python plot_heatmaps.py` | Generate KPI heatmap figures from all completed sweep results → `experiments/figures/` |
| `render_videos.py` | `uv run python render_videos.py runs/my_run/` | Convert raw frame dumps to MP4 videos (called automatically by SLURM jobs) |
| `suggest.py` | `uv run python suggest.py exp_dir param_space.yaml` | Print the next highest-EI parameter point to stdout |
| `train_surrogate.py` | `uv run python train_surrogate.py exp_dir model.pkl` | Train the multi-fidelity surrogate from existing run data |

See [Workflows](workflows/single-run.md) for how these scripts fit together end to end.
