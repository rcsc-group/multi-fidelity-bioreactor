# Workflow C — Batch sampling

Use this to generate a space-filling design of experiments (DoE) across the parameter
space without running the BO loop.

## Step 1 — Write a sample config

Copy and edit `config/sample_config.yaml`:

```yaml
experiment_dir: experiments/lhs_run_001

sampling:
  strategy: latin    # latin | random | grid | sobol
  n_samples: 20
  seed: 0

fidelity: 5          # 32×32 grid, ~5–15 min each

# t_end is computed per-run as t_mix(params) + t_buffer (adapts to each frequency)
t_buffer: 150

submit: true
walltime: "00:30:00"
```

## Step 2 — Submit

```bash
python scripts/sample.py config/sample_config.yaml
```

All runs are submitted in parallel. Results are stored in an
[f3dasm](https://github.com/bessagroup/f3dasm) `ExperimentData` store at
`experiment_dir/experiment_data/`.

## Step 3 — Inspect results

```python
from f3dasm import ExperimentData
data = ExperimentData.from_file("experiments/lhs_run_001/experiment_data")
print(data)
```
