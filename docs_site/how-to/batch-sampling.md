# How to run a batch design of experiments

Use this to generate a space-filling design of experiments (DoE) across the
parameter space without running the full BO loop — e.g. to seed a surrogate,
or to explore a region before deciding whether optimization is worthwhile.

## 1. Write a sample config

Copy `config/sample_config.yaml` and edit:

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

## 2. Submit

```bash
python scripts/sample.py config/sample_config.yaml
```

All runs are submitted in parallel (independent — no checkpoint restart
between them, unlike the chained sweeps). Results are stored in an
[f3dasm](https://github.com/bessagroup/f3dasm) `ExperimentData` store at
`experiment_dir/experiment_data/`.

## 3. Inspect results

```python
from f3dasm import ExperimentData
data = ExperimentData.from_file("experiments/lhs_run_001/experiment_data")
print(data)
```

## Next

- [Run or resume the BO loop](run-bo-loop.md) — this DoE store's structure is exactly what `loop.py` builds on
