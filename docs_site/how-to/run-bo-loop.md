# How to run or resume the BO loop

For how the surrogate and acquisition function actually work, see
[Multi-fidelity Bayesian optimization](../explanation/multi-fidelity-bo.md).
For a zero-cost way to try the mechanics first, see
[Your first optimization loop](../tutorials/first-optimization-loop.md).
This page is just running it for real, against Basilisk.

## 1. Write a BO config

Copy `config/bo_config.yaml` and edit:

```yaml
experiment_dir: experiments/bo_run_001

lf_fidelity: 5   # 32×32, ~5–15 min
hf_fidelity: 7   # 128×128, ~1–2 h

n_lf_init: 8     # LF samples for initial DoE (rule of thumb: ≈10× free params)
n_hf_init: 3     # HF samples for initial DoE (seeds the transfer-learning correlation)
n_iter: 20       # BO iterations after DoE

kla_key: kLa_25  # objective: kLa_10 | kLa_25 | kLa_50
n_candidates: 2000

walltime: "02:00:00"
job_timeout: 7200    # seconds to wait for results.json per job
t_buffer: 150
```

Parameter bounds come from `config/param_space.yaml` — see
[params.json reference](../reference/params.md).

## 2. Run it

```bash
python scripts/loop.py config/bo_config.yaml
```

The loop is **resumable**: it checks how many runs already exist in
`experiment_dir` and skips the DoE phase if it's already complete. Kill and
restart freely — nothing needs to be cleaned up first.

## Standalone surrogate tools

Useful when you want to inspect or reuse the surrogate without running the
whole loop.

**Train manually:**
```bash
python scripts/train_surrogate.py experiments/bo_run_001 surrogate/model.pkl \
    --lf-fidelity 5 --hf-fidelity 7 --kla-key kLa_25
```
Writes a pickled `model.pkl` with a `.predict(X) -> (mean, var)` interface.

**Query the next suggested point:**
```bash
python scripts/suggest.py experiments/bo_run_001 config/param_space.yaml \
    --model-path surrogate/model.pkl --kla-key kLa_25 --n-candidates 2000
```
Prints a JSON params dict to stdout — the highest-EI candidate.
