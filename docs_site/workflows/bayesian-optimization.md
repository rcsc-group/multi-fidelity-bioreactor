# Workflow D — Multi-fidelity Bayesian optimisation

Maximises kLa over the parameter space using a multi-fidelity surrogate
(kernel-ridge regression + linear regression + Gaussian process).
Low-fidelity (fidelity 5) runs build the correlation structure cheaply;
high-fidelity (fidelity 7) runs target the optimum.

```
Phase 1 — Initial DoE (parallel)
  ├── n_lf_init LHS samples at lf_fidelity
  └── n_hf_init LHS samples at hf_fidelity

Phase 2 — BO loop (n_iter sequential iterations)
  ├── train MF surrogate on all existing data
  ├── maximise Expected Improvement over n_candidates random candidates
  └── submit one HF run at the highest-EI point
```

## Step 1 — Write a BO config

Copy and edit `config/bo_config.yaml`:

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

The parameter bounds are defined in `config/param_space.yaml`.

## Step 2 — Run the loop

```bash
python scripts/loop.py config/bo_config.yaml
```

The loop is **resumable**: it checks how many runs already exist in
`experiment_dir` and skips the DoE phase if complete. Kill and restart freely.

## Surrogate tools (standalone)

**Train the surrogate manually:**
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
Prints a JSON params dict to stdout (the highest-EI candidate).

## Testing the algorithm cheaply

`tests/test_loop_integration.py` exercises the full DoE + BO loop end to end
against a synthetic benchmark objective instead of a real simulation — only
the SLURM submission/wait boundary (`scripts.simulate.submit_slurm` /
`wait_for_result`) is mocked, so the real surrogate training and EI
acquisition run for real, in seconds, with a known optimum to check
convergence against.
