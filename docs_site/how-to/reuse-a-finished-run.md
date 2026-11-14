# How to reuse a finished run

Every submission through `scripts/simulate.py` (also `make submit`, sweeps, chains and the campaign
scripts) gets a **run key**. Two runs with the same key give the same results, so the second one is
not submitted.

## What the key contains

The key is a sha256 of three things (`scripts/run_key.py`):

| Part | What | Not in the key |
|---|---|---|
| Physics | every entry of `params.json` | `run_id`, and keys that start with `_` (bookkeeping) |
| Binary | sha256 of the binary file the job runs (`_binary`, else the template's default) | its name or folder |
| Parent | for a restart: the key of the run whose `checkpoint.dump` it reads | |

Floats count to 12 significant digits, and `5.0` is the same as `5`.

The key goes into the run's `params.json` as `_run_key`, and into the index `runs/_index.jsonl`
(one line per run).

## When a submission stops

| Case | What you get |
|---|---|
| A run with the same key has `results.json` | `DuplicateRun` naming that run: use its results |
| A run with the same key is still pending or running in SLURM | `DuplicateRun` naming that run: wait for it |
| A run with the same key failed (no results, job gone) | nothing: the submission goes ahead |
| The run folder already has `results.json` from other params | `ValueError`: pick a new `run_id` |
| The binary is not on disk | `ValueError`: build it or set `_binary` |

From the command line, `simulate.py --slurm` prints the existing run and exits with code 3.

In Python:

```python
from scripts.simulate import submit_slurm, DuplicateRun

try:
    job = submit_slurm(params)
except DuplicateRun as e:
    run_dir = e.run_dir          # read e.run_dir / "results.json" when e.finished
```

## A deliberate repeat

A determinism or rank-invariance check repeats a run on purpose:

```python
submit_slurm(params, allow_duplicate=True)
```

or `simulate.py PARAMS --slurm --allow-duplicate`.

## Old runs and the index

`python -m scripts.run_key index runs` rebuilds the index from every `runs/*/params.json`. Run it
as a SLURM job: it reads every run folder and hashes every binary.

An old run's key is usable for this check only if

- its binary is still on disk and did not change after the run (a binary rebuilt under the same name
  is a different code), and
- for a restart, `params.json` names its parent (`_parent_run` or `_source_run`) and the parent's key
  is usable.

Other old runs are in the index with `"dedupe": false`. They never stop a submission.
