# How to run a single simulation

For running the binary directly and seeing every intermediate file, see
[Your first simulation](../tutorials/first-simulation.md). This page is the
shortcuts for doing it routinely.

## Locally (blocks your terminal — use fidelity 3–4 only)

```bash
make run PARAMS=runs/my_run/params.json
```

Equivalent to `cd runs/my_run && ../../build/BioReactor params.json`, with
one difference worth knowing: `make run` depends on `build/BioReactor` and
will rebuild it first if the source is stale.

## Via SLURM (non-blocking — production)

```bash
make submit PARAMS=runs/my_run/params.json
```

Or, with Python — which also waits for `results.json` and returns it:

```bash
uv run python scripts/simulate.py runs/my_run/params.json --slurm --wait --walltime 04:00:00
```

Preview the `sbatch` command without submitting anything:

```bash
make submit PARAMS=runs/my_run/params.json DRYRUN=1
```

## Set up the run directory without running or submitting anything

`launch.py` writes `params.json` and a SLURM script into the run directory
and stops there:

```bash
python scripts/launch.py path/to/params.json [runs_root]
# returns {"run_id": "...", "run_dir": "...", "slurm_script": "..."}
```

Useful when you want to inspect or hand-edit the generated SLURM script
before it runs.
