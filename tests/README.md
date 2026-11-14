# tests/

```bash
uv run python -m pytest tests/        # fast unit tests only: no binary, no SLURM, about a minute
```

`pyproject.toml` excludes the slower groups by default. They have markers:

| Marker | What it does | Where to run it |
|---|---|---|
| (none) | Pure Python checks of parsing, postprocessing, surrogates, campaign rules. | Anywhere. |
| `medium` | Builds the solver and runs a short level-3 simulation (~2 min). | CI (GitHub Actions). |
| `hpc` | Submits real SLURM jobs at level 5-7 and waits for them. | On the cluster, on purpose only. |
| `slow` | Takes more than a minute. | When you change the code it covers. |

| Folder | Contents |
|---|---|
| `tests/test_*.py` | Unit tests. `test_src_events.py` checks the solver source: no bare `return;` in an event. |
| `verification/` | Physics and numerics checks: mass conservation, boundedness, restart equals continuous run, grid convergence, comparison with Kim et al. Most are `medium` or `hpc`. |
| `integration/` | End-to-end runs of the binary. |
| `fixtures/` | Small reference data (Kim et al. upstream outputs, one geometry case). |

## Note: `scripts/test_*.py` are not tests

Files named `scripts/test_*.py` are numerical experiments (one per hypothesis in `diary.md`); many submit
SLURM jobs. pytest must never collect them, so `testpaths = ["tests"]` and `scripts` is in
`norecursedirs`. Do not move them into `tests/`. `tests/test_test_*.py` are unit tests *of* those
experiment scripts.

## Solver changes

A change to `src/` must also pass the golden-run gate: `tools/golden/` (see `src/README.md`).
