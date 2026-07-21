# How to choose a fidelity level

The full runtime table is in the [Fidelity guide](../reference/fidelity-guide.md).
This page is the decision, not the table.

## Default recommendations

| What you're doing | Fidelity | Why |
|---|---|---|
| Writing or debugging a new config/script | 3 | Seconds per run. Correctness only — don't read the physics |
| Smoke-testing a sweep before submitting long jobs | 3 | See [Your first sweep](../tutorials/first-sweep.md) |
| BO loop low-fidelity screening | 5 | Fast enough to run 8–10× more samples than HF; cheap enough that a systematic LF bias doesn't matter — the multi-fidelity surrogate corrects for it |
| BO loop high-fidelity target, standard production sweeps | 7 | The point where the physics is trustworthy and cost is still tractable at scale |
| Grid-convergence check | 8 | Only to confirm 7 hasn't left resolution artifacts — not a routine fidelity |
| Publication-grade reference | 9–10 | Days per condition; reserve for a small number of validated conditions |

## The one thing that actually matters: don't skip the smoke test

Every fidelity ≥ 7 job you submit should have already been dry-run at
fidelity 3 with the same config, checking that:

- restart segments' `logstats.dat` starts at `t > 0` (see [Diagnose a stalled or failed chain](diagnose-stalled-chain.md))
- `results.json` has finite kLa/tau values, not NaN

A stale binary or a bad config produces the exact same failure mode at
fidelity 3 in seconds as it would at fidelity 7 after two hours — there's no
reason to find out the expensive way.

## If you're not sure fidelity 7 is actually converged

Don't guess — run the grid-convergence test:

```bash
uv run python -m pytest tests/verification/test_grid_convergence.py -m hpc
```

This compares velocity RMS at fidelity 5 vs. fidelity 6 and needs a real
SLURM allocation (~20–40 min). It's marked `hpc`, so it never runs in CI —
see [Test suite](../testing.md).
