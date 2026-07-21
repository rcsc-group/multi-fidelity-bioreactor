# Your first optimization loop

The full multi-fidelity Bayesian optimization loop — DoE, surrogate training,
Expected Improvement acquisition, repeat — needs no Basilisk build and no
SLURM access to try. [`examples/synthetic_bo_demo.py`](https://github.com/rcsc-group/multi-fidelity-bioreactor/blob/main/examples/synthetic_bo_demo.py)
runs the *real* `loop.py`, `train_surrogate.py`, and `suggest.py` — only the
simulation itself is swapped for a cheap synthetic function with a known
answer, so you can watch the algorithm work in seconds and check whether it
actually found the right thing.

## Run it

```bash
uv run python examples/synthetic_bo_demo.py
```

## What's being optimized

The script's objective is a smooth bump in `(omega_b, fill_level)` peaking at
`omega_b=4.0, fill_level=0.5` (`kLa_25=1.0`), evaluated at two fidelities —
low-fidelity readings are deliberately biased `-0.15` below the true value,
so the multi-fidelity surrogate's bias correction has something real to
correct.

## Expected output (abridged)

```
Phase 1: initial DoE (4 LF + 2 HF runs)
  DoE run 1/6: fidelity=5, omega_b=2.69, t_end=197.1
    kLa_25=0.2494
  DoE run 2/6: fidelity=5, omega_b=2.07, t_end=191.8
    kLa_25=-0.0782
  DoE run 3/6: fidelity=5, omega_b=3.87, t_end=195.6
    kLa_25=0.6949
  DoE run 4/6: fidelity=5, omega_b=5.20, t_end=192.4
    kLa_25=-0.0555
  DoE run 5/6: fidelity=7, omega_b=4.55, t_end=192.1
    kLa_25=0.5164
  DoE run 6/6: fidelity=7, omega_b=6.25, t_end=199.4
    kLa_25=0.0209

Phase 2: 4 BO iterations (maximising kLa_25)

  Iteration 1/4  (best so far: 0.69494)
  Suggested: omega_b=3.86, fill_level=0.32, theta_max[0]=6.7
  Result: kLa_25=0.20827
  ...
============================================================
Optimisation complete.  Best kLa_25 = 0.69494
```

## The honest part

Look closely: the best value found (`0.69494`) came from DoE run 3, not from
any of the four BO iterations that followed. That's real, seeded,
deterministic behavior, not a bug — with only `n_candidates=200` and 4
iterations, Expected Improvement is still exploring (iteration 2 tried
`omega_b=1.76`, nowhere near the optimum) rather than having converged. This
is exactly the exploration/exploitation trade-off EI is designed to make —
see [Multi-fidelity Bayesian optimization](../explanation/multi-fidelity-bo.md)
for why. Try raising `n_iter` and `n_candidates` in the script and re-running
to see the best value actually improve past the DoE.

## Why this testbed exists

`loop.py`'s DoE-then-BO orchestration, `train_surrogate.py`'s KRR-LR-GPR
fitting, and `suggest.py`'s EI acquisition had never actually been exercised
together end to end before this testbed was written — writing it surfaced
two real bugs (an `ExperimentData` domain mismatch on every append, and a
crash in the final summary printout) that unit tests on the individual
pieces had missed. If you're changing anything in the BO loop, run this
first — it's `tests/test_loop_integration.py` under the hood, and it's part
of the default fast test suite.

## Next

- [Run or resume the BO loop](../how-to/run-bo-loop.md) — doing this for real, against Basilisk
- [Multi-fidelity Bayesian optimization](../explanation/multi-fidelity-bo.md) — how the surrogate and acquisition function actually work
