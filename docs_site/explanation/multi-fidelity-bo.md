# Multi-fidelity Bayesian optimization

## The problem it's solving

A high-fidelity run (fidelity 7) is accurate but expensive — hours per
condition. A low-fidelity run (fidelity 5) is cheap — minutes — but
systematically biased relative to the true, converged answer, not just
noisier. Running only HF samples wastes the cheap ones; running only LF
samples never finds the true optimum, since you're optimizing the wrong
(biased) function. Multi-fidelity Bayesian optimization uses many cheap LF
samples to learn the *shape* of the response surface, and a smaller number
of expensive HF samples to learn the *correction* — where and how much the
LF surface is wrong.

## The surrogate: KRR-LR-GPR

`train_surrogate.py` fits a **K**ernel **R**idge **R**egression + **L**inear
**R**egression + **G**aussian **P**rocess **R**egression model (Yi et al.)
to the accumulated LF and HF data:

```python
model.train(X=[X_hf, X_lf], Y=[y_hf, y_lf])
```

Conceptually: KRR learns the LF response surface cheaply from the many LF
points; a linear correction term is fit from the (few) HF points against
that LF surface; the GP layer provides calibrated uncertainty (`mean`,
`var`) at any query point, which is what makes an acquisition function like
EI possible at all — without a variance estimate, there's no way to
quantify "how much could this point still surprise us."

The vendored copy of this model
(`scripts/mfbml_local/krr_lr_gpr.py`) is deliberately numpy-only, with no
torch dependency, even though `torch` is listed in `pyproject.toml`'s main
dependencies (pulled in transitively by `f3dasm`/`mfbml`/`mfpml`, not used
by this project's own training code path).

## The acquisition function: Expected Improvement

`suggest.py._expected_improvement` is the standard maximization form:

```python
improvement = mean - y_best - xi
EI = improvement * Φ(Z) + std * φ(Z)     (0 where std == 0)
```

`y_best` is the best value observed so far — the incumbent EI is trying to
beat. `xi` (currently `0.01`) is an exploration margin: a candidate has to
promise *more* than a small margin over the incumbent before it's rewarded
purely for its predicted mean, which keeps EI from over-exploiting a single
lucky point. Given this project's actual kLa ranges (roughly 0.002–0.02 in
the synthetic testbed's units), `xi=0.01` is a substantial fraction of the
whole output range — worth an explicit ablation before trusting BO runs
where exploration/exploitation balance matters a lot.

## A concrete illustration of why "incumbent" needs care

Until recently, `y_best` had a real bug: when no HF observations existed
yet, it fell back to the raw *LF*-observed maximum as the incumbent. That's
exactly the wrong thing to do — LF and HF are expected to disagree
systematically (that's the entire premise of training a bias correction in
the first place), so using an LF-scale value as an HF-scale incumbent
silently corrupts the very first acquisition decisions. The fix: when no HF
data exists, evaluate the *already-trained surrogate* (not the raw
observations) at the LF-observed inputs, and use its bias-corrected
prediction as the incumbent instead. This is a good example of a
multi-fidelity-specific failure mode that a single-fidelity BO
implementation would never surface — "best observed value" stops being a
simple `max()` the moment two different fidelities enter the picture.

## Testing this cheaply

See [Your first optimization loop](../tutorials/first-optimization-loop.md)
— the whole DoE→train→acquire→repeat loop runs against a synthetic
objective with a known optimum and a deliberate LF/HF bias, in seconds, with
no Basilisk build required. Two more real bugs (an `ExperimentData` schema
mismatch, and a crash in the final summary printout) were found by actually
running this loop end to end for the first time, rather than only testing
its pieces in isolation.
