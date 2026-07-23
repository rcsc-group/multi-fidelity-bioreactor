# Checkpoint restart and warm-start chains

## Two different things share the word "checkpoint"

This project uses checkpoint restart for two purposes that look identical at
the file level (a `checkpoint.dump`, restored via `argv[2]`) but mean
something different physically:

**Segment** — restarting the *same* condition purely because SLURM's
walltime cut off a single job before the simulation finished. Nothing about
the physics changes; the run is mathematically meant to continue exactly as
if it had never been interrupted. `omega_b_prev` is never set. This is what
[Your first sweep](../tutorials/first-sweep.md) and
[How to sweep one parameter](../how-to/sweep-one-parameter.md) use.

**Warm-start** — seeding a *different* condition (a different `omega_b`,
`theta_max`, etc.) from another condition's already-developed flow field,
to skip that condition's own cold-start transient. `omega_b_prev` (and the
other `*_prev` fields) are set to whatever condition actually produced the
checkpoint, and a smooth-step ramp carries the forcing from the old
condition to the new one.

Conflating these two — assuming a same-condition segment restart is
"basically free" the way it is for pure wall-clock continuation — is exactly
where this gets subtle, below.

## How the restart ramp actually works

`BioReactor.c`'s restart path:

```c
if (params.t_checkpoint > 0.0) {
    restart_file = argv[2];
    // Smooth-step interpolation starts AT the checkpoint and runs
    // N_RAMP_CYCLES forward. alpha goes 0→1 over
    // [t_checkpoint, t_checkpoint + N_RAMP_CYCLES*T_per_st].
    t_ramp_start = params.t_checkpoint;
    ...
}
```

This branch fires on **any** restart where `t_checkpoint > 0` — it has no
condition checking whether `omega_b` (or anything else) actually changed
from the segment that wrote the checkpoint. A 3-period forcing ramp is
re-triggered every single time a segment restarts, whether it's a genuine
warm-start into a new condition or just a same-condition continuation split
across two SLURM jobs for wall-time reasons alone.

For warm-starts, that's exactly the intended behavior — you *want* the
forcing to ramp smoothly from the old condition to the new one. For a
same-condition segment restart, it's an open question whether this
introduces a spurious transient that a single, uninterrupted run at the
same total duration would never see.

## Why this matters for postprocessing

`postprocess.py`'s quasi-steady-state (QSS) window is `t > t_ramp`, where
`t_ramp` is computed once from the very first ramp (`3 × T_per_nd`,
measured from `t=0`). It has no knowledge of *later* ramps that occur at
each subsequent segment boundary in a multi-segment chain — those all fall
well inside what the QSS window considers "already settled," so if a
restart-ramp transient exists, nothing currently excludes it from KPIs like
`tau_100_max` that are explicitly a *max* over the QSS window (and therefore
maximally sensitive to a brief spike, however small).

![tau_98 across a real 2-segment smoke-test chain: no visible discontinuity at the t≈13.2 restart boundary.](../assets/img/first-sweep-tau98.png)

At fidelity 3 over a few rocking cycles, the same-condition restart in
[Your first sweep](../tutorials/first-sweep.md) shows no visible
discontinuity — which is what "checkpointing is basically free" would
predict. Whether that holds at production fidelity, over the many-period
durations `tau_100_max` is actually computed over, is exactly what the
isolating experiment below is checking — a clean restart at fidelity 3 over
one period doesn't rule out a small ramp transient getting captured by a
*max* statistic over a much longer QSS window.

## Testing the hypothesis

The isolating experiment: run one condition at a fidelity that already has
a *verified-clean* single-shot baseline (fidelity 9), deliberately split
into the same same-condition segment structure a real multi-segment chain
uses, and compare against that baseline.

!!! danger "The first attempt at this used a corrupted baseline — its result doesn't count"
    17.5 RPM was the original choice, compared against L9 run `44133566` as
    the "trusted single-shot baseline." That baseline wasn't trusted enough
    — `44133566` turned out to have been overwritten by an abandoned
    cross-condition warm-start pilot (`t_checkpoint=18.85`,
    `omega_b_prev=2.356194`, seeded from 22.5 RPM), not a genuine cold
    start. The −2.4%/−17.1% numbers from that comparison are retracted, not
    reported here. See [Validating against Kim et al. (2024)](kim-et-al-validation.md)
    for the full story of how this was found. Every condition used anywhere
    else on this site was individually re-verified (`t_checkpoint=None`,
    raw data starting at `t=0`) specifically because this happened —
    reusing an existing "baseline" without checking its actual params.json
    first is exactly the mistake to avoid.

## Resolved: not a real effect

Redone on 30.0 RPM, whose L9 baseline (`488db14b`) was verified clean
*before* use this time: `t_checkpoint=None`, `omega_b_prev=None`,
`shear_stress.dat` starting at `t=0`.

| Metric | Clean single-shot baseline | Same-fidelity, 3-segment chain | Difference |
|---|---|---|---|
| `tau_100_max` | 0.13659 | 0.13245 | **−3.0%** |
| `tau_mean_max` | 0.0010078 | 0.0009818 | **−2.6%** |

Both small — neither remotely close to the −17.1% `tau_mean_max` gap the
retracted (corrupted-baseline) 17.5 RPM experiment showed. That result was
almost certainly an artifact of comparing against a cross-condition
warm-started run, not evidence of genuine restart contamination. On a valid
comparison, **same-condition checkpoint restart segmenting does not
meaningfully perturb either metric** — both differences here are consistent
with ordinary run-to-run/restart numerics, not a systematic effect.

That rules checkpointing out as the explanation for the L9-vs-L10
`tau_100_max` sign flip (see [Validating against Kim et al. (2024)](kim-et-al-validation.md))
more firmly than before. The mesh-fidelity change (9 → 10) itself, or
something else specific to L10, remains the open explanation.
`experiments/l9_l10_checkpoint_isolation_test_30rpm/` has the full manifest
and raw results.

## `n_mix_cycles` vs `n_transition_cycles`

Related but separate: fresh runs (segment 0) use `n_mix_cycles` (typically
80) rocking cycles before oxygen injection, to let the flow field develop
from rest. Restart segments use the much shorter `n_transition_cycles`
(typically 10) instead, because the flow is *already* developed — that
assumption is exactly what makes chained sweeps 70–90% cheaper than running
every condition cold. It's also exactly the assumption that same-condition
restart-ramp contamination would undermine if it turns out to be real.
