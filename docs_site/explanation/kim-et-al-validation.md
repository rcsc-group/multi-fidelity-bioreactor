# Validating against Kim et al. (2024)

## What Fig. 13a actually reports

Kim et al.'s Fig. 13a reports two shear-stress quantities per rocking
frequency, both taken as a **max over one period in the quasi-steady
regime** — deliberately excluding the startup transient, not a global max
over the whole run:

- **τ_liq_max** (solid circles) — the absolute max of shear stress over
  space *and* time.
- **τ_liq_mean** (hollow circles) — the max *over time* of the
  spatially-averaged shear stress (i.e. `max_t ⟨τ(t)⟩_space`, not a mean
  over time).

This project's direct analogs are `tau_100_max` and `tau_mean_max` in
`results.json` (see [Output files reference](../reference/output-files.md)),
both restricted to the same quasi-steady window — see
[Non-dimensionalization](non-dimensionalization.md) for how the physical
unit conversion works, and
[Checkpoint restart and warm-start chains](checkpoint-restart.md) for how
that window is currently defined and its blind spot in multi-segment chains.

## Where we actually stand

This is written straight, not smoothed over: as of this writing, **neither
our L9 nor L10 sweep reproduces Kim et al.'s reported values within any
reasonable tolerance**, and the discrepancy doesn't even have a consistent
sign.

!!! danger "17.5 RPM L9 data was found compromised while writing this page — excluded"
    `44133566` (17.5 RPM's L9 run) was the pilot for an abandoned
    cross-condition warm-start rerun plan from earlier in this project's
    history. Its `params.json` shows `t_checkpoint=18.85`,
    `omega_b_prev=2.356194` (22.5 RPM's own frequency) — it was actually
    re-executed as a warm-start from 22.5 RPM's checkpoint, which truncated
    and overwrote its original cold-start `shear_stress.dat`
    (`fopen(..., "w")` on every invocation). The genuine cold-start raw data
    for this condition no longer exists. Every other L9 condition used
    below was individually re-verified clean (`t_checkpoint=None`,
    `shear_stress.dat` starting at `t=0`) before being trusted — this one
    wasn't caught until [Checkpoint restart and warm-start chains](checkpoint-restart.md)'s
    isolation experiment needed a clean 17.5 RPM baseline and didn't have
    one. The isolation experiment has been redone on 30.0 RPM instead,
    which passed the same verification.

**L9** (fidelity 9, single monolithic run per condition — no checkpoint
restart of any kind, cold-started, run straight to `t_end=18.243`),
re-postprocessed with the current code and verified against each
condition's own raw `shear_stress.dat`:

| RPM | `tau_100_max` error vs. Kim | `tau_mean_max` error vs. Kim |
|---|---|---|
| 17.5 | *(data compromised — see above)* | *(data compromised — see above)* |
| 27.5 | −28.9% | −39.9% |
| 30.0 | −21.3% | −37.4% |
| 32.5 | −36.2% | −43.7% |
| 35.0 | −50.8% | −50.7% |
| 37.5 | −87.2% | −62.3% |

Under-prediction, worsening at higher RPM, for both metrics. (20.0 and 25.0
RPM have no recoverable raw data; the 22.5 RPM value that once looked
promising — `tau_100_max=0.2376`, ~+3.4% — was traced to a `rsync`
race condition in the SLURM template and is not reproducible from its own
raw data. It's excluded here, not because it's inconvenient, but because
it's been directly falsified.)

**L10** (fidelity 10, same condition chained across 3–4 same-omega_b
segments, same total duration, final segment only), same code — this data
is unaffected by the L9 contamination above (different run_ids entirely):

| RPM | `tau_100_max` error vs. Kim | `tau_mean_max` error vs. Kim |
|---|---|---|
| 17.5 | +30.7% | −35.8% |
| 30.0 | +52.4% | −42.4% |
| 32.5 | +22.4% | −46.9% |

![tau_100_max vs. rocking speed: Kim et al.'s reported curve in black, our L9 (fidelity 9, single-shot) in blue consistently below it, our L10 (fidelity 10, checkpoint-chained) in red consistently above it -- the 17.5 RPM L9 point is the compromised one described above and should be read with that caveat.](../assets/img/kim-validation-tau100max.png)

## The part that "smells"

Going from L9 to L10 should mean *only* a mesh refinement — same total
simulated duration, same physical condition, higher resolution. Instead, at
the two conditions with a verified-clean L9 baseline (30.0 and 32.5 RPM),
`tau_100_max` flips sign: L9 under-predicts, L10 over-predicts.
`tau_mean_max` doesn't flip sign, but does get somewhat worse. A genuine
mesh refinement shouldn't do this on its own.

Ruled out already, with real experiments, not assumption:

- **Stale-metric-definition mismatch** (comparing an old whole-run-max
  `tau_100_max` against a new QSS-window-restricted one) — L9 was
  re-postprocessed with the exact same current code as L10, so this isn't a
  bookkeeping artifact.
- **Simple run-to-run nondeterminism** — the same environment, rerun
  identically 3×, always gives the same result.
- **Basilisk source-version difference** between the persistent OSCAR build
  and a fresh tarball build — a fresh build still reproduces L9's result.

Being tested: whether **same-condition checkpoint restart itself** explains
the `tau_100_max` sign flip. A first attempt at this used 17.5 RPM as the
isolation experiment's baseline — exactly the condition later discovered to
be compromised above, so that result doesn't count and isn't reported here.
The experiment has been redone on 30.0 RPM (verified clean beforehand) — see
[Checkpoint restart and warm-start chains](checkpoint-restart.md#testing-the-hypothesis)
for the current status.

## What this means if you're using these numbers

Don't treat either sweep's `tau_100_max`/`tau_mean_max` as validated against
Kim et al. yet. If you need a shear-stress KPI for a real decision today,
prefer `tau_98_qss`/`tau_100_qss` (the median-over-QSS-window metrics,
robust to a single transient spike) over the `_max` variants — the `_max`
variants are the ones with an unresolved, condition-flipping discrepancy.
