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

**L9** (fidelity 9, single monolithic run per condition — no checkpoint
restart of any kind, cold-started, run straight to `t_end=18.243`),
re-postprocessed with the current code and verified against each
condition's own raw `shear_stress.dat`:

| RPM | `tau_100_max` error vs. Kim | `tau_mean_max` error vs. Kim |
|---|---|---|
| 17.5 | −15.2% | −16.1% |
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
segments, same total duration, final segment only), same code:

| RPM | `tau_100_max` error vs. Kim | `tau_mean_max` error vs. Kim |
|---|---|---|
| 17.5 | +30.7% | −35.8% |
| 30.0 | +52.4% | −42.4% |
| 32.5 | +22.4% | −46.9% |

![tau_100_max vs. rocking speed: Kim et al.'s reported curve in black, our L9 (fidelity 9, single-shot) in blue consistently below it, our L10 (fidelity 10, checkpoint-chained) in red consistently above it at every RPM where both exist -- the dotted lines connecting L9 and L10 at the same RPM show the gap and its sign flip directly.](../assets/img/kim-validation-tau100max.png)

## The part that "smells"

Going from L9 to L10 should mean *only* a mesh refinement — same total
simulated duration, same physical condition, higher resolution. Instead,
`tau_100_max` flips sign at every directly-comparable condition: L9
under-predicts, L10 over-predicts. `tau_mean_max` doesn't flip sign, but
does get somewhat worse. A genuine mesh refinement shouldn't do this on its
own.

Ruled out already, with real experiments, not assumption:

- **Stale-metric-definition mismatch** (comparing an old whole-run-max
  `tau_100_max` against a new QSS-window-restricted one) — L9 was
  re-postprocessed with the exact same current code as L10, so this isn't a
  bookkeeping artifact.
- **Simple run-to-run nondeterminism** — the same environment, rerun
  identically 3×, always gives the same result.
- **Basilisk source-version difference** between the persistent OSCAR build
  and a fresh tarball build — a fresh build still reproduces L9's result.
- **Same-condition checkpoint restart itself, for `tau_100_max`** — tested
  directly: one condition (17.5 RPM) held at L9's own fidelity, deliberately
  split into the same 4-segment structure L10's chains use. Result:
  `tau_100_max` moved by only −2.4% relative to the untouched single-shot
  baseline (0.07901 vs. 0.08097) — nowhere near enough to explain L9
  under-predicting Kim by 15% while L10 over-predicts by 30%+ at the same
  condition. See [Checkpoint restart and warm-start chains](checkpoint-restart.md#resolved-a-real-but-partial-effect)
  for the full numbers.

Partially confirmed: the same experiment found checkpoint restart **does**
measurably degrade `tau_mean_max` (−17.1% at the same fidelity, same
condition) — a real contaminant, just not the one responsible for the
`tau_100_max` sign flip. That leaves the mesh-fidelity change itself (or
something else specific to L10 that isn't checkpointing) as the leading
open explanation for `tau_100_max` — not yet identified.

## What this means if you're using these numbers

Don't treat either sweep's `tau_100_max`/`tau_mean_max` as validated against
Kim et al. yet. If you need a shear-stress KPI for a real decision today,
prefer `tau_98_qss`/`tau_100_qss` (the median-over-QSS-window metrics,
robust to a single transient spike) over the `_max` variants — the `_max`
variants are the ones with an unresolved, condition-flipping discrepancy.
