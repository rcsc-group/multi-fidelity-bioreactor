# APS DFD 2026 plan (v2): why CFD mixing time and kLa misbehave in a rocking bioreactor, and multi-fidelity design that survives it

Status 2026-10-10. Talk ~2026-11-22. v1 was rejected by an adversarial review (diary 2026-10-10): v1's C2/C3
were known methods (Ulam period map from CFD: Kluenker et al. 2603.13996, 2112.11497; wave-bioreactor kLa
closure with a per-bag fitted constant: Piontek et al., PMC12868168), several tests could not fail, and L10 is
not a converged reference for the scalar QoIs. v2 builds the talk on questions our data already pose.

## Evidence (diary.md, pre-registered tests)

- E1. tau_mean, tau_95, EDR: smooth in rpm (T5: no 0.1- or 0.625-rpm jump >= 10%), change modestly L8 -> L10,
  and transfer L8 -> L10 with mfbml (6-10% / ~20% at 3 L10 points, coverage 97%, T2b).
- E2. dtmix and kLa from one run: rough in rpm (T5: +-0.1 rpm moves dtmix_0.95 by 20-40%) and strongly
  level-dependent at 32.5 rpm: dtmix_0.95 68/128/191 s, kLa_1T_10 47/29/22 h^-1 (L8/L9/L10). mfbml fails on
  them (T7).
- E3. The scalars are unresolved at every level: tracer D = 0.44e-9 m^2/s (Sc ~ 2300), oxygen D = 1.9e-9
  (Sc ~ 500). chi is computed from the cell-scale variance, i.e. at the grid scale, which changes with level.
- E4. Direct dtmix rises and kLa falls with refinement: the sign expected if the scalar's numerical diffusion
  (which falls with the grid) controls them.

## Three hypotheses, each with a decisive test

H_rep (reproducibility): a single run is not reproducible: an exact rerun changes dtmix by more than the
  release-cycle replicate sd (0.02 in log).
  Test R: 8 exact reruns already submitted (rec_l8_*, jobs 7276608-15; identical params except snapshot output).
  Reading fixed: |log rerun - log original| > 0.02 at >= 2 of 8 -> H_rep holds (round-off chaos is a noise source).

H_flow (where the roughness lives): the rpm roughness is Lagrangian structure of the flow, not of the tracer run.
  Test S (same-run): Lagrangian estimator (below) on the recorded flow of each rec_l8 run vs that run's own direct
  chi(t), and across the +-0.1 rpm pairs. Reading fixed:
  - Estimator validity: on the same run, estimator chi(t) at the matched effective diffusivity (below) reproduces
    the direct dtmix_0.50/0.75 within 10% at >= 6 of 8 runs. If not, the estimator is not trusted and S stops.
  - H_flow holds if the estimator at the PHYSICAL diffusivity reproduces >= 3 of the 4 direct +-0.1-rpm jumps in
    sign and within a factor 2 in size; H_flow rejected if its +-0.1-rpm jumps are all < 0.05 in log.

H_nd (what controls the level dependence): the L8 -> L10 change of dtmix comes from the scalar's numerical
  diffusion, not from the flow.
  Test D: the same estimator at the physical diffusivity and a fixed physical measurement scale l (independent of
  the grid) applied to the L8, L9 and L10 flow at 32.5 rpm (and 25 rpm). Reading fixed: H_nd holds if
  |log dtmix_0.95(L10 flow) - log dtmix_0.95(L9 flow)| <= 1/3 of the direct |log 191 - log 128| = 0.40, i.e.
  <= 0.13, AND the estimator's L8->L10 change is < 1/2 of the direct one.

## Estimator (the method we build; prior art cited, not claimed)

Lagrangian particles with a random walk at a prescribed diffusivity D in the recorded velocity field (tank frame;
u, v, f snapshots, 100 per period, time-interpolated): the stochastic form of advection-diffusion, with no grid
diffusion. Tracer initial condition = the direct run's (same release). chi(t) from the particle concentration
coarse-grained at a FIXED physical scale l, defined like the code's chi otherwise. Flow recorded for a few periods
and reused periodically (the Ulam / mapping-method idea of 2112.11497, 2603.13996, here used only to extend time;
the cycle-to-cycle spread is measured by using different recorded periods).
- Matched effective diffusivity: D_eff of the direct run at each level, estimated from the decay of the direct
  tracer variance in the first cycle after release vs the estimator's for a sweep of D; used only for validity.
- Checks (with acceptance thresholds, fixed before use): particle leakage out of the liquid < 1% per period;
  chi(t) change < 5% when particles x4, snapshot rate /2, and l x2 (l reported as a definitional choice).
- Benchmark first: a periodically forced double gyre with known Ulam results (2112.11497) and a pure-diffusion
  case with an analytic answer.

## Talk structure (what each outcome gives)

1. Eulerian vs Lagrangian design objectives: tau_95 and EDR are smooth and transfer across grids (E1, mfbml);
   dtmix and kLa from single runs are rough and grid-dependent (E2). Quantified.
2. Diagnosis by tests R, S, D. Every outcome is a result:
   - H_rep true: single-run CFD mixing times are not reproducible; numbers in the literature need ensembles.
   - H_flow true: the roughness is physical Lagrangian structure (islands/barriers switching with rpm) while the
     Eulerian flow is smooth: mixing time is an intrinsically non-smooth design objective -> robust design
     (expected dtmix over a +-delta rpm band) is the meaningful objective.
   - H_nd true: the grid dependence of CFD mixing time is scalar numerical diffusion; the estimator gives a mixing
     time at the physical diffusivity that converges with the flow level.
3. Multi-fidelity design at the finest level: mfbml (L8 -> L10) over rpm x theta (T6: 50 L8 points; 6-8 L10
   warm-start points) for tau_95, EDR, and the estimator's dtmix if H_nd holds (else the robust ensemble
   dtmix at L8 only, stated). Baselines: L8 x ratio, L8 with fitted affine correction, HF-only thin-plate spline,
   HF-only GP (integrated). Leave-one-out over L10 points; report interval width with coverage.
4. kLa: shown as the convergence ladder (47/29/22) with the Sc ~ 500 explanation. A Lagrangian kLa at the physical
   diffusivity is a stretch goal only, not promised.

## Compute and schedule

- Now: rec_l8 (8 runs, ~450 core-h); T6 (11 left); L8 -> L10 warm-start validation (3 L10 runs).
- Week 1: estimator code + benchmarks (Sonnet agents, Python/JAX, GPU); tests R and S at L8.
- Week 2: build a lean binary from HEAD (cross-level warm start + snapshots; f1c11e0 lacks cross-level);
  L9/L10 recording runs at 32.5 and 25 rpm (warm start, ~6 + 4 cycles); test D.
- Weeks 3-4: L10 points over rpm x theta (~5k core-h), mfbml/MF-BO, figures.
- Weeks 5-6: talk.

## Risks and fallbacks

- Estimator fails validity (S, first bullet): talk = E1 + E2 + test R + mfbml on tau/EDR over rpm x theta.
- L8 -> L10 warm start fails: L8 -> L9 -> L10 chains (about 2x cost).
- 2-D vs 3-D: all claims are about this 2-D model (as Kim et al.); stated on the first slide.
