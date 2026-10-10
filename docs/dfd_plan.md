# APS DFD 2026 plan: L10-accurate tau_95, EDR, mixing time and kLa over the design space

Status 2026-10-10. Talk ~2026-11-22. Draft for adversarial review.

## Evidence we start from (all in diary.md, pre-registered)

- tau_mean, tau_95, EDR are smooth in rpm (T5: no 0.1-rpm or 0.625-rpm jump >= 10%) and transfer L8 -> L10 with
  vanilla mfbml (Yi et al. 2407.15110, our fully Bayesian residual GP): 6-10% (tau), ~20% (EDR) at 3 L10 points,
  95% coverage 97% (T2b). Baselines: L8 x ratio 17-32%, L10-only GP 25-27%.
- Mixing time and kLa from a single run are NOT smooth functions of the design (T5, L8, theta 7, 25 runs):
  +-0.1 rpm moves dtmix_0.95 by 20-40% (log 0.18-0.35); kLa similar. Release-cycle replicates (sd 0.02 in log)
  understate this 10-20x. So mfbml fails on them (T7: 23-180% error, coverage 72%).
- kLa decreases with refinement at 32.5 rpm (kLa_1T_10: L8 47, L9 29, L10 22 h^-1). Oxygen Sc ~ 500
  (D = 1.9e-9 m^2/s): the interfacial concentration boundary layer is likely unresolved even at L10, so the direct
  L10 kLa is itself a resolution-dependent number (hypothesis, to test).
- Direct L10 mixing/kLa costs ~3600 core-h per design point (cold start + 150 post-release cycles). L10
  hydrodynamics warm-started from a settled lower level costs ~300-550 core-h (6 cycles); L8 -> L10 warm start
  is being validated (jobs 7276284-6).
- Solver already writes uniform-grid snapshots of u, v, vorticity, f, cs, tracers, oxygen in a cycle window
  (params snap_start_cycle / snap_end_cycle, every dt_video). No solver change needed to record the flow.

## Thesis of the talk

"Mixing time and kLa of a rocking bioreactor are properties of its periodic flow, not of one finite scalar run.
Estimate them from the flow, and multi-fidelity design becomes possible at the finest grid for all four QoIs."

## Contributions

C1 (finding). Single-run CFD mixing time and kLa are ill-conditioned design objectives: quantified sensitivity
(T5) vs the smooth hydrodynamic QoIs. Most CFD bioreactor studies (e.g. Kim et al. 2504.05421) report one run per
condition.

C2 (method). Transfer-operator (Ulam / mapping-method) estimator of mixing from recorded periods of the flow:
seed particles in the liquid at phase 0, advect one period through the recorded velocity snapshots, bin the end
points -> period map P (sparse Markov matrix). Mixing curve chi(n) = from P^n applied to the same initial tracer
as the direct run; asymptotic mixing rate from the second eigenvalue |lambda_2| of P. Averaging P over several
recorded periods handles cycle-to-cycle variability. Deterministic given the flow, cheap (GPU post-processing),
and with its own fidelity knobs (flow level, operator cell size, number of periods).

C3 (method). kLa from hydrodynamics: kLa = kL * a with a = interface length per liquid area (from f, exact) and kL
from a surface-renewal / surface-divergence closure evaluated on the recorded interface velocity. This replaces the
unresolved Sc ~ 500 boundary layer by a closure from the gas-transfer literature (constants taken from the
literature, not fitted to our runs). Alternative inside C2: the same operator with an interface exchange term.

C4 (demonstration). mfbml (L8 -> L10) on operator/closure QoIs over rpm x theta (50 L8 points from T6, 6-8 L10
warm-start points), giving L10-accurate maps of tau_95, EDR, dtmix, kLa with calibrated bands, and a
multi-fidelity BO choice of the operating point (max kLa subject to tau_95 bound), against single-fidelity BO.

## Pre-registered falsification tests (criteria fixed before running)

F1 (L8, cheap): record 5 periods at L8 at 5 rpm (theta 7). Operator dtmix_0.50/0.75/0.95 vs the local mean of
  direct L8 runs within +-1.25 rpm (T5 dense data, 4-5 runs each). PASS: operator inside the 95% band of the local
  mean at >= 4/5 rpm.
F2 (smoothness): operator dtmix on the T5 dense rpm set passes T5 criterion (a) (|second difference| <= 3x
  replicate sd except <= 2 narrow features).
F3 (variance): sd of the operator estimate across disjoint recorded-period subsets <= 1/3 of the direct run-to-run
  sd from T5.
F4 (L10 anchor): operator dtmix_0.95 at 32.5 rpm vs direct L10 (190.7 s, l10c_rpm32.5_seg1) and Kim (184 s):
  within 20%.
F5 (kLa): closure kLa vs (a) the direct kLa refinement ladder L8 -> L10 (does the direct value approach the
  closure as the grid refines?) and (b) published experimental kLa for wave/rocking bioreactors if one matches our
  geometry. Criterion to be fixed after the literature is read, before running.
F6 (MF): mfbml L8 -> L10 on the new QoIs over rpm x theta: median rel RMSE below L8 x ratio and L10-only GP, pooled
  95% coverage >= 85%.
F7 (cost): L10 cost per design point with C2/C3 <= 1/5 of a direct L10 mixing/kLa run.

## Compute and schedule

- Week 1: F1-F3 at L8 (recording runs ~16 core-h each; operator code in Python/JAX, GPU). L8->L10 warm-start
  validation result.
- Week 2: F4 at L10 (one warm-start recording run, ~700 core-h). Literature for C3 read; F5 criterion fixed.
- Weeks 3-4: L10 recording runs at 6-8 rpm x theta points (~5k core-h), F5-F7, MF-BO.
- Weeks 5-6: figures, talk.

## Risks

- Ulam coarse-graining adds diffusion; the asymptotic rate must converge with operator cell size (checked).
- Interface motion: the liquid region changes within a period; the period map from phase 0 to phase 0 maps the
  liquid onto itself, but particles must stay in the liquid (projection near the interface).
- If F1 fails, C2 is dropped and the talk falls back to C1 + ensemble-averaged direct runs at L8 and C4 on
  tau/EDR only.
