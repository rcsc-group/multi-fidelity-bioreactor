# APS DFD 2026 plan (v4): what multi-fidelity CFD can and cannot learn about a rocking bioreactor

Status 2026-10-21. Talk about 2026-11-22. v3 (2026-10-12) said "scalars are numerical diffusion". The PI meeting of 2026-10-15
rejected that framing: the goal is to predict L10 quantities from lower levels at a fraction of the cost, convergence to the
physical value is secondary. This version rewrites the plan around that goal. Every claim names the test that supports it
(diary.md, section and date). A claim without a test is marked HYPOTHESIS. Numbers are from the diary; confounds are listed.

## One-sentence thesis

In a 2-D VOF model of a rocking bioreactor, a chain of multi-fidelity models (mfbml, Yi et al. 2407.15110) predicts the
L10 wall shear and dissipation to 2-7 % on points it never saw, provided an L9 run exists at the same point; the scalar
quantities (kLa, mixing time) are rough in the rocking speed, are not converged with the grid, and can be predicted from L9
to about 17-19 % by a simple ratio as well as by the chain.

## Evidence

### Hydrodynamics (shear, dissipation)
- H-1. Chain L8 -> L9 -> L10 on ln(QoI), 13 L10 runs, leave-one-out: relative RMSE tau_95 4.6 %, tau_mean 3.0 %, EDR 4.1 %; 12 of 13
  inside the 95 % interval (diary 2026-10-14, H3). The two-level L8 -> L10 chain FAILS (15.6 / 24.3 %).
- H-2. The published interval was too narrow because it ignored the L9 uncertainty. With it propagated: half-widths 16 / 23 / 47 %
  (diary 2026-10-19, band check). Without the L9 run at the held-out point the EDR error is 28 % (tau_95 7.8 %, tau_mean 11 %).
- H-3. TRUE TEST SET: 5 of 6 new L10 runs at theta 5.5 (between the training angles), no refit: RMSE tau_95 2.0 %, tau_mean 3.1 %,
  EDR 7.3 %, coverage 5/5 (interim, 2026-10-21). The simple baseline "L9 x one ratio" gives 3.4 / 3.0 / 4.0 %: it matches or beats
  the chain except for tau_95. The 6th run is resubmitted (timeout at 13 h).
- H-4. There is a narrow dissipation bump near 22.5 rpm (theta 7): present at L9 and L10, spike of ~0.6 rpm width at L8
  (outlier screen and width test, 2026-10-18/19). At L9 the pair 22.4 / 22.6 rpm differs by 30 % in EDR and 3.6x in late mixing rate.

### Scalars (kLa, mixing time)
- S-1. The scatter between runs 0.2 rpm apart is in the simulated flow: not our estimator (test S), not our definitions (Kim's
  MATLAB on our fields, KP), and the same jumps appear in Kim's own pipeline (KS) (diary 2026-10-15). Median deviation from a
  smooth trend after removing 5 outliers of 63: 3-12 % (diary 2026-10-19).
- S-2. kLa at C* = 0.10 has the same L8 scatter at theta 4, 7, 9 (0.11 / 0.095 / 0.11 in ln), releases the tracer at cycle 12 or 80
  with no change (+0.6 +- 3 %), responds smoothly to tiny rpm changes. The 95 % mixing time does not: scatter 1.6-1.9x larger off theta 7,
  -43 % with a short spin-up, two families of late mixing rates at theta 9 (diary 2026-10-17/19). HERO TARGETS: kLa(C* = 0.10) and the 50 % mixing time.
- S-3. L10 kLa(C* = 0.10) at 8 + 7 points: a constant L8 ratio FAILS (RMS ln 0.82); L10/L8 is 0.14-0.55, depends on rpm and theta; at
  theta 7 and 17.5 rpm kLa falls 37 -> 15 -> 8.3 1/h from L8 to L10 (no convergence) (diary 2026-10-19/20).
- S-4. L8 -> L9 (theta 7, 14 L9 runs, leave-the-rpm-group-out): chain RMS ln 0.22 (kLa) and 0.22 (50 % mixing time); baselines L9-only
  0.26 / 0.38, L8 x ratio 0.25 / 0.33 (diary 2026-10-19; the first version leaked and is disclosed).
- S-5. L9 -> L10 (theta 7, 7 L10 runs): chain 0.19, L9 x ratio 0.17: the chain does not beat the ratio (diary 2026-10-21).
- S-6. FIXED-FLOW SCALAR-MESH TEST (exact record and replay on an L8 flow, 32.5 rpm): scalar mesh L8 / L9 / L10 gives 50 % mixing time
  11.8 / 19.0 / 24.7 s and kLa(0.25) 28.6 / 38.7 / 51.6: +30 % per level, NOT converged (diary 2026-10-19). Flow resolution and
  scalar resolution both drive the coupled ladder.
- S-7. L9 binary and protocol: same-protocol L9 runs reproduce the older cold runs to 0.01 % (diary 2026-10-20). Closed.

## Open tests (pre-registered, diary 2026-10-21)
- O-1. Two-dimensional scalar chain: new L9 scalar runs at theta 4, 5.5, 9 (12) and L10 kLa runs at theta 5.5 (5-6) as a TRUE test set; PASS if RMS
  ln error <= 0.25 and below the L9 x ratio and the L10-only baselines. Queue: 17 jobs pending, first start 2026-10-21 14:37.
- O-2. Hydrodynamic test set, 6th point (resubmitted, 39 h).

## Talk structure (12 minutes)
1. Setup (1 min): 2-D VOF rocking bag, levels L6-L10, the cost ladder (L8 ~0.5 core-h per cycle ... L10 ~30-80 core-h per cycle), the design space (rpm, theta).
2. Shear and dissipation transfer (3 min): H-1, H-3 figure (prediction and test set), the honest result H-3 (L9 is the key ingredient).
3. Scalars: what transfers and what does not (5 min): S-2 (hero targets), S-1 (roughness is in the flow, same in Kim's code), S-3 / S-4 / S-5 (predictability
   ladder), S-6 (flow vs scalar resolution: the record-and-replay tool).
4. Cost and consequences (2 min): cost per level and per point, the learning curve (TO DO), what an optimiser can trust.
5. Limits (1 min): 2-D, one geometry, 7 L10 kLa runs at theta 7, warm-start protocol at L10.

## Remaining work (about 32 days left)
1. O-1 and O-2 (queue, 2-4 days). 2. Learning curve: error vs number of L10 runs and core-hours (zero compute, after O-1). 3. Cost table per level from logs.
4. Adversarial Sonnet 5.5 review of this plan and the figures (start now on v4, again before the slides). 5. Figures in the figstyle grammar: levels, scatter, prediction + test set,
   predictability ladder, fixed-flow scalar mesh. 6. Slides and a 12-minute rehearsal.

## Risks and known confounds
- L10 kLa runs are warm restarts released after 1 cycle (binary lean-eec1f15); L8 and L9 are cold, released at cycle 80. Direct L10 check: 32.5 rpm theta 7, 22.05 vs 22.43 1/h only.
- Few L10 scalar points (15 now, 6 more coming); run-to-run spread at L10 is not measured (one run per point).
- The chain's advantage over "L9 x ratio" is small or absent; the talk must not oversell mfbml.
- The "narrow regime at 22.5 rpm" mechanism is a hypothesis (resonance); only its width at L8 and L9 is measured.
- Reviewer objection "2-D": all claims are for this 2-D model.
