# APS DFD 2026 plan (v3): what multi-fidelity CFD can and cannot learn about a rocking bioreactor

Status 2026-10-12. Talk ~2026-11-22. v2 (2026-10-10) built the talk on three hypotheses (H_rep, H_flow, H_nd) and a
Lagrangian estimator. Since then the pre-registered tests R, P, Q, I, K, L, G and the L8 -> L10 warm start have run
(diary.md, 2026-10-11/12). This version replaces v2. Every claim below cites a test; anything without a test is
marked HYPOTHESIS.

## One-sentence thesis

In a 2-D VOF model of a rocking bioreactor, the hydrodynamic design objectives (wall shear, dissipation) are smooth
in the operating parameters and transfer across grid levels, so a multi-fidelity surrogate (mfbml) maps them at the
finest level cheaply; the scalar objectives (mixing time, kLa) are set by numerical diffusion at every affordable
level, so any surrogate of them, single- or multi-fidelity, learns the numerics, not the physics.

## Evidence (all pre-registered; power stated in diary.md)

Hydrodynamics
- E1. tau_mean, tau_95, EDR are smooth in rpm at L8 (T5). mfbml L8 -> L10 (T2b, full-Bayes residual GP, 3 L10
  points): rel. RMSE tau_mean 9.7%, tau_95 5.7%, EDR ~20%; 95% coverage 97%; baselines L8 x ratio 17-32%, L10-only
  GP 25-27%. (T2b missed its own RMSE rule by a hair: 1.22 vs 1.2; stated.)
- E2. Cheap L10 points: L8 -> L10 warm start passes at 25 and 32.5 rpm (all QoIs <= 3% vs the L9-started
  reference; at 32.5 rpm also vs a cold L10). It FAILED at 17.5 rpm on tau_mean (+7.3%) because tau_mean was still
  falling after 6 cycles; tau_95 and EDR pass everywhere. At 25 rpm a 6-cycle continuation shows the fall stops at
  cycle ~7 and settles at -2% (x8e). The 17.5 rpm continuation decides the cycle count (pending, Mon 2026-10-12).

Scalars (mixing time dtmix_chi at chi = 0.50/0.75/0.95; kLa_1T_25)
- E3. Grid ladder at 32.5 rpm, physical D: dtmix_0.95 14 / 28 / 68 / 128 / 191 s (L6-L10); kLa_1T_25
  183 / 80 / 29 / 17 / 12.7 1/h (L6-L10; L10 = l10c, binary f1c11e0); kLa_1T_10 47 / 29 / 22 1/h (L8-L10).
  No sign of convergence.
- E4. Test K (D x 1 ... 1e5 at L6-L9, 22 runs): a 10x change of the physical scalar diffusivity moves dtmix by
  <= 1.7% at L7, L8 AND L9 (a physics-controlled dtmix would move ~13%/decade) and kLa_1T_25 by <= 3.7%
  (exponent <= 0.016; penetration theory 0.5). [Corrected 2026-10-12: an earlier version quoted the phase-noisy
  5-sample kLa_25 under this label.] Formal verdict INCONCLUSIVE: the high-D positive
  control (levels collapse at 1e5 D) holds for chi 0.50/0.75 but fails for chi 0.95 at L9 (unexplained stall;
  detached-pocket hypothesis tested and rejected). Stated as such in the talk.
- E5. Test L: measuring chi at a fixed physical scale (3.9-15.6 mm boxes) keeps the full L8 -> L9 gap
  (ratio 0.94-1.10, 4 cases): the grid dependence is in the transport, not in how chi is measured.
- E6. Test G: the solver's scalar microscale is 3-4 mm (3-6 cells), smooth (93-100% of liquid cells at
  0.1 < c < 0.9); the physical Batchelor scale is 4-10 um. Formal verdict NEITHER (bands too narrow; stated).
- E7. Roughness in rpm: deterministic (R: reruns reproduce dtmix_0.50 exactly, dtmix_0.95 within 0.2-1%; P: 1e-6
  rpm reproduces to 4e-8) but steep: +0.01 rpm switches chi(t) between two families (dtmix_0.95 128 vs 90 s at
  22.5 rpm). Not explained by variance trapped in a persistent region at L8 (I: NOT ISLAND). Periodicity (Q):
  37.5 rpm aperiodic; 22.5/30/32.5 consistent with periodic; low power elsewhere.

## Talk structure (12 min)

1. Setup (1 min): 2-D VOF rocking bag (Kim/Harris/Cimpeanu solver), levels L6-L10, design space rpm x theta.
   All claims are about this 2-D model; stated on the first slide.
2. Hydrodynamics transfer (3 min): E1 + E2. mfbml (Yi et al. 2407.15110, KRR-LR-GPR, Bessa group) L8 -> L10 over
   rpm x theta: maps of tau_95 and EDR at L10 with intervals; leave-one-out over L10 points; baselines L8 x ratio,
   L10-only GP. This is the constructive result and the mfbml content.
3. Scalars do not (5 min): E3 ladder; E4 the D-scan as the central figure (dtmix flat in D, kLa exponent far below
   0.5, both up to L9); E5 and E6 rule out the measurement scale and show the resolved scalar is mm-scale. The
   message: a surrogate of CFD dtmix or kLa at L10 is a surrogate of numerical diffusion.
4. Roughness (2 min): E7. Deterministic steep structure in rpm, not noise and not trapped islands at L8; the
   0.95 threshold is the fragile one (P, I and the K control all fail there).
5. Consequences (1 min): which objectives can be optimised with MF-CFD in this class of flows; what a scalar
   objective needs (a validated sub-grid transfer closure or a D-resolved method), with the kLa exponent as the
   measurable acceptance test for any such method.

Required L10 numbers (success criterion: tau_95, EDR, mixing time, kLa at L10):
- tau_95, EDR at L10: mfbml map over rpm x theta (HF points below) + the 9-rpm L10 set at theta 7.
- dtmix, kLa at L10: the existing l10c (32.5 rpm, cold, released at cycle 80) values, reported WITH E3/E4 as their
  uncertainty statement (they are numerical-diffusion numbers; the talk says so).

## Remaining work (in order)

1. x8e 17.5 rpm verdict (Mon). Decides the HF cycle count: 12 cycles if TRANSIENT; else L8 -> L9 -> L10 chains
   (about 2x cost) or a cold L10 reference at 17.5 rpm.
2. HF design over rpm x theta (T6: 50 L8 points; rpm 17.5-37.5, theta 2-10): 6-8 L10 warm-start points chosen by
   max posterior variance of the L8 -> L10 mfbml residual (sequential, 2 at a time). Budget ~6k core-h. Pre-register
   the leave-one-out acceptance (rel. RMSE and 95% coverage) BEFORE the first HF point. Exclusions stated: theta 2
   (low-frequency variability) and 35-37.5 rpm (EDR still moving after 6 L10 cycles) unless the continuation
   shows they settle.
3. Figures (figstyle.py): ladder E3, D-scan E4 (the key figure), test L, test G, the chi families (E7), mfbml maps.
4. Adversarial review (Sonnet 5.5) of this plan + the evidence, AFTER items 1-2 are pre-registered; fix or drop
   every claim it breaks.
5. Optional, only if cheap: an L10 D x 10 restart from the l10c release checkpoint up to chi = 0.75 (HYPOTHESIS
   that L10 is still flat in D; cost to be estimated from l10c timings before any decision).

## Risks

- Warm start fails at low rpm: chains or cold reference; or restrict the HF map to rpm >= 25 (stated).
- mfbml on the rpm x theta map misses its acceptance rule: report it as measured; E1 at theta 7 stands.
- Reviewer objection "2-D": the scalar result is a statement about grid-scale numerics of VOF + advection-diffusion,
  which carries to 3-D at equal resolution; the hydrodynamic numbers are 2-D only. Both stated.
- Reviewer objection "known that Sc ~ 500 is unresolved": the contribution is the measured consequence (E4: D
  insensitivity through L9, kLa exponent ladder), not the scale estimate; and that MF surrogates inherit it.
