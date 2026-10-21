# Meeting 2026-10-22

Figures are COPIES made on 2026-10-21 (night). Scripts, data and the pre-registered rules are in `diary.md` (entries 2026-10-19 to 2026-10-21). Level colours: L8 blue, L9 green, L10 pink. Conventions: markers show the quantity; filled = a maximum-type statistic, hollow = a mean-type statistic.

Requested by the user for this meeting: (1) the 5x3 chain L8 -> L9 -> L10, (2) the 5x3 usual two-level L9 -> L10, (3) the C* extrapolation problem with mfbml against a GPR with a good prior, (4) a "proper" mfbml where the low fidelity is one converged flow period repeated. Items 5-11 are leads.

## 00  How kLa and the mixing time are defined (answer to "black box")
`00_definitions.png`. One L8 run (25 rpm, 7 deg). Left: C*(t) = oxygen in the liquid / saturation; markers at C* = 0.10, 0.25, 0.50. Right: chi(t) = 1 - sigma^2(t)/sigma^2(0), the variance of the passive tracer in the liquid; markers at chi = 0.50, 0.75, 0.95.
- kLa(C* = x) = the slope of ln(1 - C*) against t over ONE rocking period centred on the first crossing of x (`scripts/postprocess.py`, `_kla_period_at_threshold`). The oxygen curve is not one exponential, so kLa depends on x: median kLa(0.50)/kLa(0.10) = 0.51, kLa(0.25)/kLa(0.10) = 0.71 (88 L8 runs).
- dtmix(chi = y) = first time chi reaches y.
- Checked against the PIs' worry "is it a postprocessing bug?": see card 11.

## 01  Chain L8 -> L9 -> L10, 5 x 3
`01_chain_L8_L9_L10_5x3.png`. Rows: tau_95, mean tau, dissipation rate (EDR), kLa at C* = 0.10, mixing time at chi = 0.50. Columns: rocking angle 4, 7, 9 deg. Markers: simulations. Dashed blue: L8 fit (KRR). Dashed green: L9 fit. Grey line and band: L10 prediction with 95 % range, trained on ALL L10 runs (mfbml, Yi et al. 2407.15110, three levels in a chain, on the log of the quantity; the L9 uncertainty is passed on to the L10 band, `scripts/h3_band_check.py`). Small bars on L10 shear and dissipation markers: spread of the last three cycles.
- Shear and dissipation: 13 L10 runs (7 at 7 deg, 3 at 4 deg, 3 at 9 deg). Leave-one-out error 4.6 / 3.0 / 4.1 % (tau_95 / mean tau / EDR), with an L9 run at the same point; true test set (5 new L10 runs at 5.5 deg) 2.0 / 3.1 / 7.3 %.
- Drawing rule (same in every panel): all level markers are drawn; the fitted curves and the L10 prediction span 17.5-32.5 rpm (the L10 range) and appear only where L9 has >= 3 runs at that angle and the row has >= 5 L10 runs. kLa row: 13 L10 runs, but L9 kLa exists only at 7 deg, so only that column has curves and a prediction.
- Mixing-time row: L10 has only 3 runs (the 50 % mixing time needs more cycles than the L10 runs have; seven runs are being continued, ready in about 2 days). No prediction is drawn.
- Confound: L10 shear and dissipation come from warm starts (from L8 for most points, from L9 at 7 deg); 12 of 13 are settled to within 6 % (check against the end of the continuation runs), the 22.5 rpm run at 7 deg is not (+19 % EDR later).

## 02  Usual two-level mfbml L9 -> L10, 5 x 3
`02_two_level_L9_L10_5x3.png`. Same layout. The low fidelity is a KRR trained on the L9 runs (dashed green); the high fidelity is the L10 GP correction (Yi et al.). No L8 information.
- Same data gaps and the same drawing rule as card 01.
- Honest reading (cards 07 and 03): a regression ratio fitted on L9 is as accurate as the chain; mfbml adds an uncertainty band.

## 03  The C* extrapolation problem (PI idea 1) as a multi-fidelity problem
`03_threshold_extrapolation.png`. Left: C*(t) with the four thresholds (0.10, 0.25, 0.50, 0.75). Right: four pairs "lower threshold -> higher threshold". Fidelity = the C* at which a run may be stopped. LF = kLa at the lower threshold, known at all runs; HF = kLa at the higher threshold, known at n runs (x axis). y = error on the other runs (RMS of ln ratio). 79 L8 runs, 15 random draws per n.
- Lines: mfbml (blue); GPR with a linear mean prior in (rpm, angle) (black); GPR whose prior mean is the lower-threshold value plus one offset (grey dashed); the lower threshold times one constant (green dotted).
- Result: the lower threshold carries the information. One constant (or the GPR with that prior) reaches 0.20-0.27 from n = 4 runs on, for 0.25 -> 0.50, 0.10 -> 0.50 and 0.50 -> 0.75; 0.25 -> 0.75 stays at 0.4. mfbml is equal or slightly worse at small n. A GPR on the higher threshold alone needs more than 32 runs to reach 0.32-0.35. The run-to-run floor of kLa at L8 is 0.10-0.19 per run.
- Practical: stopping a run at C* = 0.25 instead of 0.50 saves about 4x in time (L8, 25 rpm: 42 s vs 173 s) and predicts kLa(0.50) to about 20 %.
- Correction made before the meeting: the first version used a quadratic-mean GPR (6 coefficients), which is degenerate at n < 6; its "mfbml wins" verdict is recorded in the diary as a flaw of the baseline.
- The mixing-time version (chi up to 0.5 -> time to 0.95) was tested earlier (diary 2026-10-08): about 34 % error, the early curve does not carry the late rate.

## 04  PROPER mfbml with a converged cycle repeated as the low fidelity (PI idea 2)
Figure NOT ready on 2026-10-21. What exists:
- `10_one_period_repeated.png`: one case (L8, 32.5 rpm, 7 deg): looping the first recorded flow period reproduces kLa(0.10) and mixing time (0.50) of the coupled run to 0.6 % and 0.0 %; kLa(0.25) is off by 15-51 % (more periods help). Cost at L8: 4 min vs 2 h.
- Running: a pipeline that records ONE period after a restart from each L8 grid run and repeats it 22 cycles, for up to 40 (rpm, angle) conditions (about 1 core-hour each). When it ends, the mfbml (LF = repeated period, HF = the normal L8 run) is built exactly like card 03 (learning curve, same baselines).
- The idea in mfbml terms: LF is cheap (settle the flow, record one period, repeat it); HF is the coupled run.

## 05  A hydrodynamics-informed law for kLa (new lead)
`05_kla_law_dissipation.png`. Left: kLa(0.10) against the dissipation rate (EDR) of the same flow at L9 (green) and L10 (pink), straight lines in the log-log plot. Middle: the 50 % mixing time against EDR at L8 (blue) and L9 (green). Right: error for kLa at L10 by method.
- kLa ~ EDR^0.6 at L10 (13 runs, correlation 0.97), ~ EDR^0.5 at L9 (10 runs, 0.93); no relation at L8 (0.11). The interface area per liquid area is the same at all levels (28.6-28.9 1/m), so the law is a law for the transfer velocity kL. Small-eddy theory gives 0.25: the exponent is about twice that; WE DO NOT YET EXPLAIN IT.
- Beyond the rpm and angle trend: partial correlation 0.81 at L10 (13 runs).
- Error at L10 (leave-one-out): EDR law 0.16, L9 times ratio 0.17, chain 0.19, rpm and angle regression 0.21, L10 runs only 0.23.
- BLIND TEST QUEUED: predictions for the five 5.5 deg L10 kLa runs are frozen in `experiments/multifidelity/physics_law_frozen_theta55.json` (6.1, 7.3, 9.7, 12.0, 14.2 1/h at 17.5, 20, 22.5, 27.5, 30 rpm); the runs are in the queue.
- Mixing time: slope -0.62 (L8) and -0.66 (L9) on ln EDR, shared across levels, but a model with that slope and one L9 offset has error 0.35 (chain 0.22).

## 06  Warm starting: does it work, how much does it save?
`06_warm_start.png`. A: change of condition at the same level, cycles to settle (arrows: never settled in the recorded cycles; dotted line: the 80 cycles of a cold start). B: change of level, difference to a cold L10 run at 32.5 rpm. C: is the training flow settled? difference between the training window and the end of the continuation, for 13 L10 runs. D: L10 cycles needed before the release (cold 80, from L8 12, from L9 6).
- Resume at the same condition: validated at L8 (continued run = uninterrupted run to 0.00 % in kLa(0.10), kLa(0.25), 50 % mixing time).
- Cross-parameter transitions are the least reliable (many never settled in 60 cycles).

## 07  Cost against error for the hydrodynamic quantities
`07_cost_vs_error.png`. Error on unseen L10 runs against the cost of the training data (L10 run 935 core-h, L9 run 74 core-h). Pink: L10 runs only (GPR); green: L9 times a fitted ratio; black: three-level chain. 13 L10 runs, 15 random subsets per size.
- L9 at every design point plus 3-4 L10 runs gives 4-6 % at about 4,000 core-hours; L10 runs alone give 13-190 % at that cost and need all 13 runs to reach 13-21 %.
- The chain is NOT more accurate than the fitted ratio (it is better for tau_95 only).

## 08  An experiment, not a claim: same flow, finer scalar mesh
`08_fixed_flow_experiment.png`. One recorded L8 flow (32.5 rpm, 7 deg, 22 cycles after the release), oxygen and tracer solved on the L8, L9 and L10 mesh (exact record and replay, flux made divergence-free on the finer mesh). Values relative to the L8 mesh.
- 50 % mixing time 11.8 / 19.0 / 24.7 s; kLa(0.25) 28.6 / 38.7 / 51.6 1/h; kLa(0.10) 46.9 / 38.3 / 50.2 1/h.
- This states only what happens to the numbers; it does not say which value is physical. Caveats: one flow and one rpm; the finer-mesh flux is the coarse flux interpolated (about a 5 % systematic); the flow itself is the L8 flow.

## 09  Runs 0.2 rpm apart
`09_run_pairs.png`. Difference between two runs 0.2 rpm apart, at L8 and L9 (same binary and protocol).
- At 30 rpm L9 pairs differ by 2-3 % and L8 pairs by 12-37 %; at 22.5 rpm both levels differ a lot (L9: 28 % in kLa, 23 % in the 50 % mixing time, at least 68 % in the 95 % mixing time).
- At 22.4/22.6 rpm the dissipation differs by 30 % and the late mixing rate by 3.6x between the two L9 runs: a narrow flow regime near 22.5 rpm that appears at L9 and L10 and not at L8 (card 01, EDR row).

## 10  One recorded period repeated (see card 04).

## 11  Is it a bug? (the postprocessing question)
`11_is_it_a_bug.png`. Left: the same run, run again (L9 at three rpm; reruns at L8): 0.01-1 % difference. Middle: Kim's MATLAB code on our fields: his default mixing time is 8 % (L8) and 3 % (L9) slower because he counts air-filled surface cells; with those excluded it is within 2 %; his kLa fits agree within 1.6 % at C* = 0.25 and 0.50 (14 % at 0.10, his window includes the release). Right: a different estimator keeps the scatter (+3 %, +14 %). The jumps between runs 0.2 rpm apart are in the simulated flow and Kim's code sees them too.

## Questions for the PIs
1. Which kLa definition do the experiments use (C* range, a fit over which time window)? Would they accept kLa at C* = 0.10 and the 50 % mixing time as targets?
2. Is there any measured kLa or mixing time of this bag (even rough) to anchor the level dependence?
3. The exponent in kL ~ EDR^0.5-0.6: is there a mechanism that gives about 0.5 (surface renewal with numerical diffusion, wave breaking)?
4. What would they accept as a convergence criterion for a scalar quantity in this flow?

## Feedback after the meeting
(to be filled in)
