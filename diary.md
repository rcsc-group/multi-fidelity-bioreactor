# Lab diary (multi-fidelity-bioreactor)

How to read: CURRENT STATE is authoritative and self-consistent. Dated entries follow in order, from 2026-10-08.
Entries before 2026-10-08 are kept verbatim in diary_archive_to_2026-10-07.md (split on 2026-10-10; full copy of the
pre-split file at /oscar/scratch/eaguerov/tmp/diary_backup_20261010.md). The archive is NOT re-audited line by line:
where it conflicts with CURRENT STATE, CURRENT STATE wins, and the known conflicts are listed under CORRECTIONS.

## CURRENT STATE (updated 2026-10-10)

Goal and constraints
- APS DFD 2026 talk (Orlando, ~Nov 22-24). Must keep a form of mfbml (Bessa is the PI). Pivot 2026-10-08: vanilla
  mfbml (Yi et al. 2407.15110, KRR-LR-GPR) with the finest level as target; gcbml shelved.
- Working rules: pre-register, state each test's power, run my own test before any adversarial review; Sonnet agents
  implement (<= 2 in parallel), I decide and read the literature in full.

Model and data facts
- 2-D Basilisk two-phase rocking bag; levels L6-L10 (NN = 2^L); theta 7 deg unless stated.
- Kim's postprocessing (BioReactor3D dev/postprocessing, read 2026-10-10): chi = 1 - sigma^2/sigma0^2 over liquid
  cells at cell scale, optional block projection; kLa from an exp fit of 1 - C/C* from release.
- Tracer D = 0.44e-9 m^2/s (Pe ~ 1e7, Batchelor scale ~50 um [superseded: 4-10 um, diag_batchelor.py] vs L9 cell 0.49 mm); oxygen D = 1.9e-9. Neither scalar
  is resolved at any level.
- Measured run costs: L8 kmix protocol (80-cycle spin-up + 150 cycles, 16 ranks) ~57 core-h (kmix_l8_rpm20 12,761 s);
  L9 fig9 protocol (32 ranks, prod binary) ~580 core-h (fig9_l9_rpm32.5 18.2 h); for Kim's case: L8 ~16, L9 ~410,
  L10 ~3600 core-h. L10 6-cycle warm start: ~300-550 core-h (estimate from cost_model rates, not yet measured).

Results in force
- tau_mean, tau_95, EDR: smooth in rpm at L8 (T5). Level changes: tau_mean grows 8 -> 18% per level from L6 to L9,
  0.4-9% from L9 to L10 (9 rpm); grid convergence is NOT established (order p not identified with 3 levels).
- mfbml L8 -> L10 (T2b, full-Bayes residual GP, 3 L10 points): tau_mean 9.7%, tau_95 5.7%, EDR ~20%; 95% coverage
  97%. Baselines: L8 x ratio 17-32%, L10-only GP (integrated, fair) 25-27%. With L9 as LF the ratio is as good.
  T2b missed its pre-registered RMSE rule by a hair (1.22 vs 1.2).
- dtmix and kLa: rough in rpm at L8 (T5: +-0.1 rpm moves dtmix_0.95 by 0.18-0.35 in log); do not converge L6-L10
  (dtmix_0.95 at 32.5 rpm 68/128/191 s for L8/L9/L10; kLa_1T_10 47/29/22 h^-1). mfbml fails on them (T7, L8 -> L9).
- At fixed rpm the flow can be phase-locked: release at cycle 80/82/85 gives identical dtmix at 32.5 rpm (L6-L8),
  but CV up to 17% at 17.5 rpm (L6). Release after 33 or 50 cycles changes dtmix by 12-23% (spin-up matters).
- Dead ends: truncated-chi time-horizon fidelity (test B); EDR as a kLa proxy (Spearman 0.36-0.90); T4 retro BO
  inconclusive (9-point pool).
- Prior art (lit v2 + reviews): no grid-level MF surrogate/BO on a rocking bioreactor found (12 S2 queries failed:
  not beyond reasonable doubt). Known neighbours: Ulam/transfer-operator mixing from CFD (2603.13996, 2112.11497),
  Lagrangian scalar models (2509.25030), non-reproducible LES mixing times (Haringa 2503.20622), wave-bag kLa closure
  with per-bag fitted constant (Piontek, PMC12868168), MF-OT-ROM for two-phase flows (2603.04232).

Open tests (pre-registered in the 2026-10-10 entries; power stated there)
- P (perturbation-size scan, L8): is the dtmix roughness noise or smooth-steep structure? Jobs 7277523-30.
- L (fixed-scale chi, L8/L9 at 32.5/37.5 rpm): is the level dependence of dtmix a measurement-scale artefact?
  Jobs 7276682-5.
- R (exact reruns): determinism check only. Jobs 7276608-15 (also the flow recordings for the estimator).
- L8 -> L10 warm start (x8_l10, 3 runs): jobs 7276284-6; reference xlevel L10 and, at 32.5 rpm, cold l10c.
- T6 rpm x theta L8 grid: 49/50 done or running (11 resubmitted 7276248-58).
- Q (periodicity per rpm, 10 short L8 recordings flow_l8_*): jobs 7277961-70. (Test M, a dimensionality-reduction
  test on the same runs, was withdrawn before any data.)
- Question chain for the talk: Q1 fixed-scale mixing (test L), Q2 noise vs structure (tests P, Q), Q3 kLa convergence.
- Plan: docs/dfd_plan.md v2 is STALE on test R and test L rules; v3 after P and L.

CORRECTIONS (conflicting statements in the archive or earlier entries, resolved)
1. "the dtmix rpm roughness is reproducible" (2026-10-08): withdrawn the same day; T5 shows it is effectively noise
   at L8 (cause open: test P).
2. "MF beats HF-only GP 5-10x" (tests A, T1, T3): the HF-only GP was a collapsed MLE fit; the fair figures are in T2b.
3. "Hydrodynamic QoIs converge (tau_mean) on L6-L10" (archive, 2026-10-05 heading): overstated. 2026-10-07 shows the
   increments grow through L9; only L9 -> L10 is small. Convergence is not established.
4. L10 xlevel runs described as "3 cycles after an rpm change" (archive 2026-10-05, 2026-10-07; flag
   "chained_rpm_change" in scripts/hydro_dataset_v2.py): inaccurate. Each L10 point restarts at the SAME rpm from an
   L9 state that ran 25 cycles after a 2.5-rpm change (L9 chain), then runs 6 L10 cycles; the tau window is the last 3
   of them. At 32.5 rpm it agrees with the cold L10 run within 0.5% (tau_mean) and 2% (tau_95) (archive 2026-10-05).
5. Costs "L8 ~16 core-h" refer to Kim's case (archive 2026-10-03), not to the kmix protocol (~57 core-h).
6. dfd_plan v1 (operator/closure) rejected by review; v2 partly superseded by the 2026-10-10 power amendments.

## 2026-10-08 — Pivot after advisor meeting: vanilla mfbml, HF target. Two pre-registered tests.
- Context: APS DFD 2026-11-25. Advisors: gcbml too experimental for DFD. Free choice of canonical case IF prior art in bioreactors is absent beyond reasonable doubt (lit agent dispatched).
- TEST A (zero compute), PRE-REGISTERED before any run. Script scripts/test_mf_l10_holdout.py.
  - Data hydro_dataset_v2.json. QoIs tau_mean_t, tau_95_t. LF = L8 (variant L9), all 10 rpm. HF = L10 (9 rpm, 17.5-37.5).
  - Designs: n_HF = 3 (17.5, 37.5 + each interior rpm: 7 designs) and n_HF = 4 (endpoints + 2 interior: 21 designs). Score on the held-out L10 rpm.
  - Methods: MF = Yi KRR-LR-GPR (scripts/mfbml_local, noise estimated by the method); HF-only GP (same GPR, constant basis); LF x ratio (mean HF/LF of training); LF alone.
  - Metric: relative RMSE on held-out points, median over designs. Secondary (not a pass criterion): 95% coverage.
  - PASS = MF median rel RMSE below BOTH HF-only GP and LF x ratio in >= 3 of the 4 (QoI x LF level) cases at n_HF = 3.
  - Caveat stated now: L10 runs are chained (3 cycles after an rpm change) except 32.5 rpm.
- TEST B (zero compute), PRE-REGISTERED. Truncated chi(t) -> dtmix_0.95, Bayesian, L9 rpm sweep (10 runs).
  - Observed: chi(t) up to chi = 0.50 (the curve). Model: t_0.95 = t_0.50 + ln(10)/(r lam_e), lam_e = early decay rate fitted on chi in [0.25, 0.50] (posterior from the linear fit), r = late/early rate ratio.
  - Priors on log r: P1 physics, N(0, 0.6^2) (r 95% in [0.3, 3.2]); P2 cheap-level population, N(m, s^2) fitted on full L6-L8 curves (all rpm); P3 grid-MF, centred on log r of L8 at the same rpm, sd = spread of log r(L8) - log r(L7) over rpm.
  - PASS (per prior) = 95% coverage >= 8/10 AND median |rel err| of the median < 20% AND median relative 95% width (hi-lo)/median <= 0.6.
  - Reference: deterministic single exponential (pilot_truncated_chi.py) median |err| 34% at L9.
  - Cost caveat stated now: stopping at chi 0.50 saves only 1.5-2.2x at L9 (spin-up, diary 2026-10-03).
- DATA QUALITY, dtmix (scripts/diag_dtmix_quality.py -> dtmix_quality.json). Roughness in rpm r = log y - mean(log neighbours):
  - median |r|: L6 0.19-0.32, L7 0.16-0.18, L8 0.12-0.17, L9 0.19-0.27, Kim 0.05-0.09 (chi 0.50/0.75/0.95).
  - replicate sd of log y: 0.002-0.03 at L7/L8 and at L6 for 25-32.5 rpm; 0.06-0.17 at L6/L7 17.5 rpm.
  - So the rpm wiggles at L8/L9 are 5-50x the run noise: deterministic, not noise. Their sign often differs between levels at the same rpm (chi 0.95, 30 rpm: L8 +0.40, L9 -0.27); Kim shares only the 22.5 rpm dip (-0.28 vs L9 -0.57). Reading: grid-dependent, non-transferring structure (unconverged data), not bad runs.
  - L9 replicate rep_l9_rpm25_rel82 never ran (params.json only). No L9 noise estimate.
- TEST B RESULT (scripts/test_truncated_chi_bayes.py), L9, 9 rpm with measured t_95, 411-790 curve samples up to chi 0.50 per run:
  - log r of the full curves: L6-L8 population mean +0.01, sd 0.40; L9 values -1.0..+0.57.
  - P1: coverage 9/9, median |err| 34.9%, median width 2.25 -> FAIL (width).
  - P2: 7/9, 33.9%, 1.36 -> FAIL. P3: 6/9, 33.5%, 1.42 -> FAIL. Single exponential: 34.8%, 0/9.
  - Reading: the curve up to chi 0.50 does not contain the late decay rate (consistent with diary 2026-10-02). A calibrated prior gives an honest band of a factor ~2-3; L8 at the same rpm does not narrow it (late/early ratio does not transfer from L8 to L9).
- TEST A RESULT (scripts/test_mf_l10_holdout.py, Sonnet agent, reviewed): PASS, 3 of 4 cases at n_HF = 3.
  - Median held-out rel. RMSE, n_HF = 3: tau_mean L8->L10: MF 8.0% vs ratio 17.3% vs HF-only GP 53.9%; tau_95 L8: 5.5 vs 18.8 vs 52.6; tau_95 L9: 6.7 vs 9.2 vs 52.6; tau_mean L9: MF 7.4 LOSES to ratio 4.9 (L9 already within ~5% of L10).
  - 95% coverage of the MF band: 0-57% (overconfident). HF-only GP 69-85%.
  - Cause (read in code + Yi 2407.15110 App. B eq 23): the fit maximises eq 23 with unit signal variance on standardised y (code matches the paper; the BACKLOG "likelihood bug" item is NOT a bug). But predict() scales the epistemic variance by the plug-in sigma2 = r'K^-1 r / n. With n_HF = 3 and a 2-term linear transfer, r has 1 degree of freedom, so sigma2 -> ~0 and the band collapses. Not tested further; a fix needs a prior on the residual amplitude or more HF points.

## 2026-10-08 (cont.) — User decision: show our sharpest tools; tests 1-6 approved. PRE-REGISTRATION.
- Abstract (DFD 2026, user-pasted): Basilisk two-phase NS + multi-fidelity BO for rocking bioreactors, kLa vs shear. User: DFD is relaxed about abstracts; do not drift; show the sharpest tools.
- Data note: clean cold L10 at 32.5 rpm (l10c_rpm32.5_seg1) dtmix_0.95 = 190.7 s vs Kim 184 s vs L9 128 s. Unchecked.
- Concession: "the dtmix rpm roughness is reproducible" was too strong. Replicates exist only at 17.5/25/32.5 rpm, L6-L8, and vary only the release cycle on a settled (often phase-locked) flow. They do not test smoothness in rpm or sensitivity to the state.
- Binary f1c11e0 (lean, the kmix_l8 binary) copied to /oscar/data/dharri15/eaguerov/bin/ (scratch purge risk; file dated 2026-09-21).
- T1 (zero compute): test A protocol unchanged, QoI = EDR (mean energy dissipation rate, same window as hydro_dataset_v2: last 3 cycles before release), LF L8/L9 -> HF L10 (9 rpm). PASS as test A: MF median rel RMSE below both HF-only GP and LF x ratio in >= 1 of the 2 LF cases at n_HF = 3 and not worse than 1.5x the best baseline in the other.
- T2 (zero compute): prior on the residual amplitude in KRR-LR-GPR (fully Bayesian over the GP amplitude, replacing the plug-in sigma2 that collapses with 1 residual dof). PASS = pooled 95% coverage >= 85% on test A (tau, all 4 cases, n_HF = 3) AND median rel RMSE <= 1.2x the test-A value in every case.
- T3 (zero compute): learning curve, n_HF = 2..6 (all designs with both endpoints when feasible, else 20 random designs), tau_mean/tau_95/EDR, MF vs HF-only GP vs LF x ratio. Descriptive; no pass criterion. Decides L9 vs L10 as HF for T6.
- T4 (zero compute): retrospective pool-based BO on the 10-rpm grid. Objective: maximise EDR at L10 subject to tau_95(L10) <= the median L10 tau_95 (a stated proxy for "mixing vs shear"; kLa at L10 has 3 points only). MF-BO (L8 known at all rpm, L10 acquired one at a time from 2 initial endpoints) vs single-fidelity BO (L10 only, same start). Metric: L10 evaluations to find the true best feasible rpm, over all start pairs/seeds. Descriptive.
- T5 (compute, ~100-300 core-h): roughness test at L8, same binary/protocol as kmix_l8. Dense rpm 20.625..31.875 step 0.625 (15 new runs; with existing 20..32.5 at 2.5) + rpm perturbation 22.4/22.6/29.9/30.1 (4 runs). Reading fixed now:
  (a) REAL structure if the dense dtmix_0.95(rpm) curve has |second difference| at 0.625 spacing <= 3x the replicate log-sd (~0.02) except across <= 2 narrow features;
  (b) EFFECTIVELY NOISY if |log y(rpm +- 0.1) - log y(rpm)| or dense-neighbour jumps are >= 0.10 (comparable to the 2.5-rpm jumps 0.12-0.27).
- T6 (compute): second design parameter, rpm x angle. LF = L8 dense: theta in {2, 4, 5.5, 9} deg x 10 rpm (40 new runs; theta 7 exists). Same protocol (80 cycles spin-up + 150 post). HF (L9 or L10, 6-8 points) chosen after T3.
- SUBMITTED T5 + T6 (scripts/submit_l8_t5_t6.py): 59 L8 jobs 7156006-7156065 (16 ranks, 8 h cap). --check vs kmix_l8_rpm25/32.5: params identical except the binary path (cmp: same file). Run ids t5_l8_rpm<r>_th7, t6_l8_rpm<r>_th<theta>.
- Zero-compute T1-T4 dispatched to a Sonnet agent (spec = the pre-registration above). Prior-art v2 (Semantic Scholar free tier, raw JSON kept) + broadened canonical-case search dispatched in parallel.
- T1-T4 RESULTS (Sonnet agent; reviewed: figures opened, unit tests 3/3 pass, old dataset keys bit-identical).
  - T1 EDR (edr_mean_t, shear_stress.dat col 9, W/m^3, same window): LF L8: MF 20.0% vs ratio 32.4% vs HF-only 221.7% (win). LF L9: MF 22.3% vs ratio 8.3% (L9 alone 6.4%) -> 2.7x worse -> FAIL.
  - T2 Bayesian amplitude (Student-t, IG(2, 0.1)): point predictions identical; pooled coverage 11.3% (19/168) -> FAIL. Root cause found by the agent and confirmed in the figures: the eq-23 MLE puts the noise ratio at its 1e-5 bound and the kernel parameter at 0.01 (near-interpolation), so the kriging variance factor is ~1e-9. The amplitude is not the limiting factor; the hyperparameter MLE at n_HF = 3 is.
  - The same MLE pathology breaks the HF-only GP baseline: in test_t1_edr.png it is a flat line with spikes at the training points (length scale collapsed). So "MF beats HF-only GP by 5-10x" in tests A/T1/T3 is partly a broken baseline. LF x ratio is the honest baseline.
  - T3 learning curve: with LF L8, MF beats ratio at every n_HF = 2..6 for tau_mean (8-12% vs 16-19%), tau_95 (6-8% vs 17-19%), EDR (20-23% vs 24-33%). With LF L9, ratio wins for tau_mean and EDR (L9 ~ L10); MF wins for tau_95 up to n = 5. MF coverage rises with n (up to 88% at n = 6) but is < 60% at n <= 4.
  - Reading: MF pays when the LF is cheap and its discrepancy is not a constant factor (L8); when LF ~ HF (L9), a ratio suffices.
  - T4 retro BO (9-rpm pool, 5 feasible, best 22.5 rpm): L10 evaluations to the best, median: MF-Bayes L8 5, HF-only GP 4, MF-Bayes L9 5, random 6.3. No MF advantage. Inconclusive by design: a 9-point pool with 8/36 start pairs already containing the optimum cannot separate the arms.
- T2b PRE-REGISTERED (own follow-up of the T2 root cause; T2 approved by the user): fully Bayesian KRR-LR-GPR residual GP.
  - Keep: LF KRR, basis [1, f_l], rho flat, amplitude IG(2, 0.1) integrated analytically (as T2).
  - New: integrate the kernel length scale and the noise ratio on a 2-D grid (60 x 40) with priors: length scale in scaled input units log l ~ N(log 0.3, 0.75^2); noise ratio (noise sd / HF std) log-uniform on [1e-3, 0.3]. Marginal likelihood per grid node in closed form: |K|^-1/2 |F'K^-1F|^-1/2 (b0 + S/2)^-(a0 + (n-p)/2). Predictive = mixture of Student-t over the grid.
  - The SAME priors and machinery for the HF-only GP (constant basis), so the baseline is no longer a collapsed MLE fit.
  - PASS = on test A (tau_mean, tau_95; LF L8, L9; n_HF = 3): pooled 95% coverage >= 85% AND MF median rel RMSE <= 1.2x the test-A value in every case. Reported, not criteria: the fair HF-only GP, T1 (EDR), and T3 at n_HF = 2..6.
- T2b RESULT (scripts/test_t2b_fullbayes.py; module scripts/mfbml_local/krr_lr_gpr_fullbayes.py; 5 unit tests pass): FAIL by a hair on the pre-registered rule.
  - Coverage PASS: pooled MF 95% coverage 97.0% (163/168), was 11.3% (T2).
  - RMSE: ratio to test A 1.22 / 0.98 / 1.04 / 1.02 (tau_mean L8 / tau_mean L9 / tau_95 L8 / tau_95 L9); 1.22 > 1.2 -> FAIL (tau_mean L8: 9.7% vs 8.0%). Criterion not changed.
  - Fair HF-only GP (same priors, integrated): 25.7-27.3% at n_HF = 3, coverage 96-100%. MF beats it 3-5x; MF still beats LF x ratio with LF L8 (9.7 vs 17.3; 5.7 vs 18.8) and loses with LF L9 for tau_mean (7.3 vs 4.9).
  - Posteriors at n = 3 are broad: l ~ 0.4 [0.1, 1.2], eta barely moves from its prior.
  - Learning curve: MF coverage falls with n (tau_95 L8 100% at n = 2 -> 79% at n = 6); fair HF-only GP catches up by n = 6 for tau_95 (8.3% vs MF 8.4%).
  - Figure legend fixed (L8 was missing); rerun gives identical numbers.
- User (2026-10-08): "tau mean, tau_95, EDR, mixing time, kLa. They all should be considered ... What about mixing time and kLa?" Correct: mfbml was never tested on dtmix or kLa.
  - Why: L10 has no usable dtmix/kLa curve in rpm. Runs with tracer + oxygen at theta 7: l10c_rpm32.5_seg1 (cold, 80-cycle spin-up, the kmix protocol), fig11_l10_rpm35/37.5 (warm start, release after 5 cycles; rpm 35 ended before chi = 0.75), fig8_hist_l10b / fig9_l10* at 32.5 (other protocols). So <= 3 rpm, mixed protocols: no held-out test is possible at L10.
  - L8 (kmix_l8_rpm*) and L9 (fig9_l9_rpm*) both have all 10 rpm with the same protocol (80-cycle spin-up, release, then 150 cycles).
- T7 PRE-REGISTERED (zero compute): test A protocol with L8 -> L9 as a proxy for L8 -> L10, on dtmix and kLa.
  - QoIs: dtmix_0.50, dtmix_0.75, dtmix_0.95, kLa_1T_10, kLa_1T_25, kLa_1T_50 (whole-period kLa, as Figs 11/12).
  - LF = kmix_l8_rpm* (10 rpm). HF = fig9_l9_rpm* (10 rpm); dtmix_0.95 at 37.5 from fig9_l9_rpm37.5_ext1; kLa_1T_50 at 30 rpm is nan -> that rpm is dropped for that QoI only.
  - Fit on log y (positive QoIs that span up to 60x); score on exp of the predictive median. Same metric otherwise (rel RMSE over held-out L9 rpm; 95% coverage from mixture quantiles).
  - Model: T2b full-Bayes MF (FullBayesKRRLRGP, basis linear) and the same machinery with basis ordinary (HF-only GP); LF x ratio (ratio = mean of yh/yl over training rpm) as the honest baseline. n_HF = 3, designs = both endpoints (15, 37.5) + one interior.
  - PASS per QoI: MF median rel RMSE < LF x ratio AND < HF-only GP. Overall pass = >= 4 of 6 QoIs pass AND pooled MF coverage >= 85%.
  - Reported, not criteria: n_HF = 2..6 learning curve; the raw L8 vs L9 values (does L8 even rank the rpm like L9?), Spearman rank correlation L8 vs L9 per QoI.
  - Caveat fixed now: a pass says mfbml can carry dtmix/kLa from L8 to L9, not to L10. L8/L9/L10 dtmix_0.95 at 32.5 rpm = 68/128/191 s: not converged, so L9 is not a stand-in for L10 in value.
- T7 RESULT (scripts/test_t7_mix_kla_l8_l9.py, Sonnet agent; reviewed: figure opened, markers fixed to the figstyle chi_0.50/0.75 markers (my brief said chi_0.95 for all), rerun json identical except runtime): FAIL.
  - Per QoI at n_HF = 3 (MF / HF-only / L8 x ratio median rel RMSE; MF coverage): dtmix_0.50 31.7/42.6/30.0, 77% FAIL; dtmix_0.75 34.7/41.4/34.0, 71% FAIL; dtmix_0.95 45.0/46.1/61.5, 82% pass; kLa_1T_10 22.9/42.4/26.2, 89% pass; kLa_1T_25 45.0/78.7/72.5, 55% pass; kLa_1T_50 180/218/348, 55% pass. 4/6 pass, but pooled MF coverage 72.4% < 85% -> FAIL.
  - Spearman L8 vs L9 over rpm: dtmix 0.79/0.81/0.56; kLa_1T 0.53/0.06/0.22. For kLa at C* = 25/50% L8 does not even rank the rpm like L9: no transfer method can work on these QoIs as measured.
  - Reading: on dtmix the L8 shape is useful at chi 0.50/0.75 (ratio ~ MF, both beat HF-only), but nothing reaches < 30%. The "passes" on kLa 25/50 are wins among methods that all fail (45-348%). kLa_1T_10 is the only clean MF result (23% -> 14% from n = 3 to 6).
  - Consistent with the roughness question (T5): if dtmix/kLa(rpm) at one level is jagged at the 20-50% scale, a smooth transfer cannot fit it. T5 decides whether the jaggedness is signal or noise.
- Prior-art v2 (Semantic Scholar free tier + citation chasing; docs/lit_review_mf_bioreactor_v2.md; raw JSON in /oscar/scratch/eaguerov/tmp/ss_raw/): no prior art for grid-level MF surrogate/BO on a rocking/shaken/single-use cell-culture bioreactor with kLa/shear objectives. 23 queries with data, 1,865 papers, 8 seeds chased. 12 queries failed with 429 (incl. "multi-fidelity bioreactor", "multi-fidelity sloshing", all single-use variants): NOT beyond reasonable doubt until those run with the user's key. Eskandari 2311.05776 read in full: Hartmann6D test + plan only, no bioreactor result. Canonical ranking: SPH dam-break 4TU (downloads under maintenance), Hysing rising bubble (1.7 MB zip, 5 grid levels, published QoIs), SPHERIC Test 10 sloshing. Setinek 2511.01830: wall shear shows no LF->HF transfer (relevant to our shear objective).

## 2026-10-10 — T6 resubmission; EDR as a kLa proxy (pre-registered)
- 11 T6 jobs were "CANCELLED by 0" (uid 0 = system) at start on 2026-10-08 16:10-17:12, not by me. Resubmitted with the new --only flag of submit_l8_t5_t6.py: jobs 7276248-7276258. T5 (19) all complete; T6 29/40 complete.
- Binary BioReactor-mpi-xlevel-chain (scratch, dated 2026-09-15, near the ~30-day purge) copied to /oscar/data/dharri15/eaguerov/bin/ (cmp identical).
- User: "i want to push harder the mfbml vanilla idea, because I NEED to present something to dfd".
- PRE-REGISTERED (zero compute): EDR (= power per volume) as a proxy for kLa, so the talk's objective can use QoIs that transfer (tau, EDR). PASS if Spearman(EDR, kLa_1T_c) >= 0.8 at both L8 and L9 for c = 10 and 25 (10 rpm, theta 7). Script scripts/diag_edr_kla_proxy.py.
- EDR-kLa proxy RESULT: FAIL. Spearman L8: 0.66/0.61/0.70 (kLa_1T_10/25/50), L9: 0.90/0.36/0.37. Only L9 kLa_1T_10 tracks EDR (log-log slope 0.50). EDR is not a kLa proxy on our data; the talk objective stays on tau and EDR as themselves.
- Plan to push vanilla mfbml for DFD (user request): hero result = 2-D design map (rpm x theta) of tau_95 and EDR at L10 accuracy, LF = 50 L8 runs (T6 + theta 7), HF = 6-8 L10 runs, plus MF-BO picking the HF points. Enabler: L10 warm-started from the settled L8 state (~6 L10 cycles, ~300-550 core-h per point instead of ~3600 cold).
- PRE-REGISTERED (compute, 3 L10 runs, ~1,500 core-h): L8 -> L10 two-level warm start (scripts/submit_l10_from_l8.py; same mechanics and N_CYCLES = 6 as submit_l10_fig13a_xlevel.py, which starts from L9). theta 7, rpm 17.5/25/32.5, source kmix_l8_rpm<r>. Reference: l10_fig13a_xlevel_rpm<r> (L10 from L9). QoIs tau_mean_t, tau_95_t, edr_mean_t over the hydro_dataset_v2 window (last 3 cycles before release).
  PASS if, for every QoI and rpm, |x8 - xlevel| / xlevel <= max(5%, 2 x the cycle-to-cycle relative sd of the xlevel run over the same 3 cycles). Also reported: per-cycle tau/EDR after the restart (does it settle within 3 cycles?).
  If FAIL: use L8 -> L9 -> L10 chains for HF points (cost about 2x) or more L10 cycles.
- T5 RESULT (scripts/test_t5_roughness.py, Sonnet agent; reviewed: figure opened, the rpm +- 0.1 markers changed from filled to hollow (filled = extremum in figstyle), rerun json identical): EFFECTIVELY NOISY ((b) holds, (a) fails).
  - dtmix_0.95, L8, theta 7: 17/19 interior |second differences| > 3 x replicate sd (0.0555), in 3 features (> 2). rpm +- 0.1 changes log y by 0.18-0.35 (all 4 pairs >= 0.10); 12/20 neighbour jumps at 0.625 rpm >= 0.10.
  - Same metrics, reported only: dtmix_0.50/0.75 and kLa_1T_10/25/50 are also rough (max jumps 0.28-1.06; +-0.1 rpm pairs up to 0.49). tau_mean_t and tau_95_t are smooth (no jump >= 0.10, max 0.09/0.08); edr_mean_t nearly (2 of 24 >= 0.10).
  - Reading: a single-run dtmix or kLa is not a smooth function of rpm at L8. A 0.1-rpm change moves it by 20-40%, 10-20x the release-cycle replicate sd (those replicates sample a phase-locked flow and understate the spread). tau and EDR are smooth. This explains T7 vs T2b: mfbml transfers tau/EDR and fails on dtmix/kLa because of the QoIs, not the method.
  - My concession to the user stands and is now measured: the roughness is not reproducible structure.
  - Consequence for MF on dtmix/kLa: treat a run as a noisy sample (noise sd ~0.2 in log at L8) and average an ensemble (e.g. rpm +- small perturbations) before any transfer. Not for DFD unless cheap.
- User: plan, not summary; success = an adversarial Sonnet agent approves an above-average DFD contribution at L10 for tau_95, EDR, dtmix, kLa. Draft plan docs/dfd_plan.md (C1 noise finding, C2 transfer-operator mixing from recorded flow, C3 kLa closure from hydrodynamics, C4 mfbml + MF-BO on rpm x theta; tests F1-F7). Adversarial Sonnet review dispatched.
- Adversarial review of docs/dfd_plan.md (Sonnet): below the bar as written. Main points, accepted:
  - Prior art: Ulam period map from CFD for stirred-tank mixing times (Kluenker et al. 2603.13996; 2112.11497); wave-bioreactor kLa closure with VOF interface and fitted constant c differing ~12x between bags (Piontek et al., PMC12868168). C2/C3 are not new methods; new only in setting.
  - F1/F3/F5(a) nearly unfalsifiable; F4 passable by tuning the Ulam cell size; F7 saving comes from the warm start, not the operator.
  - L10 is not a converged reference for dtmix (68/128/191 s at L8/L9/L10) or kLa (47/29/22).
  - Sharpest open question the data already pose: is the rpm roughness of dtmix in the flow or in the tracer? Same-run operator test + diffusion control.
  - Plan to be revised and re-reviewed before large compute.
- My reading added: dtmix rises and kLa falls with refinement, both consistent with scalar numerical diffusion falling with the grid (tracer D = 0.44e-9, Sc ~ 2300; oxygen Sc ~ 500, both unresolved). Hypothesis H_nd: the scalar QoIs are limited by the scalar's numerical diffusion, not by the flow (tau/EDR change much less L8 -> L10).
- PRE-REGISTERED (compute ~450 core-h): 8 L8 flow-recording runs (scripts/submit_l8_record.py): kmix/T5 protocol and lean binary f1c11e0, + 100 snapshots/period over cycles 78-90. rpm 22.4, 22.5, 22.6, 25, 29.9, 30, 30.1, 32.5, theta 7.
  - Reproducibility check (reported): each run repeats an existing run exactly except the snapshot params (no VIDEOS in the lean build, dt_video drives only snapshots). If dtmix_0.95 differs from the original by > 0.02 in log, single runs are not reproducible (round-off chaos) and that is itself the noise source.
  - Same-run operator test criteria: fixed in the revised plan before the operator code exists.
- BACKLOG: bare `return;` in events movies_output / movies_output_tau (VIDEOS builds only, BioReactor.c ~2512, ~2623): CLAUDE.md hazard. Lean production binary unaffected (no VIDEOS).
- Round-2 adversarial review of plan v2 (Sonnet): v2 fixed the same-run test and dropped the kLa closure; remaining fatal confound: the code's chi is at the GRID scale, so a grid-dependent dtmix may be a measurement-scale artefact, not numerical diffusion; a fixed-l estimator would pass test D by definition. Cheapest fix: recompute the direct chi at a fixed scale l from tracer fields at L8/L9 (and L10) before building any estimator. Also: test R threshold too lax (fix below); prior art for non-reproducible CFD mixing times: Haringa et al. 2503.20622 (3-D LES). Lagrangian/particle mixing estimators: 2509.25030, 2603.13996.
- User (2026-10-10): adversarial review only AFTER my own experiments; I own decisions and the literature review (read in full); Sonnet implements, <= 2 in parallel. Hard constraint: the talk must keep a form of mfbml (Bessa is the PI). No experimental kLa available. GPU fine. Keep Kim's chi and add a fixed-scale version.
- Kim pushed his postprocessing to BioReactor3D (dev/postprocessing/, commit 00c0e7a). Read: bio_mix_elvis.m computes chi = 1 - sigma^2/sigma0^2 over liquid cells (variance, cell scale) and ALREADY has an optional projection (PROJECT = 1) onto N_projX x N_projY coarse boxes (1024x292 ... 16x4), i.e. a fixed-scale chi. bio_oxy.m: kLa from an exp1 fit of 1 - C/C* from release to the current time (global), plus a moving 4-point window fit; Henry coefficient 1/30.
- PRE-REGISTERED test L (scale artefact vs not), compute ~1,000 core-h: L8 and L9, theta 7, 32.5 and 37.5 rpm, lean binary f1c11e0, c snapshots 13/period over cycles 80-180 (scal_l<L>_rpm<r>_th7). chi_l(t) from c2 by Kim's PROJECT coarse-graining onto fixed boxes l = L0/16, L0/32, L0/64 (liquid cells, volume-fraction weighted), and at the native grid scale (check: native must reproduce results.json dtmix within 2%).
  Reading, per rpm and chi in {0.50, 0.75, 0.95}, d_grid = |log dtmix_L9 - log dtmix_L8| at native scale, d_l the same at scale l:
  SCALE ARTEFACT if d_l <= d_grid / 3 for l = L0/32 at >= 5 of the 6 (rpm, chi) cases; NOT an artefact if d_l >= 2 d_grid / 3 at >= 5 of 6; otherwise mixed (reported as such).
- Test R threshold REVISED before any rec_l8 result exists (disclosed): H_rep holds if the median |log rerun - log original| of dtmix_0.95 over the 8 reruns >= 0.10 (1/2 of the median +-0.1-rpm jump, ~0.20); reruns are reproducible if it is <= 0.05. Old threshold (0.02 at >= 2 of 8) dropped: it passes trivially for small effects.
- Lit (read in full, all sections + results): Khamlich, Tonicello, Pichi, Rozza, arXiv 2603.04232 (2026), MF-OT-ROM / PMF-OT-ROM. Field-level multi-fidelity for diffuse-interface two-phase flows: LF = coarser grid (128^2/256^2 -> 256^2/512^2; 40k -> 160k dof RT), residual field HF - LF interpolated in time by entropic-OT displacement interpolation between HF checkpoints of the SAME trajectory; parametric version first displacement-interpolates HF checkpoints across neighbouring parameter values. Tests: Rider-Kothe vortex (prescribed flow) and Rayleigh-Taylor (At 0.25-0.75). Key limits stated by the authors: "A key assumption ... the solution (or the residual field) varies smoothly with respect to the parameter mu. In problems where the parametric dependence exhibits bifurcations or sharp transitions in parameter space, the interpolation-based strategy may become less effective"; late-time RT breakup errors not reduced by more checkpoints; parametric version less accurate than time-only. Relevance: (i) field-level grid MF is the "spiritual mfbml" extension of scalar KRR-LR-GPR; (ii) it needs HF checkpoints on the same or neighbouring trajectories, i.e. HF runs; (iii) their own caveat matches our T5 finding (dtmix non-smooth in rpm). Interfacial area (their Fig 22: LF under-predicts) is a smooth hydrodynamic field QoI: candidate for kLa's "a" at L10. Unknown whether this is the paper Kim meant (other candidate seen: Cutforth & Mirjalili 2508.04084, autoencoders for 3-D interfacial multiphase flows; not read yet).

## 2026-10-10 (cont.) — POWER of the pending tests (user: "reflect on its power"); amendments made BEFORE any result
- Test R (exact reruns), power: WEAK, likely uninformative. Same binary, same rank count, same params: an MPI Basilisk run is probably bitwise deterministic, so |Delta| = 0 would say "reproducible" while saying nothing about sensitivity to tiny perturbations, which is the real question. R is kept only as a determinism check (does the T5 jump include run-to-run nondeterminism?). It is no longer the test of H_rep.
- NEW test P (replaces R as the H_rep test), PRE-REGISTERED: perturbation-size scan at L8, theta 7, base 22.5 and 30 rpm (T5 +-0.1 jumps 0.20-0.35 in log), offsets +1e-2, +1e-3, +1e-6 rpm (6 runs, scripts/submit_l8_perturb.py; with the T5 +-0.1 runs that gives 4 decades).
  Model: |Delta log dtmix_0.95| vs |Delta rpm|. Smooth-but-steep structure predicts slope ~1 in log-log (jump at 1e-3 ~ 1/100 of the jump at 0.1, i.e. < 0.005). Sensitive dependence (effectively noise) predicts no decay (jump at 1e-6 comparable to the jump at 0.1).
  Reading: NOISE if the median |Delta| over the two bases at 1e-6 rpm >= 0.10; SMOOTH if at 1e-3 rpm it is <= 0.02 AND the jumps fall by >= 10x per 2 decades; otherwise mixed.
  Power: the two hypotheses differ by a factor ~100 at 1e-3 rpm and ~1e5 at 1e-6; with 2 bases x 3 offsets the outcome is decisive unless the structure is a few sharp switches (then 1e-2 may straddle one; reported per offset). Caveat: if runs are bitwise deterministic, 1e-6 rpm (relative 4e-8) is still a real input change (omega_b differs at ~1e-8 relative).
  Also done for dtmix_0.50/0.75 and kLa_1T_10/25/50 (reported, not criteria).
- Test L power, amended: d_grid is small in some cases (37.5 rpm chi 0.50/0.75: 0.12/0.06 in log from existing runs), below the per-run spread (T5: ~0.2); a ratio d_l/d_grid is meaningless there. New rule: only (rpm, chi) cases with d_grid >= 0.20 in the NEW scal runs count (expected 3-4 cases: 32.5 rpm chi 0.50/0.75/0.95, 37.5 rpm chi 0.95). SCALE ARTEFACT if d_l <= d_grid/3 at l = L0/32 in all but at most one counted case; NOT if d_l >= 2 d_grid/3 in all but at most one; else mixed. Fewer than 2 counted cases -> test L inconclusive.
  Residual weakness (stated): L8 vs L9 is one realization each; if dtmix is noise (test P), d_grid itself contains ~0.2 of noise and test L cannot separate scale from noise. Test P is read first.
- L8 -> L10 warm-start validation power, amended: the reference l10_fig13a_xlevel (warm start from L9, 6 cycles) may carry the SAME transient bias as x8 (both 6 cycles after a coarse start); agreement would then not show that either is settled. Added reference: l10c_rpm32.5_seg1 (cold L10, 80 cycles) tau/EDR over the same window, at 32.5 rpm only. Reading at 32.5: x8 must also satisfy the same tolerance against l10c; if x8 matches xlevel but both miss l10c, the 6-cycle warm start is rejected for both.
- Lit (read: main text, Appendices B-D; Appendix A derivation skimmed): Lewis & Constante-Amores (UIUC), arXiv 2603.07297 (2026), exact coherent states in chaotic falling films. Two-phase dynamics reduced to the interface (2-D Kawahara/KS-type long-wave equation, Dedalus); ~2000 runs for a regime map; manifold coordinates by POD then IRMAE-WD (implicit rank-minimizing autoencoder, Zeng et al. 2024): intrinsic dimension d_M from the drop in latent singular values (d_M = 18 at L = 22, 76 at L = 30; POD decays slowly and "substantially overestimate[s]" d_M); neural ODE in manifold coordinates (DManD, Linot & Graham) used only to seed Newton-Krylov searches; 20 ECS found; leading Lyapunov exponent from twin trajectories (App. D). Limits stated: long-wave model only, not full NS/VOF; d_M scaling only at one delta.
  Relevance (my reading, untested): (i) our flow at fixed rpm can be periodic (phase-locked at 32.5 rpm) or not (CV 17% at 17.5 rpm, L6); a Lyapunov exponent / attractor dimension per rpm would say where single-run dtmix is a random variable; test P's 1e-6 rpm runs are a first twin-trajectory probe. (ii) Candidate idea, NOT pre-registered, contingent on tests P and S: mfbml on manifold coordinates of the periodic flow (LF L8 -> HF L10 per rpm), then QoIs (incl. a particle mixing estimate) computed from the reconstructed L10 flow. User: Kim's paper was from US authors near Chicago/Michigan and proposed a novel dimensionality-reduction method; this IRMAE/DManD line (UIUC / UW-Madison) is my best candidate; not confirmed.

## 2026-10-10 (cont.) — Dimensionality reduction: literature decision (provisional) and test M (pre-registered)
- User: do not hunt for Kim's specific paper; find the best-calibrated, most-adopted method as of Oct 2026 and act independently.
- Read in full (main text): Cutforth & Mirjalili 2508.04084 (KTH/Stanford): convolutional AE for 3-D interfacial fields; a moderately diffuse interface (tanh, 1-4 cells) beats sharp VOF-like and SDF inputs; a linear AE (= POD subspace) at the same compression is "substantially" worse; spectral bias loses small drops. Lewis & Constante-Amores 2603.07297: IRMAE-WD gives a calibrated intrinsic dimension where POD overestimates it.
- Seen via WebFetch summary only (NOT read by me; to read before citing): Saetta, Massa, Tognaccini, Iaccarino, Data-Centric Eng. 7 e30 (2026) "Multi-fidelity autoencoders: RANS-LES jet flow predictions": AE trained on 90 RANS, then 2 extra latent coords for 20 LES (soft penalty), RBF over design space, 50-member ensemble for model-form UQ, 20x cheaper database. Closest prior art for field-level MF with autoencoders. Also: manifold-alignment MF ROM (Proc. R. Soc. A 2022, POD); LaSDI family (LLNL).
- Provisional decision (to revisit after reading the Saetta paper in full): encoder = POD baseline vs convolutional AE on a moderately diffuse interface field, with IRMAE-WD to estimate the intrinsic dimension; transfer = KRR-LR-GPR (mfbml) on latent coordinates L8 -> L10 per design point; QoIs computed from decoded flows. This keeps mfbml central and adds Bayesian latent transfer (Saetta et al. use deterministic RBF + ensembles).
- PRE-REGISTERED test M (compute ~250 core-h): is the L8 periodic flow low-dimensional and smooth across rpm, so that a latent surrogate can interpolate it?
  Data: 10 short L8 recording runs flow_l8_rpm{15..37.5}_th7 (kmix protocol, stop at ~84.5 cycles, 100 snapshots/period over cycles 80-84) + the 8 rec_l8 runs.
  Method (Sonnet implements): per rpm, phase-align snapshots (phase from t mod T); encoders POD and CAE (+ IRMAE-WD dimension); leave-one-rpm-out for the interior rpm 17.5-35: predict the held-out rpm's flow over one period from latent coordinates interpolated in rpm (GP), decode, compare.
  Metrics on the held-out rpm: (a) EDR and interface length per phase from the decoded field vs the direct field; (b) relative L2 error of u, v.
  PASS (latent surrogate interpolates): median relative error of period-mean EDR <= 10% and of period-mean interface length <= 5% over held-out rpm, AND better than the baseline "nearest-rpm field" by >= 2x in median.
  Power: the baseline (nearest rpm, 2.5 rpm away) is a real competitor: if the flow changes smoothly, interpolation beats it; if the flow switches regime between rpm (consistent with T5 roughness of dtmix), both fail and the test says so. EDR is smooth in rpm (T5: max jump 0.21 in log at 0.625 rpm), so a 10% bar is above the per-rpm noise but below the 2.5-rpm change (EDR changes ~2x over 15-37.5 rpm). Weakness: 8 interior held-out points; reported per rpm.
  Also reported: IRMAE-WD intrinsic dimension per rpm (1 = limit cycle), and whether snapshots one period apart coincide (periodicity check, relative L2 difference).
- REVISION (same day, before any data): user: the POD/DMD/dimensionality-reduction line was third-party advice, not a requirement; use my own context. My own question chain does not need a flow surrogate:
  Q1 (test L): is the level dependence of dtmix an artefact of measuring chi at the grid scale? Experiments measure macro-mixing at a probe scale, so chi at a FIXED physical scale is the physically meaningful definition; if it converges across levels, dtmix becomes a well-posed L10 target.
  Q2 (test P): is the rpm roughness sensitive dependence (dtmix a random variable) or smooth-steep structure? If random: the design objective is its expectation and mfbml needs a noise model; cheap L8 ensembles + few L10 runs is exactly the multi-fidelity setting.
  Q3: kLa ladder 47/29/22 (ratio of successive changes 0.39) may converge geometrically (archive: "kLa convergence plausible"); test with the Sc ~ 500 caveat once P says how noisy kLa is.
  => Test M as a dimensionality-reduction test is WITHDRAWN (never run). The 10 flow_l8 runs (jobs 7277961-70, ~250 core-h) are kept for a periodicity check that feeds Q2:
  PRE-REGISTERED test Q (periodicity per rpm): relative L2 difference of (u, v) between phase-matched snapshots one period apart (linear time interpolation, 100 frames/period), median over phases and the 4 recorded periods. PERIODIC if < 1%, APERIODIC if > 10%, else quasi-periodic. Power: phase-interpolation error with 100 frames/period is ~(2 pi/100)^2 ~ 0.4% for a smooth signal, below the 1% bar (checked on the data: the error of interpolating a snapshot from its neighbours is reported as the floor; if the floor exceeds 1%, the PERIODIC bar is raised to 2x the floor and that is disclosed).
  Reading with P: periodic flow + smooth P -> the roughness is steep deterministic structure in rpm; periodic flow + noisy P -> impossible for a stable limit cycle unless multistable, check; aperiodic flow + noisy P -> flow chaos makes dtmix a random variable.
- The particle-estimator agent (already running) is kept only to finish its benchmarks; its use (test S) is decided after P and L.

## 2026-10-10 (cont.) — Is dtmix a well-posed target up to L10? Self-check, then test K (pre-registered)
- Self-check on existing data (dtmix_0.95, theta 7, same protocol): 32.5 rpm L6..L10 = 14.3, 28.5, 68.3, 128.3, 190.7 s (increments +14, +40, +60, +62). 25 rpm L6..L9 = 11.9, 29.8, 113.3, 303.6. 37.5 rpm L6..L9 = 17.5, 33.9, 29.0, 51.3. Reading: at 32.5 rpm the L8 -> L10 increments are nearly constant (+60, +62), the signature of t_mix ~ (1/lambda) ln(L/dx) for chaotic stretching (each halving of dx adds ln2/lambda); but L6 -> L8 and the 25 rpm ladder grow faster, and 37.5 rpm is non-monotone. So the log law is suggestive at one rpm, NOT established. Not presented as evidence.
- Why the grid, not physics, may set dtmix (first principles): stirring stretches tracer blobs into filaments that thin exponentially; variance is destroyed only when filaments reach the scale where diffusion acts in one stretching time, the Batchelor scale l_B = sqrt(D/gamma) [numbers here WRONG, guessed gamma; corrected below: 4-10 um from our eps]. Cells are ~0.25-1 mm at L10-L8. The cell Peclet number U dx/D ~ 1e4-1e5 (resolution needs ~2). So filaments reach the CELL size long before l_B, and the advection scheme's numerical diffusion (D_num, set by dx) destroys the variance. Prediction: dtmix depends on dx and not on D at every affordable level.
- PRE-REGISTERED test K (diffusivity scan), compute ~2.5k core-h, needs a code change (runtime multiplier tracer_D_scale on D_tracer_1, and oxy_D_scale on D_oxy_1; default 1, bit-identical results at 1 to be checked) and a lean rebuild:
  32.5 rpm, theta 7, kmix protocol; D scale s in {1, 1e1, 1e2, 1e3, 1e4, 1e5} at L6, L7, L8 (18 runs) and s in {1, 1e3, 1e5} at L9 (3 runs). Oxygen scaled by the same s in the same runs (kLa read too).
  Predictions: (i) GRID-CONTROLLED regime: for s below a level-dependent s*(L), dtmix is flat in s (changes < 10% from s = 1) and differs between levels (as now); (ii) PHYSICS-CONTROLLED regime: for s above s*(L), dtmix decreases with s and the levels AGREE (within 10% between the two finest levels at that s); (iii) s*(L) rises with refinement (D_num falls with dx).
  Verdict "dtmix is NOT a well-posed target at physical D for L <= 10" if (i) holds at s = 1 and 10 for L7-L9 AND (ii) holds at s = 1e5 for L8-L9. If dtmix at L8 changes > 30% between s = 1 and 10, the premise D_num >> D is false and the verdict is rejected.
  Power: the test has a built-in positive control: at large s the physical diffusion must take over and collapse the levels; if it does not, the test cannot distinguish "numerical diffusion" from "flow not converged" and says so. Flat-in-s at s = 1-10 is a strong, falsifiable statement: a physically controlled mixing time would change by ~ln(10)/ln(Pe) ~ 13% per decade under the log law, ~3x the 0.02 replicate log-sd at 32.5 rpm (phase-locked, CV < 1%), so a 10% bar can detect it (stated threshold, not tuned).
  Extrapolation to L10: the same rule applied with s*(L): if s*(L) grows ~2-4x per level, L10 is still far from s = 1; the minimum level where s* = 1 is estimated and reported (expected ~L15-16 from eta_B / dx, see correction below).
- CORRECTION (user challenge, same day): my Batchelor-scale numbers (20-50 um, "needs ~L13-14") were from a GUESSED strain rate (0.2-1 1/s) and are wrong. Standard definition (Batchelor 1959, as Wikipedia): eta_B = (nu D^2/eps)^(1/4) = eta_K/sqrt(Sc). It equals sqrt(D/gamma) when gamma = sqrt(eps/nu) (Kolmogorov strain rate), so the two forms agree once gamma is taken from the data, not guessed. From our own domain-mean dissipation (scripts/diag_batchelor.py -> experiments/multifidelity/diag_batchelor.json; eps = edr_mean_t/rho; gamma = sqrt(eps/nu) = 4-24 1/s):
  - eta_K = 200-500 um; dx/eta_K at L10 = 0.6-1.2: the velocity field is resolved near the Kolmogorov scale at L10 (consistent with tau/EDR changing little L9 -> L10).
  - Tracer eta_B = 4-10 um (Sc ~ 2300); oxygen eta_B = 9-22 um (Sc ~ 530). dx/eta_B(tracer) at L10 = 27-58, at L8 = 104-165. Resolving eta_B (dx ~ eta_B) would need dx ~ 4-10 um, i.e. L ~ 15-16 (not 13-14).
  - Caveats: eps is the DOMAIN mean (air included) -> liquid eps larger, scales smaller by up to ~1.4x; eps itself grows with level (32.5 rpm: 6.6e-5 at L8 -> 1.9e-4 at L10), but eta ~ eps^(-1/4) moves only ~1.3x. Local eps varies; this is an order-of-magnitude statement.
  - This is evidence that the scalar is under-resolved by 1.5 orders of magnitude at L10. It is NOT yet evidence that numerical diffusion sets dtmix; tests K and G address that.
- PRE-REGISTERED test G (scalar gradient scale, zero extra compute; data from scal_l8/scal_l9 and rec_l8 snapshots): replaces my hand-waved "sheets of pure dye" picture with a measurement. Scalar microscale lambda_c(t) = sqrt(<c'^2> / <|grad c|^2>) over liquid cells (central differences on the uniform snapshot grid), after release.
  GRID-CONTROLLED if, during the decay (chi 0.5 -> 0.95), lambda_c/dx is O(1-3) at both L8 and L9 AND lambda_c(L8)/lambda_c(L9) is within 1.5-2.5 (i.e. lambda_c scales with dx). NOT grid-controlled if lambda_c is the same physical length at both levels (ratio 0.8-1.25) and >= 5 dx.
  Power: the two outcomes differ by a factor 2 in the L8/L9 ratio, measured from thousands of cells per snapshot; snapshot-to-snapshot scatter is reported and must be below the 0.3 margin.
  Also reported: the fraction of liquid cells with 0.1 < c < 0.9 vs time (direct evidence for or against sharp striations).
- Quantisation check (chi_scale.py quant, Sonnet; spot-checked by me): solver dtmix = first logged row past the threshold, t_out = 0.02 nd (0.053-0.13 s). Rebuilt dtmix equals results.json exactly; linear interpolation changes dtmix by <= 0.09 s, i.e. <= 0.47% (chi 0.50), 0.22% (0.75), 0.18% (0.95) over the 29 kmix_l8/t5_l8 runs. Time quantisation cannot explain the T5 roughness (0.18-0.35 in log). One alternative explanation ruled out.
