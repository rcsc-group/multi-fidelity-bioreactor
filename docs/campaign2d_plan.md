# Campaign 2D: plan for approval (DRAFT, nothing submitted)

Date: 2026-11-04. Branch: campaign-2d (not pushed, not merged). Status: waiting for the lead's approval.
Draft driver: scripts/campaign2d/campaign.py (plan, submit-stage --dry-run, status, ledger run; fit and next are
stubs: fit refuses outside a compute job; submit-stage L8 --go is not executed).

## 0. Scope and the lead's answers (2026-11-04)

- Design space: rocking angle theta 2-7 deg, rocking frequency 15-37.5 rpm (Kim's range).
  Theta 9 runs exist but are outside the domain: excluded from the fit, reported only.
- QoIs (each predicted over the full 2D map at L10):
  - tau_95: wall shear 95th percentile, from the hydrodynamic window (scripts/hydro_dataset_v2.py).
  - edr_95: dissipation 95th percentile. NOT SAVED PER RUN (checked 2026-11-04). Only the
    time-series means exist (edr_mean_t, edr_max_t, ediss_mean_qss). A second histogram pass is needed
    in src/BioReactor.c (the tau_95 two-pass 200-bin histogram, about line 3190, is the template).
    Until that change exists, the fit uses edr_mean_t and labels it as a stand-in, not edr_95.
  - kLa_1T_10, kLa_1T_25, kLa_1T_50: whole-period estimator (results.json keys). kLa_c (5-sample)
    carries 10-36 % phase noise (diary 2026-10-08), so it is reported, not fitted.
  - dtmix_0.50: first crossing of chi = 0.50. Definition conflict with Kim's dtmix_strict (3-5x,
    unresolved, diary 2026-10-22). See open question Q4.
- No budget cap. Target: about one week. Replay used only where the diary supports it (section 3).

## 1. Design

### 1.1 Criterion and map
- Work on ln(QoI). The fit (FullBayesKRRLRGP, reused) gives a Student-t mixture per point with
  mean mu_g and scale sc_g per node g and weight w_g. The latent ln-variance is
  v = sum_g w_g sc_g^2 nu/(nu-2) + Var_g(mu_g).
- 2 sigma/mean on the linear scale is 2 sqrt(exp(v) - 1), about 2 sqrt(v) for small v.
  So the criterion 2 sigma/mean <= 0.30 is the same as sd_ln <= 0.15 at every point.
- Criterion on the continuous domain (amended 2026-11-04, lead): the max of 2 sigma/mean over the whole square
  theta in [2, 7] x rpm in [15, 37.5], NOT a grid max. Evaluate on a dense 101 x 101 grid, refine the top 10
  points with scipy.optimize.minimize (L-BFGS-B, bounds = domain) on -2 sigma/mean. Report the max and its
  location per QoI. Implementation: campaign.py max_crit().
- sigma here is the surrogate (latent) uncertainty, not the run-to-run scatter. See Q1.

### 1.2 L8 (LHS)
- N8 = 24 new maximin Latin-hypercube points (seed 20261104). Existing L8 sites inside the domain: 126
  (campaign.py status). The existing set is a 10-rpm x 4-theta tensor grid plus a 0.625-rpm grid at
  theta 7, so the LHS fills the theta interior (2.0-2.9, 3.4-6.9 deg).
- Why 24: the run-to-run sd at L8 is 0.095 in ln for kLa_1T_10 at theta 7 (diary 2026-10-17, test E8).
  Stop rule for L8 alone: leave-one-out RMSE in ln <= 0.10. If not met, add 8 points per batch, max 2 batches.
- Cost: 24 x 55 = 1,320 core-h (cold, 16 ranks).

### 1.3 L9 (adaptive)
- Candidate sites = L8 sites only. A warm start needs a settled L8 checkpoint at the same (rpm, theta);
  rpm-changing restarts bias dtmix (diary 2026-10-10: seam bias +28-57 %). Candidates elsewhere need a cold run.
- Score(x) = max_q [2 sigma_q(x)/mean_q(x)] / cost_L9(x), with cost_L9 = 220 (warm) or 580 (cold).
  Scalar source (amended 2026-11-04, section 3.1): kLa_1T_10, kLa_1T_25 and dtmix_0.50 may come from replay at
  rpm < 35; kLa_1T_50 and rpm >= 35 need the coupled tail. Batch size 4, 3 batches, 12 points.

### 1.4 L10 (seeded, then adaptive)
- Seed batch (4 points, fixed before L9 starts, existing L8 sites):
  (22.5, 4), (22.5, 5.5), (15, 5.5), (37.5, 5.5). Reason: the 22.5 rpm ridge (dissipation about 2.5x the
  smooth trend at L9 and L10; L10 pair EDR difference -0.26 ln, diary 2026-10-23) and the domain
  corners, which no L10 run covers.
- Adaptive batches of 2, max 3 batches (8 full points). Score as L9, cost_L10 = 1,500 (full scalar point)
  or 500 (hydro-only point: tau_95 and edr only).
- Tail: every full L10 point needs 12 warm-up cycles plus about 40 post-release cycles (dtmix_0.50 at 17.5 rpm
  needs about 38 cycles, diary 2026-10-21).

### 1.5 Stop and gates
- Gates (amended 2026-11-04, lead). Gated QoIs: tau_95, edr_95, kLa_1T_10. Each must have max 2 sigma/mean <= 0.30
  over the whole continuous domain (section 1.1). kLa_1T_25, kLa_1T_50 and dtmix_0.50 are reported only (map and
  max), not gated.
- Stop (a): the three gated QoIs have max 2 sigma/mean <= 0.30 on the L10 map (continuous max, section 1.1).
- Stop (b): an L10 batch lowers the max by less than 10 % and the max is still above 0.30. Report, do not add.
- No budget cap. Each batch needs explicit lead approval (--go). Nothing is submitted without it.

### 1.6 Expectation (to be tested by the first fit, not claimed)
- tau_95, edr: likely reachable (H3 4-7 %, H4 test set 2-7 %). kLa_1T_10: possible (W9 warm start, 3 of 3 rpm).
- kLa_1T_25, kLa_1T_50, dtmix_0.50 at L10: probably NOT reachable in one week (cold L10 point 2.5-4k core-h; W9 scatter 0.33 ln).

## 2. Parallelism, limits and wall-clock

- Partition batch. QOS priority (sacctmgr: MaxTRESPU cpu = 312). batch2 (cpu = 512) only with the lead's
  approval. mbessa-condo: forbidden (never planned, never submitted). node2356 excluded in every sbatch.
- One node per job (config/slurm_mpi_template.sh pins nodes = 1: multi-node output corruption, 2026-09-15).
- Ranks: 16 (L8), 32 (L9, L10). Concurrency cap: 9 x 32 = 288 cores (under the 312 cap).
- Walltime: Exploratory cap 48 h (memory). L8 16 h. L9 warm 40 h (run 6-7 h); L9 cold 40 h (run 18 h).
  L10: 47 h per segment. A 17.5 rpm L10 point needs about 100 h wall (62 core-h per cycle / 32 ranks x 52 cycles),
  so it runs as 3 chained segments (submit_continue.py, validated 2026-10-21).
- Queue history (sacct -S 2026-10-01, 657 bioreactor-mpi jobs with a start time):
  submit-to-start median 0.20 h, p75 1.9 h, p90 5.6 h, max 47 h.
  446 jobs started within 1 h, 152 within 6 h, 48 within 24 h, 11 after 24 h.
  Family-wide peak 832 cores in 30-min bins (includes batch2 runs); mean 157 cores while active.
  The 2026-10-21 event (18 jobs pending about 14 h, fairshare) is the worst case seen.
- Full-cluster time: 22,590 core-h / 288 cores = 78 h.
- Critical path:
  - Day 0: start the L8 LHS (24 x 16 ranks; 18 fit in the 288 cap), and the L10 seed (4 x 32 ranks, 128 cores).
    Total 384 + 128 cores, so the L10 seed waits for L8 jobs. Expect about 1 day for L8.
  - Day 1-2: L9 batch 1 (4 warm points, about 7 h each).
  - Day 2-6: L10 seed (chained segments, 3-4 days including queue).
  - Day 6-9: L10 batch 2 and L9 batch 3 (queue 0.2-5.6 h each).
  - Estimate: 7 days if the seed runs in parallel and the L9 and L10 batches overlap. 10 days if all
    adaptive batches run in sequence. The one-week target holds only with parallel seeding and at most
    two adaptive L10 batches.
- Deadline: the DFD talk is about 2026-11-22 (18 days from today). A one-week campaign approved today
  finishes about 2026-11-11 to 2026-11-14.

## 3. Speed levers and their validity

Decision key: USE = in the campaign for the QoIs listed. NO = not used. Cost saving is stated per lever.

| Lever | Evidence (diary date, numbers) | Failures | Decision |
|---|---|---|---|
| Replay, L8 (one recorded period looped, scalar only) | 2026-10-21 night (replay_lf_mf): RMS ln error at n = 4/8/16 runs: kLa_1T_10 0.08/0.08/0.07; dtmix_0.50 0.11/0.10/0.10. 2026-10-31 (chain at L8, 150 cycles): kLa_1T_10 within 10 % at 10/10 (worst +9.7 %); dtmix_0.50 9/10 (worst -10.4 %); dtmix_0.75 9/10. | 2026-10-21 (periodic, K = 1/2/4): kLa_1T_25 +51.3 / +35.1 / +15.4 % (FAIL). 2026-10-21 (L6, F2): max abs dC* 0.0375 (bar 0.02), dtmix +5.5 % (FAIL). 2026-10-28 (Kim-style rebuilt flux): kLa_1T_10 x1.23, kLa_1T_25 x2.26 (not used). The replay overestimates kLa_1T_25 by about 11 % (mean ln(HF/replay) = -0.108, 2026-10-21 night). | USE (amended 2026-11-04, lead) for kLa_1T_10, kLa_1T_25 and dtmix_0.50 at L8, rpm < 35 only, after the flow check of 3.1. NO for rpm >= 35: coupled tail. kLa_1T_50 is reported only (lead, 2026-11-04): taken where a coupled tail exists anyway, never forces one. Saving: zero under the decided rule, because tau_95 and edr_95 stay coupled-only (see 3.1). |
| Replay, L9 | rc9 records and replays done (jobs 8249920-8249940). Diary 2026-11-04 (check_replay_kla25_L9.json): kLa_1T_10 10/10 within 10 %, kLa_1T_25 9/10, dtmix_0.50 9/10 (37.5 rpm -13.2 %), kLa_1T_50 7/10 (2 NaN). Flow check: tau_95 1.0-1.7 % low at every rpm (pre-registered limit 1 %: FAIL), posY_max within 0.8 %. | Pre-registered 1 % flow limit failed (see 3.1). | USE (amended 2026-11-04, lead) for kLa_1T_10, kLa_1T_25 and dtmix_0.50 at rpm < 35, after the flow check of 3.1. NO for kLa_1T_50 and rpm >= 35. Saving against cold L9 (580) only if DNS-81 plus replay (rc9 DNS 220-370 core-h plus replay tail) is cheaper; see 3.1. |
| Replay, L10 | Job 7893240 COMPLETED (log dated 2026-10-25): kLa_1T_10 = 21.84 vs coupled 22.05 (cold) and 22.43 (warm): -0.9 % and -2.6 %. Diary never recorded this result. | dtmix NaN with kLa present (unexplained). Cost 927 core-h for 22 cycles (1 d 5 h x 32 ranks), about the same per cycle as the coupled L10 run (17-64 core-h/cycle). | NO at L10: no saving measured, and the mixing time is missing. Open: explain the NaN (Q6). |
| Checkpointing across conditions (warm start from the previous dump, W9-type) | 2026-10-24 (W9, L8 -> L9, 12 cycles, release after 12): kLa_1T_10 +4.3 / +6.2 / -8.0 % at 17.5 / 25 / 32.5 rpm (PASS, 3/3 within 15 %). 2026-10-21 (settled check, 13 L10 points): median abs difference tau_95 1.7 %, tau_mean 1.1 %, EDR 0.2 %; 12 of 13 within 6 %. 2026-10-17 (R12): release at cycle 12 vs 80: kLa_1T_10 +0.006 (se 0.032). | 2026-10-24 (W9): dtmix_0.50 -21.9 / +0.2 / -45.7 % (FAIL); kLa_1T_25 -4.0 / +8.2 / +113.7 % (FAIL). W9 replicate (2026-10-24): dtmix_0.50 ln sd 0.33 across 8/12/18 warm-up cycles (intrinsic). 2026-10-21: 22.5 rpm (6-cycle xlevel start) NOT settled, tau_95 +21 %. 2026-10-12 (CURRENT STATE): 6-cycle L8 -> L10 start FAIL (tau_mean +7.3 %). | USE for tau_95, edr and kLa_1T_10, with >= 12 warm-up cycles and the settled check. NO for dtmix_0.50, kLa_1T_25, kLa_1T_50 (cold 80-cycle protocol required, 2.5-4k core-h per L10 point). |
| Chaining (same-condition restart: restart_continue = 1, *_prev = current, tracers kept) | 2026-10-21 (L8, theta 7, 25 rpm; joined run vs uninterrupted): dtmix_0.50 0.00 % (19.1264 s); kLa_1T_10 30.2062 vs 30.2061; kLa_1T_25 14.4296 vs 14.4295; dtmix_0.95 -0.28 %; kLa_1T_50 -2.2 %. | One rpm only. 2026-10-18 note: unexplained ~8 % artifact for restart_continue = 1 (BioReactor.c ~line 1085), not closed in the diary. Missing *_prev gives seam bias dtmix +28-57 %, kLa -25 to -39 % (audit 2026-10-10). | USE to split long runs into segments (L10), with a seam gate: the overlap window must agree within 1 % on each QoI used. Not evidence for other rpm. |
| Warm start across levels (L8 -> L10, L9 -> L10) | Same numbers as the checkpointing row. 2026-10-12 (CURRENT STATE): 12-13 cycles, all QoIs within 5 % (tau_mean +3.6 % at 17.5 rpm after 13 cycles). | 6-cycle starts fail (see above). 22.5 rpm unsettled. | USE only for tau_95, edr and kLa_1T_10, with >= 12 cycles after the level change. NO for dtmix, and NO 6-cycle starts. |
| Early release (cycle 12 instead of 80) | 2026-10-17 (R12, 12 rpm, L8): kLa_1T_10 +0.006 (se 0.032); dtmix_0.50 -8 % (2.7 se). | dtmix_0.95 -43 %. | USE for kLa_1T_10 only. dtmix_0.50 keeps the 80-cycle protocol. |

### 3.1 Amendment (2026-11-04, lead decision): replay at L8 and L9

- Source rule. The looped replay is the source of kLa_1T_10, kLa_1T_25 and dtmix_0.50 at L8 and L9 for rpm < 35.
  For rpm >= 35, and for kLa_1T_50 at any rpm, the coupled DNS tail is used. Basis: diary 2026-11-04
  (experiments/multifidelity/check_replay_kla25_L8.json and _L9.json).
- Flow check amendment. Each replay run must pass a flow check before release: tau_95 and posY_max within 2 % of
  its recording DNS. The pre-registered limit was 1 %. L9 failed it (tau_95 1.0-1.7 % low at every rpm, a
  systematic offset, same sign as L8 at 0.4-0.9 %). The 2 % limit is a change made after seeing the L9 data; it is
  recorded here as an amendment, not as a pass of the original rule.
- Speed. Replay tail about 2.4x faster than the coupled tail (diary 2026-10-31).
- kLa_1T_50 (lead decision, 2026-11-04): taken only from runs that have a coupled tail anyway (rpm >= 35, L10).
  It is reported, not gated, and it does not force a coupled tail below 35 rpm. Q10 is answered by this.
- Budget consequence (campaign.py replay_saving). A point saves the replay price only if none of its QoIs
  needs the coupled tail. tau_95 and edr_95 are coupled-only in the decided rule, so every point still needs the
  coupled tail: the saving is 0 core-h and the budget is unchanged (section 4.2).
- Alternative (Q11, CLOSED: not taken, lead 2026-11-04): replay for tau_95 and edr_95 below 35 rpm would have saved
  440 core-h on the 24 new L8 points (22 below 35 rpm x 20 core-h). Not used.

Cost saving summary: checkpointing and early release save about 2x per L10 point (12 warm cycles instead of
80 cold). Replay saves nothing for the planned set (3.1). Chaining is a necessity (48 h cap), not a saving.

## 4. Cost ledger and budget

### 4.1 Measured per-point costs (sacct -X, Elapsed x AllocCPUS, 2026-11-04)
- L8 cold 230 cycles, 16 ranks: 7156006 = 52 (diary kmix ~57). L8 replay-chain DNS to cycle 81: 8184447 = 29, 8184458 = 11.
- L9 warm 52 cycles, 32 ranks: 7847552 = 191, 7847553 = 223. L9 cold 230 cycles: ~580 (diary 2026-10-19). L9 rc9 DNS 16 ranks: 8249920 = 370.
- L10 warm full, 30-40 cycles after release: 7512023 = 1,500; 7512031 = 1,162; 7720889 = 1,500; 7720883 = 892.
- L10 one-period record 7806441 = 60; L10 replay 22 cycles 7893240 = 927; per-cycle 17-64 by rpm (diary 2026-10-22).

### 4.2 Plan budget (central, from `plan`)
| Item | core-h |
|---|---|
| L8 new (24, cold) | 1,320 |
| L9 new (12, warm) | 2,640 |
| L9 cold surcharge (4 points for dtmix/kLa_25/50) | 1,440 |
| L10 full (8) | 12,000 |
| L10 hydro-only (6) | 3,000 |
| edr_95 backfill (20 L10 runs x 3 cycles, about 120 L8 runs x 0.75) | 2,190 |
| Total | 22,590 |

Range from the low and high cost bounds: about 14k (low) to about 42k (high). No cap (lead).
Replay saving in this total: 0 core-h (decided rule, section 3.1; kLa_1T_50 reported only, tau_95 and edr_95 coupled).
Q11 closed (lead): tau_95 and edr_95 stay coupled. No replay saving is taken.

### 4.3 Ledger and Kim counterfactual
- `campaign.py ledger` sums sacct over the job IDs in scripts/campaign2d/state/jobs.jsonl. That file is
  written only by a real (--go) submission, so it is empty now (ledger prints 0.0). It logs core-h per job and the running total.
- Kim counterfactual: KIM_CORE_H_PER_CASE = 120 (Kim 2025 p.6). PENDING AUDIT by another agent. Not final.
  22,590 / 120 = 188 Kim-style cases. One planned L10 full point (1,500) = 12.5 Kim cases. The two cost
  bases are different fidelities and not comparable; the number only shows scale.
- Account context: sacct family "bioreactor-mpi" since 2026-10-01: 81,452 core-h in 677 jobs. This includes
  other BioReactor work. The ledger uses job IDs, not names, so this number is not the campaign.

## 5. Risks and pre-flight checks (before any --go)

1. Binary. The kmix params point to /oscar/scratch/eaguerov/BioReactor-mpi-lean-f1c11e0, which no longer
   exists (purged). The persistent copy is /oscar/data/dharri15/eaguerov/bin/BioReactor-mpi-lean-f1c11e0
   (mtime 2026-09-21). The draft uses this copy. But src/BioReactor.c last changed 2026-10-29 (commit 61b6e87:
   Kim seed default, bulk-liquid sums). That changes dtmix by factors 0.886 / 0.931 / 0.950 (diary 2026-10-29).
   Decide the binary before submission (Q5). Rule (memory): binary mtime newer than the source commit.
2. edr_95 needs a source change and a rebuild (Basilisk, qcc from /oscar/data/dharri15/eaguerov/basilisk).
   The new histogram must not change the flow: check with cmp of tr_oxy.dat and shear_stress.dat against
   a reference run. The 2,915-run database needs the 3-cycle backfill (about 2.2k core-h).
3. Restart sanity. Face fields are not in dump(): recompute cs and fs with solid() after every restore (BioReactor3D
   CLAUDE.md). Every restart .err must contain "fs sanity [post-solid]: OK" (seen in 7893240).
4. Smoke tests before long jobs: an L6 three-cycle run for each run type, and a 3-cycle L8 restart from a
   checkpoint, before any 230-cycle or L10 job (memory: verify pipeline first).
5. Scratch purge: /oscar/scratch is purged after about 30 days. Campaign outputs and checkpoints go to
   /oscar/data. An L10 record is 324 GB per period (diary 2026-10-24); 26 TB free.
6. 37.5 rpm: aperiodic flow at L8 (diary 2026-10-12, test Q), largest replay error (2026-10-31). It is a
   corner of the L10 seed and of the criterion map. Expect it to dominate the max.
7. Theta 2 (flagged in diary 2026-10-12: low-frequency variability, cycle sd 3.5-3.75 %). Its points
   may dominate the max. Q7.
8. 22.5 rpm ridge at the scalar level (kLa_1T_10 difference 0.29 at L9, diary 2026-10-18). Seeded at L10.
9. dtmix intrinsic floor: ln sd 0.33 at L9 (W9). A latent criterion cannot be met for dtmix_0.50 without
   replicates (Q1).
10. Pre-registration: the L10 acceptance rule above must be written into the diary before the first fit.
11. Kim 120 core-h is unaudited (section 4.3).
12. L10 replay dtmix NaN is unexplained (section 3).
13. Priority QOS 312-core cap and the 2026-10-21 queue event. The seed and L9 batch share the cap.

Pre-flight checklist (all must pass before any --go): plan approved; binary decided and mtime checked;
smoke runs passed; dry-run lists node2356 excluded and no mbessa-condo; ledger reads 0.0.

## 6. Open questions for the lead

- Q1. Criterion. Is 2 sigma/mean on the latent (surrogate) field the target, or on the predictive of a
  single observed run? For dtmix_0.50 at L9 the run scatter is 0.33 ln, so the observed-run criterion
  cannot be met. Accept the latent criterion and report the floor, or fund replicates?
- Q2. L10 scalars beyond kLa_1T_10. Warm starts are not validated for kLa_1T_25, kLa_1T_50 or dtmix_0.50.
  Cold L10 points cost 2.5-4k core-h each. Fund about 4 cold L10 points (about 12k core-h), or restrict the
  L10 criterion to tau_95, edr and kLa_1T_10, and check the others at L9?
- Q3. edr_95: approve the source change (second histogram, rebuild, flow cmp, 2.2k core-h backfill)? Or
  use edr_mean_t as the stand-in and rename the QoI?
- Q4. dtmix definition: our first crossing of chi = 0.50, or Kim's dtmix_strict? The offset is 3-5x and unexplained.
- Q5. Binary: keep the f1c11e0 copy (matches the 2,915-run database, predates the 10-29 seed fix), or build
  from the current source and rescale or rerun the database dtmix?
- Q6. L10 replay (7893240): the dtmix NaN is unexplained. Confirm the L10 replay is out of scope until explained.
- Q7. Theta 2: include in the domain (flagged), or start at theta 4?
- Q8. batch2 (512-core association) for a faster week? Not assumed in the estimate.
- Q10. ANSWERED 2026-11-04 (lead): kLa_1T_50 is reported only, from runs with a coupled tail anyway (section 3.1).
- Q11. CLOSED 2026-11-04 (lead): NO. tau_95 and edr_95 stay on the coupled tail. The 440 core-h alternative
  (2 % of budget) is not taken and is not a lever for the gated QoIs.
- Q9. Deadline. The DFD talk is about 2026-11-22. A one-week campaign approved today ends about 2026-11-11 to 11-14.
  Confirm this order of work is acceptable.

## 7. Provenance (stat, 2026-11-04)
- Modified on or after 2026-10-25 (ok): diary.md (11-02); scripts/replay_chain.py (10-31); scripts/replay_l10.py (10-24);
  scripts/test_t3_learning_curve.py (10-21 17:02); src/BioReactor.c (10-29, commit 61b6e87); scripts/postprocess.py (10-29);
  scripts/fig_mf_grid.py (10-29); logs/slurm_7893240.out (10-25); eec1f15 binary (10-11, ok for L10 warm runs).
- SUSPECT (before 2026-10-25): docs/mf_problem_definition.md (10-05); scripts/simulate.py, scripts/chain.py,
  scripts/mfbml_local/kernel_ridge_regression.py (10-01); scripts/mfbml_local/krr_lr_gpr_fullbayes.py (10-08, reused);
  scripts/test_h1_mfbml.py (10-13; BOUNDS theta 4-9, not reused); scripts/replay_lf.py (10-21); scripts/submit_l10_from_l8.py
  (10-21); scripts/hydro_dataset_v2.py (10-22); experiments/multifidelity/hydro_dataset_v2.json (10-08, not used);
  experiments/multifidelity/kla_mix_dataset.csv (10-16, not used); config/slurm_mpi_template.sh (09-21, nodes = 1 rule);
  bin/BioReactor-mpi-lean-f1c11e0 (09-21, older than the 10-29 source).

## 8. edr_95 backfill (list only, nothing submitted; lead 2026-11-04)
- Until the backfill exists, the fit keeps edr_mean_t as the stand-in for every row (labelled in fit output).
- List: scripts/campaign2d/state/edr_backfill_list.csv (script scripts/campaign2d/backfill_list.py). 27 distinct
  runs the fit uses for the EDR QoI: 10 L8, 10 L9, 7 L10. All 27 have a checkpoint.dump.
- Cost (sacct per-level median rate x 3 cycles x run period): L8 6, L9 125, L10 910, total about 1,040 core-h.
  The L10 rate is a median over mixed runs and is likely high; treat 910 as an upper bound.
- Restart: each run restarts from its last checkpoint for 3 cycles, binary BioReactor-mpi-lean-2ef5b07, reads edr_95.
