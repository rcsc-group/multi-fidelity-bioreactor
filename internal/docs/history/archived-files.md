# Archived files

Files removed in the 2026-11 clean-up because no code, test or config used them. They are kept in the
git tag `archive/pre-lean-2026-11`. To get one back:

```bash
git show archive/pre-lean-2026-11:config/chain_mpi_f6.yaml > config/chain_mpi_f6.yaml
```

## config/

| File | What it was |
|---|---|
| `chain_f4_bug2_diag.yaml` | Two-segment chain that reproduced a blow-up at the mixing-time restart. |
| `chain_mpi_f6.yaml` | Fidelity-6 MPI chain, two segments with checkpoint restart. |
| `chain_mpi_fresh_diag.yaml` | Fresh MPI run past the old crash time t = 13.3. |
| `chain_mpi_smoke_diag.yaml` | The smoke-test parameters, run fresh past t = 13.3. |
| `chain_mpi_smoke_f3.yaml`, `chain_mpi_smoke_f3_4ranks.yaml` | MPI chain smoke test at fidelity 3 (16 and 4 ranks). |
| `chain_mpi_smoke_f4.yaml`, `chain_mpi_smoke_f4_n5.yaml` | MPI chain smoke test at fidelity 4 (and with 5 transition cycles). |
| `chain_restart_recovery_baseline_L8.yaml`, `chain_restart_recovery_test_L8.yaml` | Level-8 version of the restart-recovery pair (the level-7 pair is still in `config/`). |
| `sweep_fb_fill_l7.json`, `sweep_fb_fill_l7_v3_f03.json` | Early fill-fraction sweeps at level 7 (superseded by `_v2` and the `_mpi_ckpt` sweeps). |
| `tau_sweep_sentinel.sh` | Watcher script for one old shear-stress sweep. |

## scripts/

One-off launchers and checks of finished campaigns (submit, resubmit, smoke test, ETA, extension of a run). Rule: no
code, config, test or talk tool names them, and no file kept in the repo carries their name (no kept result
came from them). The diary may cite them; get them back from the tag as above (path `scripts/<name>.py`).

| Script | Last change | What it did (its own first line) |
|---|---|---|
| `_l9_l10_autopilot.py` | 2026-09-04 | Autonomous coordinator for the L9 warm-start rerun + L10 pause policy. |
| `_postprocess_l10_stalled.py` | 2026-07-21 | One-off: postprocess the 3 stalled L10 seg0 runs from already-complete raw data. |
| `_submit_restart_recovery_chains.py` | 2026-08-20 | One-off submission helper for the checkpoint-restart recovery test |
| `eta_chain_now.py` | 2026-09-09 | seg2 (25.0rpm, run bce29aa5, job 6114128): started 09:59:06, now 3:28:45 elapsed |
| `eta_check.py` | 2026-09-08 |  |
| `eta_full.py` | 2026-09-08 | Best case: cold_30 and cold_37.5 slot into the two 32-CPU vacancies as |
| `eta_reordered.py` | 2026-09-08 |  |
| `extend_censored_settling_runs_2.py` | 2026-09-12 | Second-round extension: of the 7 conditions extended by |
| `extend_fig8_late_time.py` | 2026-09-13 | Fig 8 replica investigation: our replica (t/Tp~34-37) is compared at a |
| `finalize_grid_results.py` | 2026-06-18 | Update hypothesis_ledger.json and run_provenance.csv from completed grid runs. |
| `mpi_diag_submit.py` | 2026-06-07 | Diagnostic: test whether MPI binary hangs at t_mix with 16 ranks vs 4 ranks. |
| `mpi_fix_verify.py` | 2026-06-07 | Submit verification run for the MPI movies_output fix. |
| `mpi_same_omega_restart.py` | 2026-08-03 | Falsifiable experiment: same-omega_b MPI restart (su=1, no velocity rescaling). |
| `mpi_smoke_test.py` | 2026-08-03 | 2-segment L5 MPI smoke test. |
| `resubmit_angle_chain_final.py` | 2026-09-12 | Resubmit the angle-sweep chain's final segment (theta=2deg), killed by |
| `resubmit_dtheta_decreasing_cancelled.py` | 2026-09-07 | Resubmit the 6 decreasing-direction dtheta transitions cancelled by |
| `resubmit_dtheta_decreasing_cancelled2.py` | 2026-09-07 | Second resubmission of the 2 remaining cancelled decreasing-direction |
| `resubmit_dtheta_th4.py` | 2026-09-07 | Third resubmission attempt for theta_source=4 -> target=2 (diary.md |
| `resubmit_dtheta_th4_th69.py` | 2026-09-07 | Resubmit theta_source=6.9 -> target=4 (cancelled by QOSMaxCpuPerUserLimit, |
| `resubmit_dtheta_th4_th69_exclnode.py` | 2026-09-07 | Resubmit theta_source=6.9 -> target=4 (run_id=35e2e274), excluding |
| `resubmit_dtheta_th4_th7.py` | 2026-09-07 | Resubmit theta_source=7 -> target=4 (cancelled by CPU-cap contention, |
| `resubmit_dtheta_th4_th7_exclnode.py` | 2026-09-07 | Resubmit theta_source=7 -> target=4 (run_id=c60abc7c), excluding |
| `resubmit_dtheta_th5.py` | 2026-09-07 | Resubmission for theta_source=5 -> target=2 (diary.md 2026-09-07), |
| `resubmit_dtheta_th5_exclude_node.py` | 2026-09-07 | Resubmit theta_source=5 -> target=2 (run_id=0f4b3a08), excluding |
| `resubmit_p_fix_verification.py` | 2026-09-07 | Resubmit the full-length p-fix verification (run_id=dea0e98d), |
| `resubmit_ramp_fix_test.py` | 2026-09-07 | Resubmit the ramp-fix test (run_id=c4d05e5c), excluding node2336 |
| `resubmit_rpm_chain_seg7.py` | 2026-09-10 | Resubmit the L9 rpm chain's final segment (37.5rpm), killed by TIMEOUT |
| `resubmit_testc_trivial.py` | 2026-09-07 | Resubmit Test C's trivial restart (run_id=46493d1a), excluding node2336. |
| `resubmit_testc_trivial2.py` | 2026-09-07 | Resubmit Test C's trivial restart, verification run (run_id=54c658b0), |
| `smoke_cross_level_warmstart.py` | 2026-09-10 | Tiny smoke test for cross-level (L6->L7) warm-starting before any real |
| `smoke_phase_fix.py` | 2026-09-07 | Cheap L6 smoke test of the restart phase-offset fix before spending the |
| `smoke_test_l9.py` | 2026-09-05 | L9 smoke test: measure real wall-clock cost before committing to the |
| `smoke_video_test.py` | 2026-09-08 | Smoke test the video-enabled binary (correctness fixes + VIDEOS=1) at L6 |
| `submit_chain_video_demo.py` | 2026-09-07 | Visual chaining-fidelity check (diary.md 2026-09-07): a single, 2-segment |
| `submit_compat_grid.py` | 2026-06-18 | Submit MPI checkpoint-restart compatibility grid (H12–H17). |
| `submit_cross_level_L7toL8_fsfix.py` | 2026-09-14 | L7->L8 cross-level pilot, re-run with the fs-before-refine ordering fix |
| `submit_cross_level_pilot.py` | 2026-09-11 | Real cross-level warm-start pilot: L6 (converged) -> refine to L7 -> |
| `submit_cross_level_videos.py` | 2026-09-11 | Rerun with video enabled, to visually inspect the L6->L7 cross-level |
| `submit_dtheta_L7_check.py` | 2026-09-07 | Resolution-independence check for the Delta_theta_max settling-time |
| `submit_dtheta_L7_check_v2_fixed.py` | 2026-09-07 | Re-run of the L7 Delta_theta_max check (diary.md 2026-09-07), now using |
| `submit_dtheta_th4_target.py` | 2026-09-07 | Delta-theta_max monotonicity study, decreasing-direction mirror v2, |
| `submit_extended_grid.py` | 2026-06-18 | Submit extended MPI checkpoint-restart grid (H18–H25). |
| `submit_fig13a_maskfix_reruns.py` | 2026-09-15 | Rerun the L6 and L8 Fig13a sweeps on the mask-fixed binary (diary.md |
| `submit_fig13a_rampmatch.py` | 2026-09-02 | Redo of the fig13a_redo 9-point sweep with upstream's ramp matched |
| `submit_fig13a_redo.py` | 2026-09-01 | Redo of the Fig 13a(a) replica (tau/EDR vs RPM at theta=7deg, L8), |
| `submit_fig8_hist.py` | 2026-09-21 | Re-record Kim Fig. 8(b)/(c) at L10 with enough distinct phases to find the peak. |
| `submit_figA18.py` | 2026-09-16 | The 2D/3D pair for Kim's Fig. A.18. |
| `submit_figset.py` | 2026-09-21 | One run for Kim's Figs. 2-7: the reference condition, fully instrumented. |
| `submit_kick_noise_and_longrun.py` | 2026-09-07 | Two runs that together decide whether the escape is restart-specific and |
| `submit_kick_test.py` | 2026-09-07 | Velocity-kick test (diary.md 2026-09-07): single, continuous, |
| `submit_l10_48core_probe.py` | 2026-09-10 | 48-core (single node) L10 throughput calibration before committing to a |
| `submit_l10_cold.py` | 2026-10-03 | A clean L10 point at 32.5 rpm: cold start, Kim's protocol, release at cycle 80. |
| `submit_l10_fig13a_32rank.py` | 2026-09-13 | L10 Fig13a resubmit at 32 ranks instead of 48. The original 48-rank |
| `submit_l10_fig13a_sweep.py` | 2026-09-10 | L10 Fig13a replication: 9 independent cold starts, 32.5-rpm-condition |
| `submit_l10_scaling_probe.py` | 2026-09-10 | L10 scaling probe + real progress, combined. |
| `submit_l8_coldstart_vid.py` | 2026-09-11 | Matched L8 cold-start video reference (settling_baseline_L8 predates |
| `submit_l9_angle_chain.py` | 2026-09-10 | L9 angle-sweep chain: the SECOND axis of Kim et al.'s data (this |
| `submit_l9_chain.py` | 2026-09-08 | The actual Fig 13a chained-sweep readiness test: a real ascending chain, |
| `submit_l9_rpm17.5_full.py` | 2026-09-05 | L9 sweep, condition 1 of 9 (17.5rpm, first in the low-to-high warm-start |
| `submit_l9_rpm20_checkpoint.py` | 2026-09-05 | L9 sweep, condition 2 of 9 (20rpm), CHECKPOINTED (warm-started) from |
| `submit_l9_warmcold_compare.py` | 2026-09-08 | Fig 13a readiness check at L9: 3 targets, warm (single-hop restart from |
| `submit_mf_l10_angle.py` | 2026-09-30 | L10 at the three Fig 13(b) training angles, cross-level from the L9 states. |
| `submit_mf_l7_rpm.py` | 2026-09-30 | L7 rpm sweep at 7 deg: the LF side of the Fig 13(a) multi-fidelity demo. |
| `submit_nc10.py` | 2026-05-31 | Submit restart_nc10: same as fix29_verify2 but n_mix_cycles=10. |
| `submit_replicates.py` | 2026-10-03 | Replicate and spin-up-length runs for the noise model of the MF problem definition. |
| `submit_restart_bias_multi.py` | 2026-09-03 | Does the same-condition restart-chain bias result (diary.md 2026-09-03 |
| `submit_settling_study_falsification.py` | 2026-09-07 | Falsification test for the cold-start-beats-warm-start anomaly |
| `submit_settling_study_transitions_v2.py` | 2026-09-06 | Settling-time study, corrected transitions (2026-09-06): nearest- |
| `submit_t10.py` | 2026-10-17 | Test T10 (pre-registered, diary 2026-10-17): L10 kLa_1T_10 (and dtmix_0.50 where reached) at the x8 design poi |
| `submit_test_k.py` | 2026-10-10 | Test K (pre-registered, diary 2026-10-10): diffusivity scan. Is dtmix (and kLa) set by the grid or by physics? |
| `submit_trivial_restart_regression.py` | 2026-09-14 | Regression guard for the 2026-09-14 cross-level fix: confirm an ORDINARY |
