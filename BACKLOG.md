# Backlog

Ideas and loose ends not yet scheduled. Newest first within a section. When an
item starts, give it a diary entry; when it lands, delete it here.

## APS DFD 2026-11-25 (pivot 2026-10-08: vanilla mfbml, HF target)

- H3 variance propagation (owed, zero compute): model 2 uses only the mean of model 1, so the L10 band between design points is too narrow; add rho1^2 var9 (recursive co-kriging; read the reference before citing).
- kla_mix_dataset.csv needs a binary/protocol column (L9 theta-7 fig9 runs = prod-4b3a435, cold; e9 = lean-f1c11e0, T5 protocol; L10 T10 = lean-eec1f15 warm restart).
- Tell CCV about node2356 (killed 15 jobs 1 s after start, bursts 2026-10-16 18:48 and 2026-10-17 00:30; Slurm lists it MIXED, not drained).
- e9_l9_rpm22.6_th7 stopped at the 30 h limit (chi 0.937, t = 134): dtmix_0.95 censored. Continuation needs restart_continue=1 (unexplained ~8% artifact, BioReactor.c ~l.1085) or a longer first walltime; budget 36 h for L9 T5-protocol runs.
- Same-protocol L9 rpm grid (about 5.8k core-h) only if T10/E9 make L9 necessary; E9 says L9 roughness is local in rpm (22.5 yes, 30 no).
- L10 roughness near 22.5 rpm (L9 pair split 0.29 in kLa_1T_10 there): add L10 scalar points at 22-23 rpm after T10.
- F3b follow-ups: Ls = 11 replay (~1.4k core-h), a second rpm, an L10-flow record (cost: ~26 core-h per L10 cycle for the flow, records ~16x bigger than L8).
- The time-horizon fidelity (truncated chi curve) does not carry the late rate (test B); only
  revisit with a physical late-rate model, not a prior on r.

## Methods

- **Record/replay of the passive scalars (built 2026-10-16/17, branch replay-gate merged).**
  Exact non-looped replay works (F0, F1, MPI). Looping one period FAILS (F2: the L6 flow is not
  periodic, defect 9-30%), so it is not a cost saver. Open: F3b (scalar mesh Ls = 8, 9, 10 at a fixed
  L8 flow, record written, 57 GB). The Floquet-exponent idea needs a periodic flow and is parked.
- **Heteroscedastic noise in KRR-LR-GPR.** L10 cycle-to-cycle spread of
  tau_max grows ~35x from 17.5 to 37.5 rpm; the GPR noise is homoscedastic,
  so the band is too wide at low rpm and the pooled value is set by 37.5.
- **KRR-LR-GPR band collapses at small n_HF.** The likelihood matches Yi eq 23
  (checked 2026-10-08), but predict() scales the variance by the plug-in
  sigma2 = r'K^-1 r / n; with n_HF = 3 and a linear transfer it goes to ~0
  (test A coverage 0-57%). Fix candidates: prior on the residual amplitude,
  or n_HF >= 5.

## Figures / runs

- **tau_95 histogram quantisation.** tau_95 = upper edge of 200 linear bins up to the per-step tau_max; tau_max/tau_95 grows 3.6 (L8) -> 12-13 (L10), so one bin is 2-7% of tau_95, level-dependent (archive audit 2026-10-10). Fix before final L10 runs: exact per-step percentile (sort) or >= 2000 bins; re-derive tau_95 from snapshots where available.

- Bare `return;` in VIDEOS-only events movies_output / movies_output_tau (src/BioReactor.c ~2512, ~2623): the CLAUDE.md hazard (can stop the run). Lean production builds are unaffected (no VIDEOS). Fix before any VIDEOS run.

- **Fig 8(b) tau tails ~13x Kim's** after the transient cleared. Test for
  grid-scale noise with a volume metric (EDR field spectrum / |grad u|^2
  distribution) at L8/L9/L10 before touching anything near the wall.
- `submit_fig8_hist.py` still restarts from zero forcing: set every `*_prev`
  to the current condition (as submit_fig11_l10_kla.py does). Cost of not
  doing it: ~15 cycles of L10 (~15 h) per re-record.
- Fig 9 L10 anchor (32.5 rpm): segment 2 running (job 6763888); ~4
  segments total. Add the chain to plot_fig9/plot_fig11 when it lands.
- MF tau_rpm: replace the 6-cycle L10 HF points at 35/37.5 with the ~30-40
  L10 cycles logged by fig11_l10_rpm35/37.5 (pre-release window).

- L6->L7->L8 kLa ladder at one bad L9 point (27.5 rpm or theta=4): does it
  converge toward Kim like 32.5 rpm, or stay anomalous?
- `make replicas` target wrapping the README's regenerate list.
- Unsubmitted: `submit_figset.py` (Figs 2-7), `submit_fig14.py` (14/15),
  `submit_figA18.py`. Parked by request: A.16(c).
- Reproducibility of production kLa: at L4, OpenMP round-off moves kLa_50 by 2x (diag_omp_determinism.py). Test whether L7 kLa changes with MPI rank count (8 vs 16). If it does, kLa error bars need a round-off/chaos term.
- t_end period rounding: FIXED in source (2026-09-30). Production scratch binaries still use the old rule; picked up at the next production rebuild.
- Fig 13 MF figures still use the old buggy 3-grid u_bar (0.24/0.30). Faithful E&H gives 0.5-1.5 (diag_eh_ubar.py). Decide: apply it, or get L11 grids.
- kLa: whole-period estimator adopted in Figs 11/12 (2026-10-01). Still open: re-check mixing time (dtmix) the same way; the CI restart test still uses the 5-point kLa.
- src/BioReactor.c ~l.368 comments geometry_a as 'bag half-length (m)', but snapshots show the bag spans the full L0 = 1 = 0.25 m. Check and fix the comment (2026-10-11, test I agent).
- Kim offline comparison (2026-10-28): when jobs 8051486/8051487 end, compare kim_l8_rpm32.5_th7 with lf_l8_rpm32.5_th7 (tracer total, max concentration, kLa_1T_10, dtmix_0.50); then decide whether to emulate his T/100 snapshot hold and read the rest of TracerTest_test.c.

- 2026-10-29: Kim vs ours offline-vs-DNS test pre-registered (diary.md, 2026-10-29 PRE-REGISTERED); gates G0-G2 next, then T0/T1.
