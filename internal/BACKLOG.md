# Backlog

Open work after APS DFD (2026-11-25), ranked by ambition. When an item starts, give it a diary entry; when it
lands, delete it here. DFD-era entries were moved out on 2026-11-22 (see diary.md, same date).

## Ambitious tracks (user, 2026-11-22: "profiling campaign, or a GPU equivalent")

- **GPU solver.** The production solver runs on a uniform tree (refine(level < L), no adapt_wavelet), so a
  Basilisk GPU multigrid backend (basilisk/src/grid/gpu) is in reach. Blocker: the GPU backend has no
  embed.h support, and the bag wall is an embedded boundary. Candidate: in the bag frame the wall is fixed and
  rectangular, so it can be the domain boundary (non-square multigrid domain) with the rocking as a body force,
  which the solver already applies. Gate: L8 QoIs (tau_95, EDR, kLa, t_mix) vs the embed CPU solver within the
  replicate scatter, and a measured cost per cycle.
- **Profiling campaign.** Where do the core-h go per level (Poisson/multigrid, VOF, embed, tension, scalars,
  I/O, MPI)? Strong and weak scaling per level; rank count vs cost; find the cheapest layout per level.
  Output: a cost model per level that the BO cost curves can use.

## Urgent

- The 171 runs in runs/ that swept theta on pre-label binaries stay defective (scripts/bo2d/audit_cross_theta.py);
  none is in the kLa or hydro database. build/BioReactor-mpi* date from 2026-09 (pre-label): rebuild before use.

## Models (zero compute)

- Recursive co-kriging variance propagation (add rho1^2 var9; read the reference before citing).
- Heteroscedastic noise: L10 tau_max spread grows ~35x from 17.5 to 37.5 rpm.
- KRR-LR-GPR band collapse at n_HF = 3 (plug-in sigma2); prior on the residual amplitude, or n_HF >= 5.

## Physics and data

- 2D L10 truth over theta: Kim grid theta 3/5/6 (~2k core-h each), theta 2 warm from fig10_l9_th2, T10 tracer
  continuation of x8_l10_rpm32.5_th5.5. Cost basis: L10 theta 4, 97 min/cycle x 32 ranks. Cheaper after the
  cross-theta test passes.
- tau_95 quantisation: 200 linear bins, one bin = 2-7 % of tau_95, level-dependent. Exact percentile or >= 2000 bins.
- Fig 8(b) tau tails ~13x Kim's: test grid-scale noise with a volume metric at L8/L9/L10 first.
- F3b replay follow-ups: Ls = 11 (~1.4k core-h), a second rpm.
- kLa vs MPI rank count at L7 (8 vs 16): if it moves, error bars need a round-off/chaos term.
- L10 roughness near 22.5 rpm: L10 scalar points at 22-23 rpm.
- e9_l9_rpm22.6_th7 censored at 30 h: budget 36 h for L9 T5-protocol runs.

## Clean-up

- r1_cost_cv.py: L8 must exclude the 12 fill runs (reads 47, range(35) hardcoded).
- kla_mix_dataset.csv: add a binary/protocol column.
- Bare `return;` in VIDEOS-only events movies_output / movies_output_tau (src/BioReactor.c).
- 84 scripts with Path("/oscar/...") checkout paths -> Path(__file__).resolve().parents[k].
- `make replicas` target; geometry_a comment (half vs full length).
- Tell CCV about node2356.
