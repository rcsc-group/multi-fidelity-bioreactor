# Solver source (`src/`)

Basilisk C solver for two-phase flow (water and air) in a rocking cell-culture bag, with oxygen
transfer and passive tracers. It began as the driver code of Kim et al. (2025),
<https://github.com/rcsc-group/BioReactor>.

Build with `make build` from the repository root (other variants: see the `Makefile`). Run with
`BioReactor params.json [restart.dump]` (parameters: `docs_site/reference/params.md`).

Read [`docs_site/developer/solver-pitfalls.md`](../docs_site/developer/solver-pitfalls.md) before you
change restart code, face fields, or events.

## Files

| File | What it holds |
|---|---|
| `BioReactor.c` | Start here. Compile switches, global state, `main()` (setup only), and the includes that set the event order. |
| `init.h` | `event init`: the bag geometry and initial fields, or restore from a checkpoint (with every repair a restore needs). |
| `scalars.h` | Release of the tracer (`event tracer`) and the oxygen reservoir (`event oxygen`). |
| `forcing.h` | The rocking motion as body forces (`event acceleration`), steady-streaming vorticity, droplet removal. |
| `snapshots.h` | Fields on a uniform grid for Python: snapshot windows and the velocity at the checkpoint. |
| `checkpoint.h` | `checkpoint.dump` at the end, and an early dump before the wall-clock limit. |
| `diagnostics.h` | Time series: `logstats.dat`, `normf.dat`, `vol_frac_interf.dat`, `tr_oxy.dat`, `shear_stress.dat`, `edr_pct.dat`. |
| `restart.h` | Restart checks: `assert_fs_sane()`, `restart_diagnostic_*.txt`. |
| `replay.h`, `replay_events.h` | Flux record/replay (`-DRECORD_FLUX`, `-DREPLAY_FLUX`); file format in `replay.h`. |
| `movies.h` | Video frames (`-DVIDEOS=1`), rendered by `scripts/figures/render_videos.py`. |
| `henry_oxy2.h` | Two-phase scalars with a Henry's-law jump: advection and diffusion (from Kim et al., with MPI/restart fixes). |
| `utils2.h` | `statsf2()`: statistics over the liquid (Kim et al., unchanged). |
| `params_read.h`, `jsmn.h` | `params.json` reader (jsmn, MIT licence). |

## One time step

Basilisk calls the events of a step in the order it registers them. User events run first, in the
order of the includes at the end of `BioReactor.c`; then the library's `last` events (VOF advection,
scalar diffusion, projection).

| Order | Event | Header | When |
|---|---|---|---|
| 1 | `init` | `init.h` | t = 0 (also on restart) |
| 2 | `tracer`, `oxygen` | `scalars.h` | from t_mix (oxygen: every step) |
| 3 | `streaming_vorticity`, `acceleration`, `remove_drop` | `forcing.h` | every step |
| 4 | `field_snapshots` | `snapshots.h` | snapshot window, video cadence |
| 5 | `maxruntime_stop`, `dump_checkpoint` | `checkpoint.h` | every 10 steps; t_end (ends the run) |
| 6 | `logstats`, `pressure_diagnostics`, `normcal` | `diagnostics.h` | every 0.1; every t_out = 0.02 |
| 7 | `rec_checkpoint`, `vof` (record/replay hook) | `replay_events.h` | record/replay builds only |
| 8 | `movies_output`, `movies_output_tau` | `movies.h` | `-DVIDEOS=1` only |
| last | VOF, diffusion, projection | Basilisk, `henry_oxy2.h` | every step |

## Compile switches

Set with `-DNAME=value` (the `Makefile` targets do this). Defaults in **bold**.

| Switch | Default | What it does | Used by |
|---|---|---|---|
| `EXTRA_TRACERS` | **1** | 0 drops tracers c, c1, c3 (about 40% faster) | `make build-mpi-lean`, `build-mpi-3d` |
| `VIDEOS` | **0** | write video frames | `make build-video`, `build-mpi-video` |
| `DIAGNOSTICS` | **0** | write `pressure_diag.dat` | `make build-health` |
| `RECORD_FLUX` / `REPLAY_FLUX` | **0** | flux record / replay (see `BioReactor.c` top) | `scripts/replay_gate.py`, `replay_lf.py` |
| `CROSS_LEVEL_WARMSTART` | **0** | restart from a coarser-level checkpoint | `make build-crosslevel`, `tests/verification` |
| `CROSS_LEVEL_PREREFINE_SOLID` | **1** | 0 = old, wrong order (A/B tests only) | — |
| `N_RAMP_CYCLES` | **3** | length of the forcing ramp, in periods | slow-ramp tests |
| `STREAMING_VORTICITY` | **1** | steady-streaming vorticity output | — |
| `SEED_KIM` | **1** | tracer seeded in full liquid cells only (0: also f > 0.5) | `scripts/mixfix_check.py` |
| `LEGACY_SU` | **0** | old restart velocity rescale (A/B only) | `scripts/analyse/ab_test_su_rescale.py` |
| `REC_KIM_FLUX`, `REC_KIM_MASK` | **0** | record a rebuilt flux, to measure that approach's error | `scripts/analyse/replay_kim.py` |
| `PREREG_DIAG` | **0** | per-step sums, same code as on the Kim et al. side | `scripts/prereg_offline.py`, `scripts/kim_prereg/` |
| `REC_CHECKSUM` | **0** | per-step weighted sums in record and replay, to find where they part | `scripts/replay_gate_f3/` (flags on the build command line) |

Removed switches (finished experiments; in git tag `archive/pre-lean-2026-11`): `KICK_TEST` (velocity-kick
stability probe, `scripts/submit_kick_test.py`), `G_RESTART_MODE` (modes 0/1 of the g reset on restart).

## Changing the solver

1. Make the change on a branch.
2. If it should not change results, run the golden-run gate on a compute node:
   `sbatch tools/golden/golden.sbatch <src dir> <work dir>` for the base and for your branch, then
   `python tools/golden/golden.py compare base/manifest.json new/manifest.json --drop-unstable base2/manifest.json`.
   Also `... golden.sbatch <src dir> <work dir> compile` builds every research switch.
3. If it should change results, say which ones and show it with a test (`tests/verification/`).
