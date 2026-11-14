# Solver pitfalls

Traps in Basilisk and in this solver that fail **silently**: the run finishes, the numbers look
plausible, and they are wrong. Read this page before you change restart code, face fields, events or
the statistics. Each section has an anchor; comments in `src/` point here as
`pitfall "<anchor>"`.

To check that a change to `src/` does not change results, use the golden-run gate
(`tools/golden/golden.py`, see [Solver source guide](solver-guide.md)).

---

## bare-return

**Never write a bare `return;` inside an event.** qcc compiles an event body into a function that
returns `int`, and a nonzero value tells the event loop to stop. A valueless `return;` returns an
undefined value. The run can end early and still print "Simulation complete".

Use `return 0;` to skip the rest of the event, or wrap the body in an `if`. Use `return 1;` only to
stop the run on purpose (as `event dump_checkpoint` does).

## event-order

Basilisk runs the actions of one step in the order in which it **registers** the events, that is, the
order in which it meets them in the source. In `BioReactor.c` the `#include` lines at the end of the
file set this order. Moving an include can change results. Events marked `last` (the VOF, diffusion
and projection events of the library) run after all user events of a step.

Actions that share one event name run in **reverse** order of registration. `replay_events.h` uses this
to run its `event vof` before the VOF step of the library.

## fs-after-restore

`dump()` does not write **face fields** (`fs`, `uf`, any `face vector`): `dump_list()` in
`src/output.h` skips them. After `restore()` a face field is not zero either: it is partly filled.
Measured on a level-6 restart: about half of the fully open faces (both neighbours `cs > 0.99`) had
`fs == 0`. `cs` (a plain scalar) comes back exactly.

What goes wrong: `refine_embed_linear` chooses its interpolation from the coarse `fs`. Where `fs` is 0 it
falls back to piecewise-constant values. During a refine after restore (a warm start at a finer level),
half of the cells get that fallback: a checkerboard pattern, about +30% on shear stress and +200% on
dissipation, with no error message.

Do this after every `restore()`:

1. Rebuild `cs`/`fs` with `solid()` (in `init.h`: `bag_solid()`).
2. Do it **before** any `refine()` or `adapt_wavelet()`.
3. After `solid()`, call `restriction ({cs, fs})`: `solid()` writes leaf cells only, and the coarse
   levels feed the multigrid solvers every step. (Basilisk's own pattern: `src/test/neumann3D.c`.)

`assert_fs_sane()` (in `restart.h`) stops the run if a fully open face has `fs == 0`.

## attributes-after-restore

`restore()` does not restore the **function attributes** of fields: `refine`, `prolongation`,
`restriction`, `embed_gradient`. The library sets them in `event defaults (i = 0)` and
`event metric (i = 0)`, and a restart never has `i = 0`. Without the embed-aware versions:

- the scalar multigrid spreads NaN from solid cells into fluid cells;
- the pressure gradient at the wall uses the generic gradient, and the run fails (SIGFPE) about one
  rocking period after the restart.

`reapply_field_attributes()` in `init.h` sets them again, before any `refine()`.

## restore-nodump

`p` and `pf` have `nodump = true` by default. `restore()` skips every field that is `nodump` **in the
reading process**, whatever the writing process did. Set `p.nodump = pf.nodump = false` around both
`dump()` and `restore()`, or the pressure comes back as exactly zero.

## restore-event-time

`restore()` marks every event at the restored time as done without running it. An event at exactly
that time (for example a release at `t = t_mix` when the replay starts at `t_mix`) never runs. The
replay therefore releases by a time test in an `i++` event.

## event-loop-end

`run()` stops when no event with an end condition is left. Events with no end (`i++` with no
condition) keep the loop going for ever. The coupled run ends with `return 1` from
`event dump_checkpoint`; the replay needs `event end_replay (i++; i < 2000000000)` for the same reason.

## restart-phase

The time unit `T_bio` depends on `omega_b` **and** `theta_max`. A checkpoint written at a whole number
of periods of one run is not at a whole number of periods of a run with other parameter values. The
restored fields are at phase 0, so the forcing must also be at phase 0. Since 2026-11-13 the dump's label
(`checkpoint.phase`, see restart-ramp) carries the forcing phase in periods (`cycles`), and the forcing uses
the integrated phase `psi(t)` from there (`src/forcing_state.h`). Do not change `t` itself (the scheduler
then fails in `dtnext()`).

Compute the offset from the restored time alone (the nearest period boundary), not from the ratio of
the old and new periods. The ratio is 1 for the second link of a chain that changes nothing, but that
link's time already carries the first link's offset, so it would resume off phase. Errors of tens of
degrees put the flow on a different quasi-steady state (about +25% on mean shear stress).

An early (wall-clock) dump is not on a period boundary; its label holds its phase. The solver reads the
label from its **working directory**, not from the dump's folder. A resume that runs in another folder must
copy it there first (`scripts/submit_continue.stage_phase`; `simulate.submit_slurm` copies it for MPI
jobs). Since 2026-11-13 a missing label stops the restart; before, the solver silently took the nearest
period boundary. Measured 2026-11-13 (feature bundle): a dump
at 1.84 periods was resumed 0.16 period off phase; velocity +26 %, and the run stopped at 3.84 periods
instead of 3.

## restart-ramp

**Symptom:** the flow jumps at the seam of a restart. In `fig9_l9_rpm37.5_ext1` the liquid velocity per
period went 0.64, 0.64 | 0.64, 0.69, **0.96**, 0.69, 0.68: +49 % three periods after the seam.

**Cause:** every restart starts a new forcing ramp (`t_ramp_start = t_checkpoint` in `src/BioReactor.c`):
over 3 periods the rocking goes from the `*_prev` values (`omega_b_prev`, `theta_max_prev`, ...) to the
current ones. That is how a sweep changes from one condition to the next smoothly. But when `*_prev` is
not in `params.json`, it is 0, so a plain continuation ramps from **zero** rocking while the liquid still
sloshes at full amplitude.

Found again on 2026-11-13 in all 8 L8 -> L9 campaign warm starts (`campaign.py` wrote no `*_prev`): the
L9 velocity peaked +40 % and was still +10..+17 % at the tracer release; the BO objective dtmix_0.25 came
out 14.2 instead of 19.9 periods at 34.8 rpm (diary 2026-11-13).

**Rule since 2026-11-13 (binaries built from this `src/`):** the dump carries its source. Every dump writes
a label, `checkpoint.phase`: the forcing at that instant (`omega_b`, `theta_max`, ... ; inside a ramp, the
ramped values), its phase in periods, and the velocity unit of the fields (`U_bio`). A restart ramps from
the label to `params.json`, in amplitude, phase **and frequency** (`src/forcing_state.h`). Same forcing
(continuation, cross-level warm start) = no ramp, with nothing to set. `params.json` holds the target only:

- no label next to the run -> the solver stops (ERROR, exit code 1);
- a label with the forcing **and** `*_prev` keys -> stops (two sources);
- an old label (phase only, written before 2026-11-13) -> needs `*_prev` (allowed, with a WARNING);
- `restart_forcing.json` records the angle and its rate at the restart; the bundle guard S13 checks that
  they equal the label's (no jump).

**Older binaries** (all production binaries before the lean branch is merged) still need `*_prev` = the
source forcing; without it they ramp from zero. `scripts/simulate.binary_reads_label` tells the two apart,
`simulate.source_forcing_keys` writes `*_prev` only for an old binary, and `submit_slurm` refuses a restart
that breaks either rule. Use `scripts/submit_continue.build_params` or `scripts/chain.py`; never write these
keys by hand. To check old runs: a restart whose `params.json` has no `theta_max_prev` and that ran on an
old binary has this defect.

Corollary for tests: a cold run is inside its own ramp for its first 3 periods. Compare a restarted run
with a cold run only after the cold run's ramp.

## qcc-dimension-check

**Symptom:** `qcc` fails with exit code 2 and prints nothing but the usual `stencils` warnings. The same
source can build in one folder and fail in another.

**Cause:** with dimensional analysis on (the default), `qcc` runs `main()` and the events in an interpreter
before it accepts the build (`ast/interpreter/`). There `params.json` is never read, `fopen` returns no file
and `fread` an undefined value (`ast/interpreter/overload.h`). Every condition on an undefined value splits
the interpreter's run, and when the run uses up its call budget (`-maxcalls`, default 20 M) the build fails
(`qcc.c`: "status = 2; // dimensional error"), with no message unless you pass `-run=2`. What costs the
budget: loops bounded by a parameter (`k <= params.n_harmonics`), loops over parsed tokens, iterative solves,
and any branch-heavy function called from `main()` or called again from a late event. Found 2026-11-13 with
the first version of `src/forcing_state.h` (a 200-step bisection over the forcing; a jsmn token loop; the
forcing recomputed in the checkpoint event).

**Rule:** to see why a build fails, add `-run=2` and look for "reached maximum number of interpreter calls".
In code that `main()` or the events run: bound loops by constants (`N_MAX`) and switch terms off with 0/1
weights; prefer `fmin`/`fmax` to `if`/`?:`; compute the forcing once per step (the acceleration event) and
reuse it; never call `exit()` on a path the interpreter can reach (a missing file is such a path): set a
flag and `return` from `main()`. Then build in two different folders before you trust the result.
A larger `-maxcalls` does not help (100 M failed like 20 M: the split paths grow too fast), and a cache test on
`t` (`if (t_cached != t)`) splits the path too, because `t` is undefined there (both seen 2026-11-13).

**Corollary, the forcing cache is one step old at a dump:** the dump events run before the acceleration
event in a step, and after the last step no acceleration runs. So `_forcing_last` holds `t - dt` when the
label is written. The label takes the phase at `t` from the light `psi_at()` and writes `Th`, `Th_d`,
`Th_2d` with their own time `t_forcing`. Before this fix every restart lagged one step in phase and stopped
one step late (guard S12).

## restore-first-step

**Symptom:** after a restart, the first step that applies the rocking force is not at the dump time
`t_checkpoint` but one short step later (`restart_forcing.json`: `t` > the label's `t`).

**Cause:** `restore()` (Basilisk `src/output.h`) brings back `t` and the step number `i`, then runs the
events with `events (false)` to catch them up. Its last call is at the restored step, so it marks every
`(i++)` event of that step as done. The time loop then goes on after `init` in the same pass: only
`properties`, `stability` and events with their own condition run; `vof`, advection, `acceleration` and
`projection` do not. `stability` sets `dt`, the loop moves `t` by it, and the fields stay as dumped. The step
size starts small after a restore (here 6e-4 of a period at L6, then x1.5 per step), so the phase error is
small. Seen with `qcc -events` (event trace) on the bundle's `cont` case, 2026-11-13.

**Rule:** do not expect a restarted run to repeat the unbroken run bit for bit. Guard S13 allows the first
forcing step to come at most 1e-3 period after the dump.

## record-landing

A record build (`-DRECORD_FLUX=1`) has two empty events (`src/replay.h`) whose only job is to make a time
step end exactly on each end of the record window. Basilisk shortens the step to land on any event time,
so after the window ends the record build takes other time steps than a default build. Both are correct;
they are not bitwise equal after the window (≈0.1 % at the next output). Compare a record build with a
default build only up to the end of its record window; for longer comparisons use two runs of the same
build.

## ramp-derivatives

When a restart ramps the amplitude from an old value to a new one, the angle is
`Th = A(t) sin(...)` with `A` changing in time. The angular velocity and acceleration in the body force
must include the `A'(t)` and `A''(t)` terms. Without them the forcing does not match the angle during
the ramp, and the restarted run settles in a different state from a cold start.
Since 2026-11-13 the frequency ramps too: the phase is the integral of `w(t)`, never `w(t) t` (that gives
a chirp, a rate `w + w' t`), and the angular acceleration has the extra term `A k w'(t) cos(...)`.

## scalar-continuation

On a continuation (`restart_continue = 1`) the tracer and oxygen fields are the experiment: do not reset
them and do not release the tracer again. On every restart, call `boundary()` and `restriction()` on the
scalars: coarse levels and MPI ghost cells keep stale values from the checkpoint otherwise, and the
multigrid diverges slowly.

## mask-outside-bag

`f` is initialised over the whole box, also outside the embedded bag, where the velocity is 0. A
"liquid" mask `f > 0.5` alone counts those cells in the area of every spatial mean (3.5x too much area at
level 10). Use `cs > 0 && f > 0.5`, and weight cut cells by their fluid area `cs*Delta^2`. Maxima are not
affected, which makes the error easy to miss.

## volume-metric

Check a fix with a **volume** quantity, not only a wall or surface quantity. Grid-scale noise spread
through the domain hardly moves the wall shear stress but strongly moves the viscous dissipation
`~|grad u|^2`. Example: the same defect gave +30% on shear stress and +200% on dissipation.

## openmp-arrays

Inside `foreach()`, qcc adds OpenMP reductions for scalar accumulators, not for arrays. `bins[b]++` is a
data race under OpenMP (two identical runs gave percentiles 38% apart). An `omp` pragma inside
`foreach()` does not compile. Give each thread its own bins and add them after the loop
(`merge_bins()` in `diagnostics.h`). MPI builds have no OpenMP, so they are not affected.

## frame-cadence

A sample interval that divides the rocking period exactly (`T/N`) samples the same `N` phases for ever,
and misses the peaks. The video and snapshot interval is `T/(N + 0.618...)`: successive frames fall at
evenly spread phases.

## stale-header-copies

A local copy of a Basilisk header drifts from the installed Basilisk API. If a local `.h` file does not
compile, check `git log` for it first. If it was only imported, compare it function by function with
the canonical file, replace it with the canonical file, and add back only the project additions. Do not
patch it line by line.
