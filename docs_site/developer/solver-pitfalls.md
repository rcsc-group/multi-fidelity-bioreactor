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
restored fields are at phase 0, so the forcing must also be at phase 0: add an offset to `t` in the
forcing (`t_phase_offset`). Do not change `t` itself (the scheduler then fails in `dtnext()`).

Compute the offset from the restored time alone (the nearest period boundary), not from the ratio of
the old and new periods. The ratio is 1 for the second link of a chain that changes nothing, but that
link's time already carries the first link's offset, so it would resume off phase. Errors of tens of
degrees put the flow on a different quasi-steady state (about +25% on mean shear stress).

An early (wall-clock) dump is not on a period boundary; it writes its phase to `checkpoint.phase`.

## ramp-derivatives

When a restart ramps the amplitude from an old value to a new one, the angle is
`Th = A(t) sin(...)` with `A` changing in time. The angular velocity and acceleration in the body force
must include the `A'(t)` and `A''(t)` terms. Without them the forcing does not match the angle during
the ramp, and the restarted run settles in a different state from a cold start.

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
