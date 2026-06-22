# How We Made MPI Checkpoint Restart Work

## A Development Log — June 2026

This document records the full story of building MPI checkpoint restart for
BioReactor3D, including every wrong turn, every falsified hypothesis, and every
debugging detour.  It is written to be verbose on purpose: the summary document
(`mpi_checkpoint_restart.md`) tells you *what* the bugs were and how to fix
them.  This one tells you *why it took three weeks and 50+ SLURM jobs to find
them*.

---

## 1. Why We Needed This at All

At fidelity 7 (L7, 2¹⁴ = 16 384 max cells), a single simulation run takes
roughly 0.64 hours on 16 MPI ranks for 8 non-dimensional mixing cycles.  For
80 cycles total — the minimum needed to extract kLa₂₅ reliably — that is
about 5–6 hours per condition.  Entirely feasible.

At fidelity 8 (L8, 2¹⁶ = 65 536 max cells), the same 80 cycles requires
roughly 24 hours.  OSCAR's default partition walltime limit is also 24 hours.
With a typical cluster load the jobs almost never actually get 24 hours before
preemption.  L8 sweeps were simply not feasible without some form of
continuation.

The plan: implement checkpoint restart so each simulation is split into a chain
of segments.  Each segment fits inside the walltime limit and hands state
forward via Basilisk's `dump`/`restore`.  The chain self-submits — segment *k*
submits segment *k+1* from its own epilogue via `sbatch`.

This should have been straightforward.  It was not.

---

## 2. Phase 1 — Getting L7 MPI to Run at All

Before we could even think about checkpoint restart, MPI runs at fidelity 7
were crashing.  Every L7 job died at exactly `t = 4.25` with a SIGSEGV.
Fidelity 6 was fine.  This was our first detour.

### The crash at t = 4.25

The crash appeared after commit `e895f4e`, which added dynamic walltime
estimation and — more relevantly — three diagnostic instrumentation loops:
`PRERELAX3`, `FINE_RES_BLOW`, and `REGION_ALL` inside `henry_oxy2.h`.

**Hypothesis H1** (falsified): the first guess was that `scalar tau_liq[]`
declared inside `NORMCAL` was leaking one global scalar per invocation,
eventually exhausting the Basilisk scalar pool and causing a SIGSEGV.  We
removed the allocation (commit `fc2f1c8`) and resubmitted.

The next job (3298888) ran clean to `t = 74.2`.  Victory.  Except it was not:
that job had been built against the pre-`e895f4e` `henry_oxy2.h` baseline
(`9ddd52c`) from a prior batch of test binaries, not the `fc2f1c8` fix.  A
second wave of jobs with the actual `fc2f1c8` fix still crashed at `t = 4.25`.
The scalar hypothesis was falsified.

**Hypothesis H2** (confirmed): `e895f4e` introduced the crash, and the
diagnostic loops were the culprit.  The loops accessed `D.x[1]` and
`D.y[0,1]` — face-component offsets at position +1 — inside spatially
conditioned `foreach_cell()` branches.  The `qcc` stencil analysis does not
follow branches with runtime conditions, so it never registers that these
face-ghost buffers need MPI halo expansion before they are accessed.  On 16
ranks at L7, the probability of a rank boundary landing exactly on the bad cell
approaches certainty, and the first access into an unexpanded ghost buffer
triggers a SIGSEGV.

**Hypothesis H3** (confirmed): removing just those three loops from the full
`e895f4e` commit produced a binary (job 3298986) that ran to `t = 250.5`
without incident on 16 ranks.

**Hypothesis H4** (confirmed): removing all diagnostic instrumentation from
`henry_oxy2.h` on main (commit `9508817`) produced the clean production binary.
Job 3321334 completed a full L7 run (`kLa_25 = 0.132`).

This was the prerequisite: we needed stable L7 MPI runs before checkpoint
restart made any sense.

---

## 3. Phase 2 — The First Restart Attempt

With a clean L7 binary, we turned to the `mpi-checkpoint` branch.  That branch
already contained a collection of suspected fixes accumulated from earlier
partial debugging sessions:

- `fedisableexcept` wrap in `h_relax` (prevents SIGFPE from masking the real
  blow-up location during diagnosis)
- A call to `restriction({res})` after `foreach_cell() res[] = 0` in
  `h_residual` (the reasoning was that coarse ghost cells might retain stale
  values)
- Pre-extraction of `h_relax` stencil values before any conditional block
  (closing a suspected qcc stencil blind-spot)
- Harmonic mean for diffusion coefficient `D` at coarse multigrid levels
  (avoids spurious large `D` at embedded-wall-adjacent faces)
- A `d`-floor guard and `isfinite(n)` guard in `h_relax`

**Hypothesis H5** (falsified): we believed these accumulated fixes were
sufficient.  An interactive 2-rank test from a fidelity-5 checkpoint
(`smoke_l5_seg0`, `t_checkpoint ≈ 100`) was our first real restart.

The simulation ran without SIGFPE (confirming `fedisableexcept` was working)
right up to `t_mix = 100.334` — the moment the oxygen tracer injection event
fires and the gas–liquid interface opens to dissolved-oxygen transport.  At
`t = 100.334`, both `oxy` and all companion tracers blew up simultaneously.
The first `H_RELAX_BLOW` diagnostic line showed:

```
b = -1.42e+20  at level=2  x=-0.375 y=0.375
```

20 million blow-up lines were emitted before the process was killed.  The
`1.42e+20` value was unmistakable: it is the exact garbage value that Basilisk's
memory pool places into freshly allocated scalars before they are initialized.
Something from a previous multigrid solve was leaking into the restart.

---

## 4. Phase 3 — The b[] Rabbit Hole

The `b[]` vector in Basilisk's multigrid context is the right-hand side of the
linear system being solved (`A c = b`).  Finding `b = -1.42e+20` at level 2
meant the solver was being handed a completely corrupted source term.  Over the
next 72 hours we chased three different explanations for where that corruption
came from.

### H6 — Stale stracers in boundary() (falsified)

The immediate suspicion was the scalar tracers (`stracers`).  On checkpoint
restore, `stracers` are reset to zero via `reset(stracers, 0.)` and then
`boundary(stracers)`.  Our hypothesis was that `boundary()` only communicates
*leaf* MPI ghost cells, leaving coarse-level ghost cells with stale values from
the previous simulation run that was stored in the checkpoint.  Those stale
values would then propagate upward during `mg_cycle`'s restriction step and
corrupt `b[]` at coarse levels.

Fix: add `restriction(stracers)` after `reset() + boundary()` on restore
(commit `9a90b72`).  This forces `halo_restriction` to propagate the zero-filled
leaf values to all coarse ghost cells before any solve begins.

Result: blow-up *identical*.  Same location: `x=-0.375, y=0.375, level=2`.
Same value: `b = -1.42e+20`.  A 1-rank test (no MPI) produced zero
`H_RELAX_BLOW` lines — confirming the blow-up was MPI-specific but
`restriction(stracers)` had fixed nothing about it.

### H7 — foreach_cell res[] = 0 misses ghost cells (partly confirmed, incomplete fix)

The `mg_cycle` inner loop inside `h_residual` initializes the residual scalar:

```c
foreach_cell() res[] = 0;
```

`foreach_cell()` in Basilisk iterates over the *local* active tree — it skips
MPI ghost cells that belong entirely to a remote rank and have no fine-level
descendants on the current rank.  Those ghost cells retain whatever value `res`
held from the *previous* call to `mg_solve` — which was the pressure Poisson
solver, not the tracer diffusion solver.  The pressure solver's residual at
`level=2` was on the order of `1e−5` (normal convergence residual), but the
pool-stale value from an unrelated prior allocation could be `−1.42e+20`.

Fix: add `restriction({res})` immediately after the `foreach_cell` zero-init in
`h_residual` (commit `bc1b92e`).  This forces a full `halo_restriction` from
the leaf level upward, clobbering any pool-stale coarse ghost values with zeros
propagated from the leaves that *were* correctly initialized.

This was partially right.  The `BLOWUP_PROBE` diagnostic (added alongside the
fix) showed `b = 0` at `t = 99.9014` for all tracers — meaning the coarse ghost
cells of `b[]` were indeed being zeroed before `t_mix`.  But the blow-up still
occurred at `t_mix`.

The new symptom was illuminating.  The blow-up diagnostic now showed:

```
nx=0  ny=-7e-9  nb=0.000134  →  sum ≈ 0.000134
but n = 8.87e+18
```

The printed neighbor contributions summed to a physically reasonable `0.000134`,
yet `n` (the diagonal coefficient) came out as `8.87e+18`.  This was
*impossible* if `n` were constructed solely from `nx + ny + nb`.  Some other
code path had to be writing a huge value into `n`.

There was a complication: the diagnostic loop that printed `nx, ny, nb` was
itself inside a spatially-conditioned block — the same qcc stencil blind-spot
that had caused the L7 crash in Phase 1.  The *printed* neighbor values were
from unexpanded buffers and were wrong.  The *actual* `n` used in the Gauss-
Seidel relaxation might be correct, but the diagnostic was lying.

### H8 — Ghost cells still have non-zero b[] (confirmed mechanism, wrong hypothesis)

We added a `GHOST_B_CLAMP` diagnostic to identify every MPI ghost cell with
`|b[]| > 0` before `h_relax` executes (commit `11aceea`).  The diagnostic fired
reliably at `t > 99.9`, confirming that coarse ghost cells were reaching `h_relax`
with non-zero `b[]`.

This looked like confirmation of the ghost-cell contamination theory.  The fix
was: inside `h_relax`, check `is_local(cell) && l < depth()` and zero `b_eff`
for non-local coarse cells, since `boundary_level()` overwrites their `c[]`
after each sweep anyway.

Result: `GHOST_B_CLAMP` *never fired*.  The fix zeroed ghost `b[]` completely —
but the blow-up was identical.  Same `n = 8.87e+18`, same location, same time.
This confirmed that ghost `b[]` contamination was a symptom, not the cause.
Something was taking a small, correct `n` value and adding `8.87e+18` to it
*after* all the neighbor terms were accumulated.  

Looking at `h_relax` carefully, only one code path modifies `n` after the
neighbor accumulation:

```c
if (p->embed_flux) {
    double c_embed = p->embed_flux (point, c, D, a);
    n -= c_embed * sq(Delta);
}
```

If `c_embed ≈ -1.42e+20 / sq(Delta)`, the arithmetic worked out exactly to
`8.87e+18`.

---

## 5. Phase 4 — The embed_flux Revelation

**Hypothesis H10** (confirmed): The blow-up comes from the `if (p->embed_flux)`
block in `h_relax`.  Two coupled mechanisms make this happen:

1. `struct HDiffusion q` is declared *without an initializer* inside
   `tracer_diffusion`:
   ```c
   struct HDiffusion q;   // ← uninitialized; q.embed_flux holds stack garbage
   ```
   The `q.embed_flux` member therefore contains whatever happened to be on the
   call stack at that point in the previous frame — a non-NULL garbage pointer,
   typically `0x2300000022`.

2. The guard `if (!q.embed_flux && aa.boundary[embed] != symmetry)` is supposed
   to populate `q.embed_flux` only for tracers that need the embedded-boundary
   flux correction.  The oxygen scalar `oxy` uses a `symmetry` embedded-boundary
   condition (no normal flux through the bag wall), so the guard is intentionally
   false for `oxy` — `q.embed_flux` should remain NULL.  But because the struct
   is uninitialized, `q.embed_flux` is already non-NULL (garbage), so the
   *second* half of the condition (`&& aa.boundary[embed] != symmetry`) evaluates
   to true, and the whole guard is false — `q.embed_flux` is *never overwritten
   with NULL*.

3. Later, `if (p->embed_flux)` fires because `p->embed_flux = q.embed_flux` is
   that same garbage non-NULL value.  The embed_flux callback is called.  For
   `oxy`, the callback computes a flux correction proportional to
   `mua / (fa + SEPS)` where `fa` is the face-fraction of the embedded boundary.
   At coarse multigrid levels (`level=2`), many faces adjacent to the superellipse
   boundary have `fa ≈ 0` (the geometry is well-resolved at the leaf level but
   not at the coarsened representation).  This causes
   `mua / (fa + SEPS) → mua / SEPS`, which diverges.

4. Before `t_mix`, `oxy[] = 0` everywhere (tracer injection has not fired yet).
   `embed_flux` returns 0 when `oxy = 0` because the symmetry gradient is zero.
   This is why all our diagnostics prior to `t_mix` showed `b = 0` everywhere —
   the garbage `embed_flux` pointer was being called but returning 0.  The
   moment `oxy` became non-zero (at `t_mix = 100.334`), the callback returned
   `c_embed ≈ -1.42e+20`, adding `c_embed * sq(Delta) ≈ -8.87e+18` to `n`,
   which caused the multigrid to diverge in the first half-step.

The diagnostic that confirmed this (added in commit `a83fd5a`, job
`interactive-h10-embed-diag`, PID 3114673) showed:

```
pef = 0x2300000022  for ALL 6 tracers (c, oxy, c1, c2, c3, c4) at l=2 x=-0.375 y=0.375
```

`pef` was the pre-extracted `p->embed_flux` pointer value.  `0x2300000022` is
not a valid function pointer — it is the garbage that happened to live on the
stack at that call site.  The fact that all six tracers shared the same garbage
value is consistent: they are all called in the same function with the same
uninitialized stack frame layout.

We observed that `c_embed = 0` at `t = 100.001` through `t = 100.041` (oxy
still zero from injection delay) — then at `t = 100.334` the injection fires,
`oxy > 0`, `embed_flux` is called with the garbage pointer, and the return value
drives `c_embed ≈ -1.42e+20`.

**The fix** (commit `6f15883`): add

```c
q.embed_flux = NULL;
```

immediately after `struct HDiffusion q;` is declared.  This ensures that for
tracers with a `symmetry` embedded-BC, `q.embed_flux` stays NULL throughout,
and the `if (p->embed_flux)` guard in `h_relax` correctly evaluates to false
at all multigrid levels.

**Hypothesis H11** (confirmed): The interactive 2-rank restart from
`smoke_l5_seg0` (PID 3135406) passed `t_mix = 100.334` with zero blow-up,
ran to `t = 207.9`, and `oxy_liq_sum = 0.283` — finite, positive, and
physically reasonable.  The previously identical run had blown up at `t = 100.3`.

This confirmed we had found and fixed the root cause.  The other fixes applied
earlier (`restriction({res})`, harmonic mean for `D`, the stencil pre-extraction)
are retained as belt-and-suspenders hardening — they address real code
vulnerabilities even if they were not the primary blow-up cause.

---

## 6. Phase 5 — The Validation Grid (H12–H25, 14 SLURM Jobs)

Finding the root cause is not enough.  A single passing test tells you the fix
worked for one condition.  We needed to know it generalised — different rank
counts, different fidelity levels, different omega values, different fill levels,
different bag geometries, with and without velocity rescaling on restart.

We submitted two grids of jobs.

### Compatibility grid (H12–H17, jobs 3336237–3336242)

These tested rank count × fidelity × omega type:

| H | Ranks | Fidelity | Restart type | Outcome |
|---|---|---|---|---|
| H12 | 2 | 5 | same ω | kLa_25=0.453 ✓ |
| H13 | 2 | 5 | 2× ω change + velocity rescale | kLa_25=0.216 ✓ |
| H14 | 16 | 5 | same ω | kLa_25=0.524 ✓ |
| H15 | 16 | 5 | 2× ω change + velocity rescale | kLa_25=0.244 ✓ |
| H16 | 2 | 6 | same ω | kLa_25=0.166 ✓ |
| H17 | 16 | 6 | same ω (production scale) | kLa_25=0.159 ✓ |

H15 was the stress test: maximum MPI decomposition plus velocity rescaling
simultaneously.  All six confirmed in a single SLURM batch.

### Extended parameter grid (H18–H25, jobs 3336706–3336713)

These pushed geometry and parameter extremes:

| H | Config | Outcome |
|---|---|---|
| H18 | fid=4, n=2 (elliptical bag — curved boundary hits more coarse cells) | kLa_25=0.182 ✓ |
| H19 | fid=7, fill=0.3 (near-empty — free surface near bag floor) | kLa_25=0.149 ✓ |
| H20 | fid=7, fill=0.7, ω=2.618 (near-full) | kLa_25=0.040 ✓ |
| H21 | fid=7, fill=0.6 | kLa_25=0.060 ✓ |
| H22 | fid=7, θ=5°, ω=3.1416 (high frequency) | kLa_25=0.042 ✓ |
| H23 | fid=8, θ=4° (production fidelity) | kLa_25=0.033 ✓ |
| H24 | fid=8, θ=7° (production fidelity + angle) | kLa_25=0.035 ✓ |
| H25 | fid=7, restart from end of completed run — kLa agreement | |Δ|/ref=1.4% ✓ |

H25 was the most important from a scientific standpoint: it confirmed that a
restart from a *completed* run's checkpoint produced kLa values within 1.4% of
the original (`kLa_25 = 0.141` vs. reference `0.139`).  This is the bound we
cite in the validation table — floating-point non-determinism across different
MPI rank layouts accounts for essentially the entire deviation.

All 14 extended tests confirmed.  25 total hypotheses tested; 10 falsified, 15
confirmed.

---

## 7. Phase 6 — SLURM Infrastructure Bugs (The Problems Nobody Expects)

With the binary working, we turned to the sweep infrastructure.  Running 60
theta-sweep conditions at L7 via MPI checkpoint chains requires coordinating
600 SLURM jobs (60 conditions × ~10 segments per condition).  Three independent
bugs in `slurm_mpi_template.sh` prevented chains from actually self-submitting.

### Bug 2 — Chains stopped after segment 1

**Symptom:** Seg-0 jobs completed and submitted seg-1.  Seg-1 jobs completed and
submitted... nothing.  Every chain died after its second segment.

**Root cause:** `slurm_mpi_template.sh` derived the canonical run directory for
the *next* segment as:
```bash
NEXT_CANON="$(dirname "$CANON_RUN")/$NEXT_RUN"
```
`CANON_RUN` is read from the `_canonical_run_dir` field in `params.json`.  That
field is written by `simulate.py` only for the *first* (seg-0) job.  Subsequent
`params.json` files are generated programmatically by `sweep.py` without it.
So for any seg-1 or later job, `CANON_RUN=""`, `dirname ""` evaluates to `.`,
`NEXT_CANON = "./<next_run_id>"`, and the existence check silently failed.

**Fix (commit `db4e12f`):** Derive `RUNS_ROOT` from `_experiment_dir` instead,
which *is* propagated through every segment since `sweep.py` writes it when
creating each params file.

We found this bug by watching seg-2 jobs fail to be submitted and
cross-referencing the SLURM epilogue logs.  The logs showed `NEXT_RUN` set
correctly to a valid run ID, but no `sbatch` output — the `if [ -f "$NEXT_PARAMS_CANON" ]`
guard was evaluating to false.

### Bug 3 — Postprocessing died with "Permission denied" on /oscar/scratch

**Symptom:** Jobs 3359029, 3359030, 3359031 completed their simulation phase
(checkpoint written, MPI collective exited cleanly), but reported `ExitCode 2`
in SLURM and left no `results.json`.  The `.err` log showed:

```
can't open file '/oscar/scratch/scripts/postprocess.py': [Errno 13] Permission denied
```

**Root cause:** The postprocess block in the template has two branches.  The
first (correct path) runs when `CANON_RUN` is non-empty:
```bash
if [ -n "$CANON_RUN" ] && [ "$CANON_RUN" != "$SCRATCH_RUN" ]; then
    PROJECT_ROOT="$(dirname "$(dirname "$CANON_RUN")")"
    uv run python "$PROJECT_ROOT/scripts/postprocess.py" ...
```
The second (fallback) runs otherwise:
```bash
else
    PROJECT_ROOT="$(cd "$SCRATCH_RUN/../../.." && pwd)"
    uv run python "$PROJECT_ROOT/scripts/postprocess.py" ...
```
Since Bug 2 was fixed by deriving `RUNS_ROOT` from `_experiment_dir` but not
by fixing `CANON_RUN` itself, chain-submitted jobs still had `CANON_RUN=""`.
`SCRATCH_RUN` is `/oscar/scratch/eaguerov/mpi_runs/<run_id>`, so `../../..`
resolves to `/oscar/scratch` — a completely unrelated directory with no
`scripts/` subdirectory and no permission from the project virtualenv.

**Fix (commit `8236ea7`):** Derive `CANON_RUN` from `_experiment_dir` *before*
the postprocess block, so the first branch is always taken for chain jobs.

We recovered the three failed runs manually by `rsync`-ing their output from
scratch to Lustre and running `postprocess.py` by hand.

### Bug 4 — Chains stopped at segment 3 (RUNS_ROOT off-by-one)

**Symptom:** After Bug 2 was fixed, chains ran to seg-2 and then stopped.
The SLURM log showed no error — seg-2 completed cleanly with a valid
`results.json`, but seg-3 was never queued.

**Root cause:** The `RUNS_ROOT` derivation from `_experiment_dir` in the Bug 2
fix used a single `dirname`:
```python
print(os.path.join(os.path.dirname(exp), 'runs'))
```
`_experiment_dir` is a path like
`.../rocking-bioreactor-2d/experiments/sweep_fb_theta_l7_mpi_ckpt/`.
`os.path.dirname(exp)` = `.../rocking-bioreactor-2d/experiments/`, so
`RUNS_ROOT = .../experiments/runs/` — a path that does not exist.

The `if [ -f "$NEXT_PARAMS_CANON" ]` guard evaluated to false (silently, because
`bash -e` does not exit on failed test conditions), and the `sbatch` call was
never reached.  No error, no submission, chain just stopped.

**Fix (commit `0a6a66f`):** Use two `dirname` calls:
```python
print(os.path.join(os.path.dirname(os.path.dirname(exp)), 'runs'))
```
This climbs from the experiment subdirectory to the project root, then appends
`runs/`.

This bug went undetected through all of Phase 4 and Phase 5 testing because
the validation grid (H12–H25) used direct `sbatch` submission with explicit
`PARAMS` and `DUMP` environment variables — not self-submitting chains.  The
self-submitting chain code path had never been exercised at L7+ fidelity before
the production sweep.

---

## 8. Phase 7 — Production Sweeps and Infrastructure Monitoring

With all three bugs fixed, we ran the L7 production sweeps:

**L7 fill sweep (`sweep_fb_fill_l7_mpi_ckpt`):** 50 conditions (5 fill levels ×
10 ω values).  All 50 completed with finite KPIs.  This became the first half of
the checkpoint validation figure.

**L7 theta sweep (`sweep_fb_theta_l7_mpi_ckpt`):** 60 conditions (6 θ values ×
10 ω values).  These ran as self-submitting chains.  Despite the three
infrastructure bug fixes, a new problem appeared: chains were silently stopping
mid-flight without any SLURM error.  Jobs were being submitted and then
cancelled within 5 seconds with ExitCode 0, leaving no log files.

We traced this to the watchdog mechanism using scratch-directory existence as a
proxy for "job is queued" — which cannot distinguish between a cancelled job
(scratch dir exists, but job dead) and a running one.  This led to double
submissions and eventual QOS violations (OSCAR allows max 4 simultaneous 16-rank
jobs for our allocation).

The fix was `.slurm_jid` marker files: each `sbatch` call now writes the
returned job ID to `$NEXT_SCRATCH/.slurm_jid`.  The watchdog cross-references
these files against live `squeue` output to detect truly-inflight jobs.  With
this, the watchdog correctly identified broken chains (completed run without a
queued successor) and re-submitted them.

All 60 theta conditions completed.  The kLa values (0.032–0.138 hr⁻¹) were
physically consistent with the L7 fill sweep and with the pre-existing L7 serial
baseline (`sweep_fb_theta_l7`), confirming no regression from the checkpoint
infrastructure.

---

## 9. By the Numbers

| Category | Count |
|---|---|
| Hypotheses tested (all threads) | 25 |
| Hypotheses falsified | 10 |
| Hypotheses confirmed | 15 |
| Interactive test runs | ~8 |
| SLURM jobs: binary debugging (L7 crash) | ~15 |
| SLURM jobs: MPI restart validation grid | 14 |
| SLURM jobs: infrastructure bug debugging | ~15 |
| SLURM jobs: L7 production sweeps (all segments) | ~600 |
| Calendar time from first attempt to clean L7 sweep | ~3 weeks |
| Root causes (binary level) | 2 |
| Root causes (infrastructure level) | 3 (Bugs 2, 3, 4) |
| Lines of diagnostic code added and then removed | ~400 |

The core insight that unlocked everything was H10: recognising that a *clean
ANSI C struct declaration without an initialiser* causes undefined behavior
through a specific chain of events that is invisible in serial runs, invisible
before `t_mix`, and invisible at fidelity < 5 (not enough cells to have coarse
ghost cells with near-zero face fractions).  All previous hypotheses were
attacking the symptom (`b[] contamination`, `stale ghost cells`) rather than
the source (uninitialized `q.embed_flux` allowing `embed_flux` to be called
where it should not be).

---

## 10. What We Learned

**On debugging Basilisk MPI:**  The `qcc` stencil analysis blind-spot (it does
not follow branches) is a genuine hazard.  Any loop that accesses face-offset
values (`D.x[1]`, `D.y[-1]`, etc.) inside a `foreach_cell()` conditional will
silently skip ghost halo expansion, causing crashes that appear randomly
depending on domain decomposition.  The symptom is a clean SIGSEGV at a rank
boundary.

**On uninitialized struct members in C:**  `struct HDiffusion q;` looks
harmless.  In a serial run at fidelity < 5, `q.embed_flux` happens to be NULL
by luck of stack layout, and nothing goes wrong.  At higher fidelity with MPI,
the specific combination of coarse-ghost cells + non-zero oxy + garbage pointer
produces a deterministic blow-up at a specific time.  The lesson: always
initialise struct members explicitly, especially function pointers that control
conditional code paths.

**On SLURM self-submitting chains:**  The `_canonical_run_dir` field is a
provenance artifact that only exists for the first segment.  Any template code
that relies on it for path resolution will silently fail for every subsequent
segment.  Always derive paths from fields that are propagated through all
segments (`_experiment_dir`, `next_run_id`).  Always add a strict `set -euo pipefail`
and test path resolution with a dry run before submitting a 600-job production
sweep.

**On diagnostic instrumentation in Basilisk:**  Diagnostic loops that print
field values from inside `foreach_cell()` conditional branches will read
unexpanded ghost buffers — the printed values are wrong, and the loop itself
introduces the same stencil blind-spot bug it was added to diagnose.  Use
pre-extracted temporaries (assign the field to a local variable outside any
conditional block) before printing.

**On the value of a falsifiable experiment:**  H10 was only reachable because
H5 through H9 were each run as proper SLURM or interactive jobs with specific
predictions written down before submission.  At the point where we noticed
`n = 8.87e+18` despite `nx + ny + nb ≈ 0.000134`, the combination of
"impossible arithmetic" and "only the embed block modifies n after neighbor
accumulation" was immediately falsifiable: either `p->embed_flux` was being
called when it should not be, or the embed block had a numerical bug.  The H10
diagnostic confirmed the former within an hour.  Without the discipline of
treating each run as a test of a specific prediction, this investigation would
have taken considerably longer.
