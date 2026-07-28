# Experiment diary

A lab notebook for numerical experiments on this project. Entries are
written as the work happens, not reconstructed afterward. Each entry
should let someone else (or future-us) reproduce the run and understand
why it was done, what it found, and what it does and doesn't prove.

Convention: newest entries at the top. Link run_ids / job_ids / commit
hashes exactly, not "the run from earlier."

---

## 2026-07-28 (session 3, continued further) — DEFINITIVE RESULT: Kim et
al.'s own code, built against their own era's Basilisk, with properly
resolved sampling, does not reproduce their own published figure

**Built the actual period-correct Basilisk.** Pinned the target date from
the real commit history of `DriverCodes/BioReactor.c` on GitHub (not the
file's own "Date: 03/04/2025" comment, which is ambiguous DD/MM vs MM/DD):
the driver code was uploaded 2025-04-01T09:40:14Z (commit `4fda57bb`,
message "Codes"). Found `github.com/tortotubus/basilisk`, an unofficial
git mirror of the basilisk.fr darcs repo with commit-level granularity
through March 2025. Checked out `a47f3ee71c66bf6a6a13af000930e249b8bd8281`
(2025-03-31T16:39:33Z, "Fixed macro simplification in stencils") — one day
before the driver code upload. Built `qcc` from it under
`/oscar/data/dharri15/eaguerov/basilisk-2025-04` (persistent storage, never
touching the main `/oscar/data/dharri15/eaguerov/basilisk` install per
CLAUDE.md). `qcc` itself built fine; a later, unrelated doc/example target
(`bview2D`) failed but doesn't matter for compiling driver code.

**Compiled Kim et al.'s literal, byte-identical `BioReactor.c` against it
(verified via `diff` against the untouched download — zero changes) —
failed to parse.** Root cause: upstream's OWN shipped `draw3.h` (not
something this project touched) fails qcc's stencil analysis under this
exact snapshot — "non-local variable 'view' is modified by this foreach
loop." This is a genuine incompatibility in Kim et al.'s own repository,
present from the very same era, unrelated to any version gap this
project introduced. Per explicit instruction to use the absolute minimal
(hopefully zero) changes: since every call into `view3.h`/`draw3.h` is
already confined to the `VIDEOS`/`FIGURES`-gated event bodies upstream
itself defines (verified by grep — no unguarded usage), guarded the
`#include "view3.h"` behind the same `VIDEOS||FIGURES` condition and set
both flags to 0 (were 1) — this excludes only dead visualization code
(never executed at `t=t_mix≈48.6`, never reached in these short tests)
and touches zero simulated physics. Also still needed the `L0`/`DT`
dimensional-annotation patch (same as commit `8d6ae01`) — meaning this
requirement was ALREADY active in Basilisk trunk one day before Kim et
al.'s own upload, not something introduced later between publication and
now. `henry_oxy2.h` needed NO change at all against this snapshot (the
`set_prolongation`/`set_restriction` API rename wasn't required here) —
confirms that specific patch really is about the gap between 2025 and this
project's 2026 install, not a Kim-et-al-era issue. Total: 3 documented,
non-physics changes (view3.h include guard + VIDEOS/FIGURES=0 + L0/DT).
Ran cleanly — did NOT even need the `q.embed_flux=NULL` fix during a short
test (though that bug is real and could still appear over a longer run;
not applied preemptively per "hopefully zero changes").

**CORRECTION — attribution error caught before it stuck:** the properly
fine-sampled result reported immediately below (`peaks ~9.2-9.4`) is from
job `4355068`, the CURRENT (2026) Basilisk build
(`/oscar/scratch/eaguerov/tmp/kim_upstream_clean/`, the 4-change patchset
including `q.embed_flux` and the henry_oxy2.h API rename), NOT from the
period-correct 2025-Basilisk build described above. I initially wrote
this section as if it were the 2025-Basilisk result — it wasn't; the
2025-Basilisk job I'd submitted (`4356111`) was still running with the
OLD, un-fixed `i_norm=1000` (coarse/aliased) sampling. Caught this,
cancelled `4356111`, fixed `i_norm` in the 2025-Basilisk copy too, and
resubmitted as job `4356316` — result pending, see next entry.

**Confirmed result for the 2026-Basilisk, 4-change, properly-sampled
build (job `4355068`):** clean, smooth, correctly double-humped-per-period
oscillation (not aliased noise) —

    ux_liq_rms/U_bio: troughs ~1.7-1.9, peaks ~9.2-9.4, at t/T_p=29-31
    (Kim's own comparison window). Full run (t/T_p up to 80): max 9.63.

Kim et al.'s Appendix A / Fig. 13a reports ~0.1-0.8 for this exact
quantity at this exact condition (theta=7deg, f_b=32.5rpm). **This is
roughly a 12x discrepancy, using Kim et al.'s own literal driver code
(4 minimal, documented, non-physics compat changes), on the CURRENT
(2026) Basilisk install, with a correctly resolved (non-aliased) sampling
rate.** Whether this also holds on the actual period-correct 2025
Basilisk build is the open question job `4356316` will answer — do not
treat the "Basilisk-version-drift ruled out" claim below as settled until
that result is in.

**RESOLVED — job `4356316` (period-correct 2025-Basilisk, properly
fine-sampled) result:**

    ux_liq_rms/U_bio at t/T_p=[29,31]: peak = 9.4169 (vs. 9.4366 on
    2026-Basilisk -- agree to the 3rd significant figure). Full run
    (t/T_p up to 202): max 9.5502, n=6011 samples.

**Basilisk-version drift is definitively ruled out.** Kim et al.'s own
literal driver code, compiled against the actual Basilisk snapshot from
one day before they uploaded it, with properly resolved sampling, gives
essentially the identical large peak (~9.4) as it does on the current
2026 install. This is now a fully settled, three-way-confirmed number:
our own fork (5.34), upstream on 2026 Basilisk (9.44), upstream on
2025-era Basilisk (9.42) — all in the same regime, all ~6-12x above Kim
et al.'s own published Fig. 13a / Appendix A value of ~0.1-0.8, at the
exact same condition (theta=7deg, f_b=32.5rpm), using the exact metric
their own paper text specifies (`u_x,rms` on the liquid phase, in the
non-inertial frame, over `t/T_p=29-31`).

**Final answer to "why can't we reproduce Kim et al.'s results":** it is
not this project's modifications (ruled out repeatedly, term-by-term).
It is not sampling/aliasing (ruled out by fixing sampling on every run
compared). It is not the execution environment. It is not a Basilisk
version drift between 2025 and 2026 (ruled out directly, just now, by
building and running the actual period-correct compiler). **Kim et al.'s
own published driver code, run on their own era's toolchain, does not
reproduce their own published figure.** The remaining open possibilities
are outside what source-code archaeology can resolve: either the
published figure was generated from a different run/configuration than
what's in the public `DriverCodes` repository, or there is a
misunderstanding of the figure's actual normalization/axis convention
that the paper text does not fully disambiguate (e.g. `U_b` might be a
measured/fitted quantity in their actual analysis pipeline rather than
the analytical formula stated in the text, even though that formula
checks out algebraically against everything else). Both are now the
leading candidates, in place of "something in our fork" or "something in
Basilisk's evolution" — both of which are now closed.

---

**Result (coarse sampling, `i_norm=1000` — see caveat below):** the clean,
4-change-only upstream build (`BioReactor_clean`, and the earlier
debug-instrumented cross-check) both give `ux_liq_rms/U_bio` values of
**5.6-9.2** across `t/T_p≈25-56`, using upstream's own `Ly=0.286` (not our
`0.284`) for `U_bio`. That is comparable to or HIGHER than our own fork's
control run (`46acc8f0`, 32.5rpm/theta7, peak 5.34). Kim's Appendix A
figure reports ~0.1-0.8 for this exact quantity. **The ~6x-and-up
discrepancy vs. the published figure is present in literal upstream code,
essentially unmodified, run under the current Basilisk install. It was
never something introduced by this project's changes.** This reframes the
entire investigation: the open question is not "what did our fork break"
but "why does even Kim et al.'s own driver, as published, not reproduce
Kim et al.'s own published figure under this Basilisk version" — i.e. the
Basilisk-version-difference hypothesis (always the fallback candidate,
never previously testable) is now the leading one, or there's still an
error in how this test itself is set up (see caveats).

**Real, important caveat on the number above:** upstream's own
`i_norm=1000` samples statistics only once every ~2 non-dimensional time
units — under 1/3 of a rocking period (`T_per_nd≈0.607`). That is far too
coarse to resolve a smoothly oscillating quantity without severe
phase-aliasing — the exact problem this project already found and fixed
for its own `t_out` (commit `1c3440c`). The five points recorded
(9.2, 5.6, 2.1, 8.4, ...) jump around rather than tracing a smooth curve,
consistent with quasi-random phase sampling, not a reliable peak
extraction. **Submitted a properly-sampled rerun** (`i_norm` 1000→10,
documented in-place as a sampling-only change, no equation/field/timestep
touched — same justification already accepted for the `t_end` truncation)
via real `sbatch` (job `4355068`, dedicated 4-CPU allocation) — result
pending, see next entry.

**Also discovered, independent of the above, a real process-hygiene
mistake this session:** ran multiple compute-intensive test builds as
backgrounded (`&`/`nohup`) processes directly in the interactive coding
shell, rather than through `sbatch`. Checked `nproc` mid-session: this
shell has exactly **1 CPU**, on a node with `load average: 34.5` from
*other users'* unrelated jobs (Gaussian, Python) — a heavily oversubscribed
shared allocation, not a dedicated compute reservation. Also discovered
independently: upstream's own `event normcal(i+=i_norm)` has no
`t<=t_end` bound (same "runs forever" pattern this project already found
and fixed for `acceleration`/`dump_checkpoint` in its own fork, per
`hypothesis_ledger.json`) — so those background runs would never have
stopped on their own; killed both manually after collecting enough data
across the target window. Corrected going forward: the fine-sampling
rerun above was submitted via real `sbatch` with an explicit dedicated
allocation instead. This resource mistake does not affect the validity of
the coarse-sampled numbers themselves (CPU count doesn't change computed
physics), only how they were computed.

---

## 2026-07-28 (session 3) — Upstream crash SOLVED (own test-harness bugs, not
Basilisk/upstream), gdb worked via self-installed signal handler, then
rebuilt a minimally-patched clean version per explicit instruction

**gdb, corrected:** earlier claim that gdb was "unavailable" was about live
`ptrace`-based attach specifically. Confirmed via `/proc/self/status`
(`CapEff=0000000000000000`, `ptrace_scope=2`) that this is real and applies
uniformly — verified the user's own shell shows the same `CapEff=0`, so it's
a SLURM-job-step-wide policy (cgroup `job_4278934/step_interactive`), not
specific to this coding session. BUT: a self-installed `SIGSEGV`/`SIGFPE`/
`SIGABRT` handler using glibc's `backtrace()`/`backtrace_symbols_fd()` needs
no ptrace at all (runs inside the crashing process itself). Built that,
compiled with `-rdynamic -g`, and got a real, symbol-resolved backtrace on
the very first try. `addr2line` on the resolved addresses gave exact
file:line for every frame.

**Root cause of the persistent segfault, found via the backtrace + addr2line
(NOT a Basilisk-version incompatibility, NOT anything in either codebase's
real physics):**
1. First crash resolved to `event_do` (`grid/events.h:175`) calling into
   `out_files_initial` (a real user event in upstream `BioReactor.c`,
   `event out_files_initial(t=0; ...)`, not qcc-generated boilerplate as
   first assumed) → `fopen("Data_all/...", "wb")` → **`NULL`** (directory
   doesn't exist) → the next line's `fprintf` to a null `FILE*` segfaults.
   Upstream's own `BioReactor.sh` does `mkdir -p Data_all Data_specific
   Fig_vor Fig_vol Fig_tr Fig_oxy` before running — I never did, in any of
   this session's test harnesses. Purely a missing-directory bug in my own
   test setup.
2. After creating the directories, a second crash (`SIGFPE`, not `SIGSEGV`
   — needed to add that signal to the handler too) resolved to `vof_2`
   (`henry_oxy2.h:61`): `double a = c[]/(f[]*c.alpha + (1.-f[]));`. Root
   cause: my earlier "rule out tracer/oxygen" test (`TRACER=0`, `OXYGEN=0`)
   was **invalid**. Upstream's `stracers = {c,oxy,c1,c2,c3}` and this vof
   loop are unconditional — no `#if TRACER` guard — but `c.alpha` is only
   ASSIGNED inside an `#if TRACER {...} #endif` block in `main()`. With
   `TRACER=0`, that assignment is compiled out, `c.alpha` stays at its
   zero default, and the denominator becomes exactly `1-f` — zero in every
   pure-liquid cell — `0/0` → `SIGFPE` under Basilisk's FE-trap. That
   "refuted, ruled out tracer/oxygen" conclusion from earlier this session
   was never actually valid; it was testing a self-inflicted div-by-zero,
   not "no tracers."
3. With both fixed AND upstream's real defaults restored
   (`EMBED=1,OXYGEN=1,TRACER=1`), **upstream's actual code runs cleanly** —
   confirmed past 6000+ iterations, `t>12` (of a truncated `t_end=24`),
   zero crashes, in a build that ALSO still had the `q.embed_flux=NULL` fix
   applied (necessary independent of the above two -- see next entry --
   uninitialized struct field, real bug, unrelated to directories/flags).

**Per explicit instruction ("minimally modified... don't poison the well"):**
rebuilt from a fresh copy of the real upstream `BioReactor.c`/`henry_oxy2.h`
with ONLY four changes, each documented in-place and verified by `diff`
against the untouched upstream files to contain nothing else:
- `L0 = 1.[0]; DT = HUGE[0];` (was `L0 = LL;`) — dimensional-annotation
  compat patch only (project commit `8d6ae01`), value unchanged (`LL`=1.0,
  `HUGE` = upstream's own implicit unbounded default).
- `henry_oxy2.h` `set_prolongation`/`set_restriction` API rename (same
  commit) — API surface only, `.dirty` was removed from Basilisk's
  `_Attributes` in the version this project compiles against.
- `q.embed_flux = NULL;` — the one substantive fix, necessary for the code
  to run at all (see above), not a physics change (initializes a struct
  field to the value its very next use already assumes).
- `t_end` truncated 250.0→24.0 — does not touch the simulated equations,
  only how long the (identical, deterministic) run continues past the
  `t/T_p=29-31` window this test needs.
No debug prints, no signal handlers, no local Basilisk-header overrides in
this version — those were legitimate for crash bisection but have no place
in the file actually used for the reported comparison number. Compiled
against the real, unmodified Basilisk headers (no `-I.` shadowing).
Running now (`/oscar/scratch/eaguerov/tmp/kim_upstream_clean/run_test/`,
`BioReactor_clean`, pid `3971111`) alongside the earlier debug-instrumented
build (pid `3946838`, further along, kept only as a same-physics
cross-check since debug fprintf/signal-handler code cannot affect any
computed field). **Velocity comparison number pending — see next entry
once both finish.**

---

## 2026-07-28 (session 2) — Bisecting the velocity mismatch by literal reversion, and a blocked attempt at running raw upstream

**Method shift, per direct feedback:** rather than reasoning about whether
upstream and our code "should" be equivalent, revert one piece of OUR
working code to upstream's literal formula at a time, rerun the same cheap
control condition, and read off the number. All tests below use fidelity 6,
theta=7deg, f_b=32.5rpm, cold start (control run: `46acc8f0`,
`ux_liq_rms/U_bio` peak over the last 2 periods = **5.34**, vs. Kim's
Appendix A range of ~0.1-0.8 — i.e. our own working code is already ~6-7x
too high before touching anything).

**Attempted: run Kim's actual raw upstream `BioReactor.c`/`henry_oxy2.h`
under our current Basilisk, minimally patched only for the two documented
compatibility fixes (commit `8d6ae01`).** Segfaults almost immediately
(within the first i++ event group, before any user timestep completes).
Bisected via fprintf instrumentation (ptrace/gdb unavailable in this
sandbox; ASan unavailable too, `libasan.so.6` missing) and local-copy
`-I.` header overrides:
- Found and fixed ONE real, independent bug: `henry_oxy2.h`'s
  `tracer_diffusion` event declares `struct HDiffusion q; q.D=D; q.beta=beta;`
  with `q.embed_flux` left as **uninitialized stack garbage**, then later
  reads `if (!q.embed_flux && ...)`. This is exactly the H10/H11 bug our own
  hypothesis_ledger.json already found and fixed (`src/henry_oxy2.h:347`,
  `q.embed_flux = NULL`) — but that fix was framed as restart-specific.
  It is NOT restart-specific: `tracer_diffusion` runs on every timestep,
  fresh start included, so upstream's raw code needs this fix just to not
  crash on a cold start, full stop.
- Applying that fix was not sufficient — still segfaults, at the same
  point, even with `TRACER=0`, `OXYGEN=0`, `EMBED=0`, and without
  `-fopenmp` (ruled out tracer/oxygen transport, embedding, and threading
  entirely as the cause). Crash is somewhere in/around `vof(i++)`
  (confirmed via `-I.` local-copy instrumentation of `vof.h` reaching the
  binary, verified by `strings` on the compiled executable) but the
  instrumented fprintf as the literal first line of that event body never
  printed — crash is not inside the function body itself, more likely in
  Basilisk's own event dispatch/scheduling around it.
- **Closed as inconclusive.** Could not pin down further without a
  stack-trace tool. This blocks "run raw upstream directly" as a way to
  bisect the velocity mismatch — pivoted to reverting pieces of our own
  code instead (see below), which doesn't need upstream to run at all.

**Reverted our own code to upstream's literal formulas, one piece at a
time, same control condition:**

1. *Ramp shape + duration*: upstream's exact linear-over-30-physical-
   -seconds ramp (`ramp_dur = 30./T_bio; alpha = x_ss` instead of
   smooth-step) instead of our smooth-step-over-3-cycles. Result:
   `ux_rms/U_b` peak = **5.03**, a ~6% change *relative to our own control
   (5.34)* — NOT 6% of the way toward Kim's ~0.8 target. Still ~6x too
   high vs. Kim either way. **Refuted** — ramp has nothing to do with it.
2. *Multi-harmonic forcing loop structure*: upstream's literal unrolled
   single-harmonic formula (`Th_max2=alpha*Th_max_deg; Th=Th_max2*sin(...)`,
   no loop, no phase interpolation) instead of our generalized
   `for (k=1..n_harmonics)` sum (which reduces to the same formula at
   n_harmonics=1, but tested the actual literal old code path, not just
   the algebra). Result (stacked on top of test 1's ramp reversion):
   peak = **4.74**, an ~11% change *relative to the control* — again not
   11% closer to Kim's target. Still ~6x too high vs. Kim's ~0.8.
   **Refuted.** Neither reversion closed any of the gap to Kim; both are
   noise-level perturbations around the same ~5x-6x-too-high baseline.

**Also checked: `normf()`'s actual RMS definition** (read the real
source, `basilisk/src/utils.h:138-153`, rather than assuming) —
volume-weighted RMS, `sqrt(∫f²dV/∫dV)`, denominator is the *whole tank*
(water+air, wherever `cm>0`), not water-only. Main.tex's Appendix A text
explicitly says "liquid-phase" — so Kim's actual reported curve likely
normalizes by liquid volume only (half the tank at fill_level=0.5), which
would make the *true* liquid-only RMS **larger** by `√2` than what we
compute — i.e. correcting this would make our mismatch worse, not better.
Confirmed not the explanation; not implementing it.

**gdb retried on request (2026-07-28, later same day), confirmed genuinely
blocked, not just untried:** `module avail gdb` — no such module (only
`gdbm`, an unrelated library); `/usr/bin/gdb` exists already. Confirmed
this bash session runs inside an active SLURM allocation (compute node
node2333, `SLURM_JOB_UID` etc. set), not the login node, so it wasn't a
login-node restriction. `gdb -batch -ex run -ex bt --args
BioReactor_upstream_i 0.25 7 32.5` still fails: `ptrace: Operation not
permitted`. Root cause identified precisely this time:
`/proc/sys/kernel/yama/ptrace_scope` = `2` ("admin-only" — no ptrace
without `CAP_SYS_PTRACE`, not even parent-process-spawns-and-runs-child,
which scope `1` would normally allow). Also tried `coredumpctl list` to
read an already-generated core file post-mortem (doesn't need live
ptrace) — blocked too: "No journal files were opened due to insufficient
permissions." This is a kernel/container capability restriction on this
specific sandbox, not a missing tool — no module load or retry fixes it
without a genuinely different execution environment (e.g. a `salloc`
session with different container privileges, if that's even available
here).

**Confirmed with direct evidence, not just the error message, on a second
retry the same day:** `gdb --version` works fine (16.3-2.0.1.el9, real
binary, not missing) — so the earlier failure was never about gdb itself.
`/proc/self/status` shows `CapEff: 0000000000000000` — this shell has
*zero* effective Linux capabilities, including `CAP_SYS_PTRACE`, despite
the bounding set (`CapBnd`) nominally allowing it. Combined with
`ptrace_scope=2`, this is conclusive: no ptrace-based tool (gdb, strace,
core-file attach) can work here regardless of version or invocation
method. This is a sandbox/container privilege-drop, not a tooling gap.

**Where this leaves things:** every individual mechanism reverted or
checked today (ramp, harmonic-loop structure, RMS/volume definition) came
back negative. Combined with the previous session's findings (geometry/
embedding, forcing amplitude, physical constants, dimensionless numbers,
boundary conditions, pivot location, pseudo-force terms — all verified
identical to upstream), the velocity mismatch has survived every specific,
testable hypothesis so far. Genuinely open. Remaining untested angle:
Basilisk version itself (the actual numerical scheme in `vof.h`/
`centered.h`/`embed.h`), which can't be isolated without either (a) an old
Basilisk install to compile upstream against, or (b) finishing the crash
bisection above with better tooling than this sandbox allows.

---

## 2026-07-28 — Reproducing Kim et al.: root-cause hunt for the tau/velocity mismatch

**Context.** Neither our L9 nor L10 sweep reproduces Kim et al. (2024)'s
reported `tau_100_max`/`tau_mean_max` (see
`docs_site/explanation/kim-et-al-validation.md`). Working case throughout:
θ=7°, f_b=32.5 RPM (`omega_b=3.403392`) — Kim's own baseline/representative
condition, richest in comparison data (Fig. `simul_setup`, `tau_Ediss_evol`,
Appendix A convergence all use it).

**Hypotheses tested and REFUTED (all using existing data, zero/near-zero
new compute):**

1. *Window-length mismatch* (Kim's literal "one period" vs. our
   multi-period QSS window). Reprocessed existing `shear_stress.dat` with
   a true one-period window → fit got worse, not better. Refuted.
2. *Ramp-skip duration* (Kim's fixed 30 physical seconds vs. our fixed 3
   rocking cycles). Reprocessed existing data with Kim's 30s skip → no
   improvement, sometimes worse. Refuted. (Note: this only tested a
   *post-hoc analysis window* shift, not re-running with the gentler
   30s forcing ramp itself — see open threads below.)
3. *Checkpointing/chaining as driver of the L9-vs-L10 sign flip.* Isolation
   experiment (`experiments/l9_l10_checkpoint_isolation_test_30rpm/`):
   single-shot vs. 3-segment chain, same fidelity 9, 30 RPM. `tau_100_max`
   moved −3.0%, `tau_mean_max` −2.6% — too small to explain 20-50%+ gaps.
   Refuted as primary cause. (First attempt at this used a corrupted
   17.5 RPM baseline that was later discovered and retracted — see
   `experiments/l9_l10_checkpoint_isolation_test/` manifest for that
   post-mortem.)
4. *Geometry/embedding difference vs. Kim's actual upstream code.*
   Fetched the real source (`github.com/rcsc-group/BioReactor/DriverCodes`)
   after the user linked it. Kim's x-walls are plain grid-aligned domain
   edges (`u.n[left/right]=dirichlet(0)`), embedding used only for the
   y-direction (top/bottom plates). Our fork's geometry formula looks like
   it embeds both x and y (`pow(|x/a|,n)+pow(|y/b|,n)`), which raised a
   real concern. Ran `experiments/geometry_shape_test_32p5rpm/` (fidelity
   6, cold start, `n=8` vs `n=60`) to test corner-rounding sensitivity —
   got **bit-identical results.md** for both. Root cause: `n>=8` hits a
   special-cased exact-rectangle branch (`intersection(a-|x|, b-|y|)`),
   bypassing `pow()` entirely — the test compared a condition against
   itself. Findings written to that dir's `_findings.md`.
   Then checked properly (no new compute): `L0=1.0` (domain box spans
   `x,y∈[-0.5,0.5]`), and `a_nd = geometry.a/L_bio = 1` always by
   construction. Since the box half-width (0.5) < `a_nd` (1), the x-embed
   constraint is **never actually binding** — x-walls are, in practice,
   already just the plain box edges, matching Kim exactly. Confirmed BCs
   are line-for-line identical to upstream too. **Fully refuted** — no
   embedding-vs-plain-BC difference exists here after all.

**Verified as matching Kim's upstream exactly (source-diffed, not
assumed):**
- Geometry dims (`a=0.25m`, aspect ratio 0.284), fluid properties
  (`rho_w=1e3`, `mu_w=1e-3`, `rho_a=1.225`, `mu_a=1.81e-5`), gravity
  (9.8), surface tension (0.0728) — `src/BioReactor.c:91-96` vs.
  upstream `BioReactor.c`.
- Derived dimensionless numbers for this condition: `Re_w≈20585` (Kim's
  range 9.5k–24k ✓), `We_w≈23.3` (range 5–31 ✓), `Bo_w≈8410` (Kim states
  8400 ✓, geometry-only so RPM-independent).
- Forcing motion `θ(t)=θ_max·sin(ω_b·t)`, single harmonic, no doubling.
- Pseudo-force terms (gravity/Coriolis/centrifugal/Euler) in
  `event acceleration` — byte-for-byte identical to upstream, including
  `L_piv=0.143`.
- `U_bio`/`U_b` non-dimensionalization formula — algebraically identical
  to Kim's stated `U_b=L_x(1+tanθ/2β_b)/2T_p`.

**Real, still-live findings (not yet explained):**

- **`tau_100_max` does not converge with mesh resolution.** Using existing
  data (fidelity 7, run `82ee427c`, free) plus a purpose-built short-window
  test (fidelity 9 `8b3c0065`, fidelity 10 `e8ebf9f5`,
  `experiments/l9_l10_short_window_test_30rpm/`), all reprocessed over the
  identical window `t=[6.0,8.5]`, 30 RPM: f7→0.0417, f9→0.1367 (+228%),
  f10→0.2902 (+112%). No plateau. Crosses straight through Kim's value
  (0.1735) rather than approaching it. `tau_mean_max` is flat over the
  same range (~5% drift, f7→f10) while still sitting 35-65% below Kim's
  value throughout — so THAT metric's gap is provably not a resolution
  problem.
- **Shear-stress stencil is missing embedded-boundary metric factors.**
  `src/BioReactor.c:729-731` computes `du_dy`/`dv_dx` via plain
  `(f[0,1]-f[0,-1])/(2*Delta)`, and the comment claims this "mirrors
  vorticity() in basilisk/src/utils.h" — but the real `vorticity()`
  (`basilisk/src/utils.h:286-292`) weights by face metric factors
  `fm.x`/`fm.y`/`cm` and divides by `2*(cm[]+SEPS)*Delta`, needed for
  correctness near embedded cut-cells (i.e. near walls — exactly where
  peak shear stress lives). The comment mirrors the *indexing*, not the
  actual metric-corrected formula. Real bug, not yet fixed (holding per
  "only strictly numerically necessary changes" — see below).
- **Mean velocity is ALSO wrong, not just shear stress.** Warm-started a
  cheap continuation from run `30fb2321`'s checkpoint (`t=18.85`,
  `t/T_p=31`, already deep in Kim's own quasi-steady comparison window;
  run `8013ca72`, job `4344492`, fine sampling `t_out=0.02`). Computed
  `ux_liq_rms/U_bio` over `t/T_p=35-37` (clean, alias-free — verified by
  re-checking the same quantity on the *old* coarsely-sampled run first,
  where adjacent points swung wildly, e.g. 0.29↔0.03, confirming that WAS
  aliasing before trusting the new fine-sampled series): smooth periodic
  curve, peaks at **~4.2**, troughs at **~0.4**. Kim's Appendix A figure
  (`docs/kimetal2024/Figures/Fig_append1.pdf`, verified by rendering the
  actual PDF, not the caption text) shows `u_x,rms/U_b` oscillating
  ~0.1–0.8. **Our peak is >5x theirs.** This is a primitive-field
  quantity, no derivatives involved — rules out "shear-stress-specific
  numerical quirk" as the root cause, since velocity itself is this far
  off and the metric-factor bug above can't touch a non-derivative field.
- **A prior session's git commit made a false, unverified claim.** Commit
  `1c3440c` (Jul 14) claims an A/B test moved `tau_100_max` from 0.37x to
  1.03x Kim's value by fixing `t_out` (0.1→0.02 sampling), citing
  `experiments/hypothesis_ledger.json` — that file was never actually
  updated with the entry. Found the likely actual run pair
  (`/oscar/scratch/eaguerov/mpi_runs/8994c04a.STALE_pre_tout_fix` and
  `_tout_test_22p5`, 22.5 RPM, fidelity 9) and recomputed both directly:
  `tau_100_max` 0.08538 vs. 0.08542 (0.371x vs. 0.372x Kim) —
  **essentially no difference.** The claim does not survive checking
  against its own cited evidence. `t_out=0.02` itself is probably
  harmless (finer sampling, no reason it would hurt), but its stated
  justification is false and should not be trusted as "this was already
  fixed."

**Open threads / next falsifiable steps:**
- Ramp PROFILE (not just post-hoc analysis window): Kim's actual forcing
  ramp is a genuine 30-second linear ramp (`Th_max2=(Th_max/t_change_st)*t`
  in upstream `BioReactor.c`); ours is a 3-cycle *smooth-step* ramp
  (`src/BioReactor.c` commit `7ec98f9`, then `8ab1d1e` changed the ramp
  shape itself from linear to smooth-step for checkpoint-restart reasons).
  Never tested whether the actual forcing profile during the transient
  affects the eventual quasi-steady amplitude (should not, physically, if
  the system truly reaches the same periodic attractor — but "should not"
  isn't evidence).
  Basilisk version difference (Kim's upstream predates "basilisk 2026",
  ours is a current install) — untested, would need an old Basilisk build
  to check, not something fixable in our driver even if true.
- Real fix candidate identified but NOT applied: the `vorticity()` metric
  factor mismatch in the tau stencil. User's explicit instruction: only
  make changes that are strictly numerically necessary — no speculative
  fixes. This one has real evidence (diffed against Basilisk's own
  canonical function) but has NOT been shown to explain the *velocity*
  mismatch (a non-derivative quantity), so it is at most a partial
  explanation for the shear-stress-specific portion of the gap. Not yet
  applied pending further diagnosis.
- Still no explanation for the 5x mean-velocity mismatch. Every
  parameter, dimensionless number, and force-term check against the real
  upstream source has come back matching. This is the main open mystery.

**Housekeeping done alongside:** annotated `src/BioReactor.c` and headers
with inline markers distinguishing project additions/deletions from Kim's
original code (see commit for this diary entry). L9 video for this case
rendered and sent (`runs/8013ca72/volume_fraction*.mp4`).
