# Experiment diary

A lab notebook for numerical experiments on this project. Entries are
written as the work happens, not reconstructed afterward. Each entry
should let someone else (or future-us) reproduce the run and understand
why it was done, what it found, and what it does and doesn't prove.

Convention: newest entries at the top. Link run_ids / job_ids / commit
hashes exactly, not "the run from earlier."

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
