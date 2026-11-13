// init.h -- event init: geometry, initial fields, or restore from a checkpoint.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

// ================================================================== //
//                      INITIAL CONDITIONS                            //
// ================================================================== //
event init (t = 0)
{
  if (restart_file) {
    // ── Checkpoint restart ────────────────────────────────────────────────
    // restore() is called here so that Basilisk's grid infrastructure is
    // fully active.  It sets t = t_checkpoint and restores all fields.
    // All timing was already computed in main() from params.t_checkpoint so
    // the event system has the correct t_dump_checkpoint before run().
    //
    // [PROJECT FIXED, 2026-09-07] restore() internally calls Basilisk's own
    // dump_list(all, true) to decide which scalars to populate, and that
    // function skips any scalar whose CURRENT (in-memory, on THIS restoring
    // process) .nodump flag is true. p/pf default to nodump=true (per the
    // comment below) in a fresh process, regardless of what the DUMPING
    // side set before calling dump(). Without this override, restore()
    // finds "p" in the file's field-name table, fails to match it against
    // its own (p-excluded) scalar list, and silently routes that column's
    // data to a discarded placeholder instead of the real p[] field --
    // meaning p (and pf) were NEVER actually restored, staying at their
    // fresh-declaration default (exactly zero), for every restart this
    // project has ever run, no exceptions. Confirmed directly (Test C,
    // diary.md 2026-09-07): u.x/u.y/f matched the pre-dump state to full
    // double precision, but p came back as EXACTLY zero (min=max=sum=
    // stddev=0) regardless of whether the restart was a trivial identity
    // restart or an already-anomalous one -- this bug is present on EVERY
    // restart, independent of theta_max/omega_b changing at all, and is
    // the leading root-cause candidate for the whole cold-start-beats-
    // warm-start investigation (a from-scratch pressure solve is forced on
    // every restart instead of the "tiny correction" the code intended,
    // directly contradicting the p.nodump=false comment at the dump call
    // site below).
#if !REPLAY_FLUX
    p.nodump = pf.nodump = false;
#endif
    if (!restore (file = restart_file)) {
      fprintf (stderr, "ERROR: restore() failed to open '%s' — aborting\n", restart_file);
      exit (1);
    }
#if !REPLAY_FLUX
    p.nodump = pf.nodump = true;
#endif
    write_restart_diagnostic ("post_restore");
#if REPLAY_FLUX && EMBED
    scalar cs_dumped[];             // the cleaned geometry the coupled run had; compared below
    foreach()
      cs_dumped[] = cs[];
    // Coarse-level metric. A fresh coupled run never restricts fs after solid(): embed.h's event metric
    // restricts fs = 1 once at i = 0, and the multigrid of the scalar diffusion (h_relax, embed_flux) reads
    // the coarse levels of fs ever after. The dump has no face fields, and the restart path's own
    // restriction({cs,fs}) after solid() would give the coarse levels a different value than the run being
    // replayed had; that changes the oxygen solve (measured: 1e-3 after one step). So reproduce the fresh
    // run: fs = 1 on all levels here, and no restriction after solid() below. The coarse cs is dumped, and
    // restriction({D, beta, cm}) in tracer_diffusion refreshes it every step as in the coupled run.
    foreach_face()
      fs.x[] = 1.;
    restriction ({fs});
#endif


    // Re-apply the prolongation/restriction setup from event defaults(i=0) in
    // henry_oxy2.h.  That event fires at i=0 on a fresh start but is skipped on
    // restart (i resumes from the checkpoint value, so i==0 is never seen again).
    // Basilisk's dump/restore does NOT serialize scalar function-pointer attributes
    // (restriction, prolongation, refine), so they revert to defaults.  Without
    // restriction_volume_average the multigrid h_relax propagates NaN from solid
    // cells (cm=0 → d=0 → c[]=n/0) into fluid cells → NaN kLa and tracer sums.
    //
    // [PROJECT MOVED, 2026-09-12] This whole block (prolongation/restriction
    // re-application for stracers, u/p/pf/g/uf, and f) used to sit AFTER the
    // cross-level warm-start's refine() call below -- meaning refine()
    // interpolated every field using Basilisk's generic, non-embed-aware
    // default prolongation instead of the embed-aware one this project's own
    // (pre-existing) comments already say is required near the embedded
    // boundary. Moved here, before refine(), so cross-level interpolation
    // uses the correct prolongation from the start. Ordering is harmless for
    // an ordinary (non-cross-level) restart, since no refine() call happens
    // there either way. Falsifiable test: does this close (or shrink) the
    // ~32%/~49% L6->L7/L7->L8 cross-level tau discrepancy (diary.md
    // 2026-09-12)? Not yet re-run.
#if TREE
    for (scalar s in stracers) {
      s.refine = refine_embed_linear;
      set_prolongation (s, refine_embed_linear);
      set_restriction (s, restriction_volume_average);
    }
#endif

    // Re-apply centered.h's EMBED defaults (event defaults(i=0) is skipped on
    // restart because i resumes from the checkpoint value).  Without these,
    // p/pf/u/g use non-EMBED prolongation/restriction and uf uses
    // refine_face_solenoidal instead of refine_face.  The pressure gradient at
    // the embedded wall falls back to the generic gradient (no
    // pressure_embed_gradient), so the Poisson solve near the solid accumulates
    // small errors each timestep.  After ~300 steps (~1 T) the velocity
    // diverges → SIGFPE in tracer_diffusion/h_residual.
    // Confirmed experimentally: fresh run to t=15 is stable; restart crashes at
    // t=t_checkpoint+1.14T with both 4 and 16 MPI ranks.
#if TREE && EMBED
#if !REPLAY_FLUX
    uf.x.refine = refine_face;
    foreach_dimension()
      uf.x.prolongation = refine_embed_face_x;
    for (scalar s in {p, pf, u, g}) {
      s.refine = refine_embed_linear;
      set_prolongation (s, refine_embed_linear);
      set_restriction (s, restriction_embed_linear);
    }
    for (scalar s in {p, pf})
      s.embed_gradient = pressure_embed_gradient;
#endif // !REPLAY_FLUX

    // [PROJECT ADDED, 2026-09-13] cs/fs themselves have exactly the same
    // restart-attribute-loss problem as u/p/pf/g/uf above, and it was never
    // fixed: embed.h's own `event metric (i=0)` is what normally assigns
    // cs.refine=embed_fraction_refine, cs.prolongation=fraction_refine,
    // fs.x.prolongation=embed_face_fraction_refine_x, and it is the ONLY
    // place in Basilisk that calls `restriction ({cs,fs})` -- both skipped
    // on restart for the same reason (i != 0). Root cause of the
    // cross-level checkerboard (diary.md 2026-09-13): cs/fs are restored
    // (cs) or recomputed (fs, via the solid() call below) at LEAF level
    // only; nothing ever pushed the corrected leaf values up to the tree's
    // coarse (non-leaf) cached cs/fs. Confirmed against Basilisk's own
    // canonical idiom for recomputing solid() outside the normal i=0 init
    // path (src/test/neumann3D.c:43-47): solid() must always be followed by
    // reasserting cs/fs's refine+prolongation attributes and an explicit
    // restriction({cs,fs}) call. None of the six fix attempts already
    // documented below touched cs/fs's own attributes or tree consistency
    // -- they all operated on u/p/g instead.
    cs.refine = embed_fraction_refine;
    cs.prolongation = fraction_refine;
    foreach_dimension()
      fs.x.prolongation = embed_face_fraction_refine_x;
#endif // TREE && EMBED
    // Re-apply vof.h defaults: fraction_refine is set for f in vof.h's
    // event defaults(i=0), skipped on restart.  Without it AMR uses bilinear
    // prolongation for f, creating non-physical VOF fractions near the interface.
#if TREE
    f.refine = fraction_refine;
    set_prolongation (f, fraction_refine);
#endif // TREE

    // [PROJECT ADDED, 2026-09-10] Cross-level warm-start pilot (diary.md
    // 2026-09-10): restore() reconstructs the checkpoint's tree EXACTLY as
    // dumped -- a checkpoint written at fidelity F restores as a uniform
    // tree of depth F, even when THIS run's own NN/init_grid was set up
    // for a finer fidelity. Gated behind -DCROSS_LEVEL_WARMSTART=1, off by
    // default and NOT part of the normal restart path.
    //
    // [CORRECTED, 2026-09-10] First attempt used adapt_wavelet() with a
    // deliberately tiny (1e-30) tolerance on every field, reasoning that
    // any nonzero wavelet detail coefficient would trigger refinement.
    // Wrong, caught by the verification below on a smoke test (not
    // assumed correct and left unverified): grid->n came back 752 against
    // a wanted 16384. adapt_wavelet refines on the wavelet DETAIL
    // coefficient, which is bit-exactly 0.0 over large flat sub-regions of
    // these fields (f exactly 0 or 1 away from the interface, quiescent
    // near-zero velocity in the air phase) -- and 0.0 is not greater than
    // any positive tolerance, however small, so those regions never
    // refined at all. An error-driven criterion cannot deliver unconditional
    // uniform refinement no matter how the tolerance is tuned.
    //
    // Fixed: Basilisk's own `refine (bool cond)` macro
    // (grid/tree-common.h) is genuinely unconditional -- it loops
    // (foreach_leaf, refine_cell) until no leaf satisfies `cond` anywhere,
    // with no error estimate involved, and refines with field list `all`
    // (interpolates every declared scalar via its own .prolongation
    // attribute, not just the ones named explicitly). `level < params.
    // fidelity` refines every leaf below the target depth, repeatedly,
    // until the whole domain is uniform at that depth -- the actual
    // "spatially interpolate the coarse field onto the finer grid
    // everywhere" operation the experiment needs.
    //
    // Verified after the call, not assumed: depth() and the leaf count
    // must both match the fully-refined uniform grid this run's own NN
    // implies. If either check fails the run aborts rather than silently
    // continuing on a partially-refined grid.
#if CROSS_LEVEL_WARMSTART
    {
      // [CORRECTED, 2026-09-10] The first two attempts checked
      // grid->n == NN*NN, assuming the checkpoint's grid is uniform across
      // the WHOLE bounding box -- wrong, and caught by looking at the
      // actual pre-refine count instead of assuming it: restore() hands
      // back only 512 active leaves for a nominally-uniform L6 grid whose
      // full bounding box would be 4096, because embed.h's own machinery
      // keeps cells entirely outside the thin bag geometry coarse (no
      // fluid physics there to resolve) -- true of every run this project
      // has ever done, at every fidelity, not a defect in this checkpoint.
      // The actual success criterion is that refine() cleanly multiplied
      // whatever active-cell count existed by exactly 4 (one uniform
      // level, one pass), which is what a leaf refining into 4 children
      // means in 2D -- not that the whole box reached maximum depth.

      // [PROJECT ADDED, 2026-09-14] Recompute cs/fs analytically on the
      // COARSE grid BEFORE refine(), so refine_embed_linear has correct fs
      // to work with when it prolongates u/p/g.
      //
      // Why this matters: refine() prolongates u/p/g with
      // refine_embed_linear (embed-tree.h:268), and EVERY interpolation
      // branch there is gated on the COARSE fs -- bilinear and triangular
      // both require `coarse(fs.x,i) && coarse(fs.y,0,j)`, diagonal
      // requires `coarse(fs.x,i) && coarse(fs.y,child.x,j)` (or the
      // transpose); if all fail it falls to "pathological cases", whose
      // two 1D-gradient corrections are ALSO fs-gated, leaving
      // `s[] = coarse(s)` (piecewise-constant injection). So whatever fs
      // holds at refine() time directly controls the quality of the
      // interpolation.
      //
      // What fs actually holds after restore(), MEASURED on the L6->L7
      // pilot rather than assumed (fs is a face field, and Basilisk's
      // dump_list() skips face fields outright -- output.h:1037,
      // `if (!s.face && ...)`):
      //   cs   restores EXACTLY          (sum 1171.46, = analytic to 5 figs)
      //   fs.x sum 594.99 vs 1189.76 analytic
      //   fs.y sum 594.32 vs 1216   analytic
      // and, decisively, among faces that should be FULLY OPEN (both
      // neighbours cs > 0.99, no geometry anywhere near) roughly HALF come
      // back as hard ZERO: 586 of 1170 in x, 576 of 1088 in y, already so
      // the instant restore() returns. Since the branch tests above are
      // truthiness tests on fs, that alternating half-zero pattern throws
      // ~50% of cells onto the injection branch while the interleaved rest
      // interpolate properly -- which IS the period-2 "checkerboard", and
      // it is imposed over the WHOLE domain, not just near the wall.
      //
      // With this call in place: 0 of 1224 open x-faces and 0 of 1152 open
      // y-faces are zero, and the real 60-cycle L6->L7 pilot lands at
      // +0.04% on tau and +0.01% on dissipation versus the L7 cold start
      // (vs +29.7% / +219.8% without it). See diary.md 2026-09-14.
      //
      // The post-refine solid() call further below is still required and
      // still runs; it makes cs/fs exactly analytic at the FINE resolution,
      // which prolongation alone only approximates.
      //
      // Toggle for A/B verification: -DCROSS_LEVEL_PREREFINE_SOLID=0
      // reproduces the old (broken) ordering.
#ifndef CROSS_LEVEL_PREREFINE_SOLID
#define CROSS_LEVEL_PREREFINE_SOLID 1
#endif
#if EMBED && CROSS_LEVEL_PREREFINE_SOLID
      {
        double a_nd = params.geometry_a / L_bio;
        double b_nd = params.geometry_b / L_bio;
        if (params.geometry_n >= 8.)
          solid (cs, fs, intersection(a_nd - fabs(x), b_nd - fabs(y)));
        else
          solid (cs, fs, 1. - pow(fabs(x/a_nd), params.geometry_n)
                            - pow(fabs(y/b_nd), params.geometry_n));
      }
#endif

      // fs must be usable HERE -- refine() below interpolates u/p/g through
      // refine_embed_linear, which reads it. This is the check that would
      // have caught the 2026-09-14 bug on day one.
      assert_fs_sane ("pre-refine");

      long n_before = grid->n;
      refine (level < params.fidelity);
      long n_after = grid->n;

      // [PROJECT GENERALIZED, 2026-09-15] This check used to demand
      // `n_after == 4 * n_before` -- a hardcoded SINGLE-level assumption.
      // refine() loops internally until no leaf satisfies the condition, so
      // a k-level jump multiplies the active leaf count by 4^k in 2D (an
      // L6->L8 jump measured 512 -> 8192 = 16x = 4^2 and was rejected as
      // "incomplete" even though it was correct, job 6402091). A k=0 jump
      // (restarting at the SAME fidelity, e.g. the 2nd and later segments
      // of a chained run) leaves the count unchanged and was likewise
      // rejected, which forced chains to swap binaries mid-run.
      //
      // The real invariant is: the grid reached the target depth, and every
      // refined leaf split cleanly, i.e. the ratio is an exact power of 4.
      long ratio = (n_before > 0) ? n_after / n_before : 0;
      bool clean_pow4 = (n_before > 0) && (n_after == ratio * n_before);
      while (clean_pow4 && ratio > 1) {
        if (ratio % 4) { clean_pow4 = false; break; }
        ratio /= 4;
      }
      clean_pow4 = clean_pow4 && (ratio == 1);
      if (pid() == 0)
        fprintf (ferr, "cross-level warmstart: refine n %ld -> %ld "
                 "(want an exact power of 4), depth %d\n",
                 n_before, n_after, depth());
      if (depth() != params.fidelity || !clean_pow4) {
        fprintf (stderr, "ERROR: cross-level refine incomplete -- "
                 "depth()=%d (want %d), n %ld -> %ld (want n_after/n_before "
                 "to be an exact power of 4: 1 for a same-level restart, 4 "
                 "per level refined)\n",
                 depth(), params.fidelity, n_before, n_after);
        exit (1);
      }
      if (pid() == 0)
        fprintf (ferr, "cross-level warmstart: refined to depth=%d, "
                 "n=%ld cells\n", depth(), (long) grid->n);
    }
#endif

#if REPLAY_FLUX
    // Scalars on a finer grid than the flow (F3): refine the restored grid to replay_scalar_level. Order as
    // in the cross-level path above (CLAUDE.md): cs/fs from solid() on the coarse grid FIRST, because
    // refine_embed_linear and embed_fraction_refine read fs; then refine(); the fine solid() and the cleanup
    // follow below. The fs = 1 coarse-level convention of a fresh run is set again on the new levels.
    if (params.replay_scalar_level > params.fidelity) {
      if (depth() != params.fidelity) {
        fprintf (stderr, "ERROR: replay: restored grid has depth %d, params.fidelity says %d\n",
                 depth(), params.fidelity);
        exit (1);
      }
#if EMBED
      {
        double a_nd = params.geometry_a / L_bio;
        double b_nd = params.geometry_b / L_bio;
        if (params.geometry_n >= 8.)
          solid (cs, fs, intersection(a_nd - fabs(x), b_nd - fabs(y)));
        else
          solid (cs, fs, 1. - pow(fabs(x/a_nd), params.geometry_n)
                            - pow(fabs(y/b_nd), params.geometry_n));
      }
      assert_fs_sane ("pre-refine");
#endif
      refine (level < params.replay_scalar_level);
      if (depth() != params.replay_scalar_level) {
        fprintf (stderr, "ERROR: replay: refine reached depth %d, wanted %d\n", depth(),
                 params.replay_scalar_level);
        exit (1);
      }
      if (pid() == 0)
        fprintf (ferr, "replay: scalar grid refined from level %d to %d, %ld cells\n", params.fidelity,
                 params.replay_scalar_level, (long) grid->n);
      if (params.replay_diag != 2) {
      foreach_face()
        fs.x[] = 1.;
      restriction ({fs});
      }
      if (params.replay_diag) {
        stats s1 = statsf (f), s2 = statsf (oxy), s3 = statsf (c2), s4 = statsf (c1), s5 = statsf (cs);
        if (pid() == 0)
        fprintf (ferr, "replay: after refine f [%g,%g] sum %.10g, oxy [%g,%g], c2 [%g,%g], c1 [%g,%g], cs [%g,%g]\n",
                 s1.min, s1.max, s1.sum, s2.min, s2.max, s3.min, s3.max, s4.min, s4.max, s5.min, s5.max);
      }
    }
#endif // REPLAY_FLUX

    // fs (embed face fractions) does not survive dump/restore intact, so it
    // must be recomputed from the same static geometry used at fresh-start
    // (rocking is body forces, not a moving solid). Without this the NS
    // solver misreads the walls and the velocity collapses on the first step.
    //
    // [CORRECTED, 2026-09-14] This comment used to say "After restore, fs=0
    // everywhere". That is FALSE and it cost real debugging time (it was
    // taken at face value and used to build a wrong root-cause story; see
    // diary.md 2026-09-14). MEASURED on the L6 restart, before any
    // recompute: fs is roughly half-correct in aggregate (fs.x sum 594.99
    // vs 1189.76 analytic; fs.y 594.32 vs 1216), and cs restores EXACTLY
    // (1171.46 both). The real defect is structural, not a blanket zero:
    // ~50% of fully-open fluid faces (both neighbours cs > 0.99, geometry
    // nowhere near) come back as hard ZERO -- 586/1170 in x, 576/1088 in y.
    // That alternating half-zero pattern is what drives refine_embed_linear
    // onto its fs-gated pathological branch for half the cells during a
    // cross-level refine(), which is the actual origin of the "checkerboard".
#if EMBED
    {
      double a_nd = params.geometry_a / L_bio;
      double b_nd = params.geometry_b / L_bio;
      if (params.geometry_n >= 8.)
        solid (cs, fs, intersection(a_nd - fabs(x), b_nd - fabs(y)));
      else
        solid (cs, fs, 1. - pow(fabs(x/a_nd), params.geometry_n)
                          - pow(fabs(y/b_nd), params.geometry_n));
      // [PROJECT ADDED, 2026-09-13] solid() only writes LEAF cells
      // (fractions.h's foreach()/foreach_face() default to leaves). On a
      // cross-level refine(), the tree's coarse (non-leaf, pre-refine)
      // cs/fs cache is never touched by solid() and is never restricted
      // anywhere else on the restart path -- it keeps the STALE, aliased
      // pre-refine values (diary.md 2026-09-13's root cause) forever,
      // feeding the embedded Poisson/diffusion operators at coarse tree
      // levels every single timestep. Canonical Basilisk idiom for calling
      // solid() outside the normal i=0 init path always follows it with
      // this restriction (src/test/neumann3D.c:43-47); our restart path
      // was missing it.
#if TREE && !REPLAY_FLUX      // a replay keeps the coarse levels it set up above, see there
      restriction ({cs, fs});
#endif
    }
    // fs must also be usable here, on the FINAL grid the solver will run on
    // -- this covers the ordinary (non-cross-level) restart path too, where
    // fs is likewise broken by restore() and only repaired by the solid()
    // call just above.
    assert_fs_sane ("post-solid");
#endif

#if REPLAY_FLUX
    // The coupled run repeats fractions_cleanup(cs, fs) every step (event adapt of centered.h, trees
    // only), so the cs/fs it advects f with are the cleaned ones, and the dump holds them; solid() above
    // has just rebuilt the raw ones. The cleanup is deterministic and reaches a fixed point, so one call
    // restores the cleaned geometry. The uf=0 on closed faces that event adapt also does is already in
    // the recorded uf.
#if TREE && EMBED
    fractions_cleanup (cs, fs);
    {
      int ndiff = 0;
      double dmax = 0.;
      foreach (reduction(+:ndiff) reduction(max:dmax)) {
        double d = fabs (cs[] - cs_dumped[]);
        if (d > 1e-12) ndiff++;
        if (d > dmax) dmax = d;
      }
      if (pid() == 0)
        fprintf (ferr, "replay: cs after solid()+fractions_cleanup differs from the dumped cs in %d cells "
                 "(max %.3g)\n", ndiff, dmax);
    }
#endif
    if (params.replay_diag == 2)
      restriction ({cs, fs});
    if (params.replay_scalar_level > params.fidelity) {
      // measured at level 7: tracer values interpolated into the cells of the wall that the fine geometry
      // calls solid (cs = 0) feed the cut cells through small face fractions and the oxygen solve diverges
      // in the first step. A fresh run has 0 there (the solver itself sets them to 0, h_relax: d <= 0).
      foreach()
        if (cs[] == 0.)
          for (scalar s in stracers)
            s[] = 0.;
      boundary (stracers);
      restriction (stracers);
    }
    replay_open ();
#endif

    // [PROJECT NOTE, 2026-09-13] A cross-level warm-start (CROSS_LEVEL_
    // WARMSTART) introduces a real, root-caused artifact: a checkerboard
    // pattern in u/p in a ~5-coarse-cell-wide band
    // next to the embedded wall, which accounts for the ENTIRE observed
    // 32-49% domain-mean tau discrepancy against a same-condition cold
    // start (measured directly by exclusion, diary.md 2026-09-13). Root
    // cause: Basilisk's refine_embed_linear (embed-tree.h:268) picks its
    // per-child interpolation formula from the COARSE parent's fs/cs
    // neighbor pattern, and the coarse source grid's own discretization of
    // the curved bag boundary has a genuine period-2 aliasing in fs.x near
    // the wall -- refine_embed_linear faithfully (and correctly, per its
    // own contract) propagates that into the fine grid. SEVEN independent
    // fixes have been tried and removed from this file after all failing on
    // the real 60-cycle L6->L7 pilot (five had zero effect, one made it
    // worse, one -- the cs/fs reapplication + restriction({cs,fs}) added
    // above, 2026-09-13, the same restart-attribute-loss bug already fixed
    // for u/p/pf/g/uf but never extended to cs/fs itself -- produced a
    // trajectory BIT-IDENTICAL to the unfixed run, max diff 7.9e-7 over the
    // whole 60-cycle run): moving this prolongation re-application before
    // refine() (kept -- see below, it's a real separate bug, just not the
    // cause of this artifact), an explicit post-refine incompressibility
    // projection, plain non-embed-aware interpolation, near-wall corrective
    // smoothing at two widths (one verified via direct measurement to
    // actually clean the field at t=0), propagating that correction to the
    // tree's coarse-level cache via restriction(), and the cs/fs restart-
    // attribute fix. The checkerboard regenerates within under 1% of one
    // rocking cycle even from a verified-clean start, and does NOT appear
    // in an ordinary (non-cross-level) restart into the same converged
    // state, nor in a cold start forced through an near-instant ramp -- so
    // it is specific to the cross-level refine() path, not to onset speed
    // or to handling a converged flow in general, and lives in the solver's
    // own per-timestep dynamics, not in the initial condition, and not in
    // anything related to cs/fs's own tree consistency either (that fix is
    // kept anyway -- it's a real, separate correctness fix matching
    // upstream's own contract, just not the cause of this artifact). See
    // diary.md 2026-09-13 for the full investigation (kept there, not
    // here, to keep this file uncluttered) before attempting further fixes.

#if !REPLAY_FLUX
    // Rescale stored velocity and pressure to the new segment's non-dim frame.
    // [PROJECT FIXED, 2026-09-07] The comment this replaces claimed "U_bio ∝
    // omega_b (fixed geometry, theta_max)" and rescaled by su =
    // omega_b_prev/omega_b alone -- but U_bio = V_bio/(H_bio/2)/T_per, and
    // V_bio = L_bio/4*(H_bio + 0.5*L_bio*tan(theta_max)) depends on
    // theta_max too. That assumption is simply wrong for exactly the case
    // this project's settling-time sweeps rely on: checkpointing across a
    // theta_max-only change (omega_b_prev == omega_b, so the old code
    // silently applied su=1, no rescale at all, regardless of how much
    // theta_max changed). Fix: compute the FULL U_bio ratio (both omega_b
    // and theta_max contributions), not just the omega_b part.
    //
    // HONEST STATUS: dimensionally required, empirically inert. This was
    // first committed (646b7b3) as a "confirmed bug fix" and as the likely
    // seed of the cold-start-beats-warm-start anomaly. Both claims were
    // wrong. The anomaly is the restart forcing-phase error fixed in
    // main(); A/B'ing this term against -DLEGACY_SU=1 with that fix in
    // place (scripts/ab_test_su_rescale.py, diary.md 2026-09-07 (4)) gives
    //     Delta_theta = 0.1 (term = 0.26% on u): 1.0001 vs 1.0005
    //     Delta_theta = 5   (term = 14%   on u): 0.9999 vs 1.0002
    // against a 0.43-0.48% cycle-to-cycle spread -- no measurable effect at
    // any Delta_theta_max tested, including the 14% case. It is kept
    // because u_nd = u_phys/U_bio and a segment that changes U_bio must
    // convert, not because any experiment supports it; the reference state
    // simply re-equilibrates within its basin. Do not cite it as evidence
    // for anything.
    if (params.omega_b_prev > 0.) {
      // H_bio is a LOCAL in main() (not a global), so it must be
      // recomputed here from the globals it derives from (L_bio, Ly),
      // which main() has already set by the time any event fires.
      double H_bio_here  = 2.*L_bio*Ly;
      double T_per_prev  = 2.*pi/params.omega_b_prev;
      double Th_max_prev = params.theta_max_prev[0] * pi / 180.;
      double V_bio_prev  = L_bio/4*(H_bio_here + 0.5*L_bio*tan(Th_max_prev));
      double U_bio_prev  = V_bio_prev/(H_bio_here*0.5)/T_per_prev;
      // [PROJECT ADDED, 2026-09-07] -DLEGACY_SU=1 restores the original
      // omega_b-only ratio, so the theta_max contribution to su can be
      // A/B'd against this fix rather than asserted. It needs A/B'ing: the
      // change was committed as a "confirmed bug fix" on dimensional
      // grounds alone, and its one verification run moved the metric from
      // 5.8% to 5.4%, i.e. not at all. For a theta 6.9->7 restart the term
      // it adds is 0.26% on velocity, and the kicknoise run showed a
      // deliberate 0.3% per-cell random velocity perturbation decaying to
      // nothing (every quantity back to 1.000 within 0.1% over 160 cycles),
      // so it cannot matter at small Delta_theta_max. At theta 2->7 it is
      // 14%, which is the case worth measuring.
#if LEGACY_SU
      double su = params.omega_b_prev / params.omega_b;
#else
      double su = U_bio_prev / U_bio;
#endif
      (void) U_bio_prev;
      foreach() {
        u.x[] *= su;
        u.y[] *= su;
        p[]   *= su * su;
        pf[]  *= su * su;   // half-step pressure (BCG tracer advection), same scale as p
      }
      // uf (face-centered velocity) must also be rescaled.  centered.h's
      // advection scheme (BCG), the VOF step, and the AMR adapt-at-t=0
      // all use uf directly before it is rebuilt from u.  Without this,
      // the first Poisson projection sees a factor-of-(1/su)^2 kinetic
      // energy mismatch and destroys the restored velocity field.
      foreach_face()
        uf.x[] *= su;
      boundary ({u.x, u.y, p, pf, uf.x, uf.y});
      // Scale g (combined pressure-gradient + acceleration term from BCG predictor).
      // g is stored in the checkpoint in seg0 non-dim units.  With ramp=1 immediately,
      // the acceleration part of g is correct (w_bio_st is constant across omega_b).
      // The pressure-gradient part scales as su² (same as p).  Scale by su² to keep
      // the pressure-gradient contribution accurate in the first BCG half-step.
      // [PROJECT ADDED, 2026-09-08] G_RESTART_MODE. The uniform su^2 below is
      // provably wrong: g combines the pressure gradient with the whole
      // acceleration field a, and those do NOT share a scaling law --
      //   -sin(Th)/Fr^2            gravity      ~ su^2
      //   2*Th_d*u                 Coriolis     ~ su   (Th_d is su-invariant)
      //   Th_d^2*(x+L_piv sin Th)  centrifugal  ~ 1
      //   Th_2d*(y+L_piv cos Th)   Euler        ~ 1
      // so scaling all of g by su^2 mis-scales three of the four terms, by
      // |su^2-1| on the su-invariant ones: 24% at theta 2->7 (su=0.874,
      // converges) but 33% at omega 37.5->32.5 (su=1.154, escapes) and 71%
      // at omega 17.5->32.5 (su=0.538, escapes). Harmless while su stayed
      // near 1 on theta-only restarts; exposed by omega restarts.
      //
      // g only feeds the first BCG predictor half-step and is rebuilt by
      // centered_gradient(p,g) at the end of it, and a COLD start runs that
      // first half-step with g==0 (centered.h's event init never sets it).
      // Mode 2 therefore reproduces cold-start behaviour exactly, mode 1
      // leaves the restored value untouched, mode 0 is the current su^2.
      // Measured on the worst case (omega 17.5->32.5, su=0.538, 71% error
      // on the su-invariant terms): all three modes give tau/tau_ref
      // 1.2273-1.2308, identical to within noise. g simply does not matter
      // -- it survives one predictor half-step. So the su^2 scaling was
      // never the cause of the omega escape (that is a basin crossing set
      // by the size of the jump), but it is still provably wrong code, and
      // the default is now mode 2: match a cold start exactly rather than
      // apply a scaling that is right for one of four terms.
#ifndef G_RESTART_MODE
#define G_RESTART_MODE 2
#endif
#if G_RESTART_MODE == 0
      foreach() {
        g.x[] *= su * su;
        g.y[] *= su * su;
      }
      boundary ({g.x, g.y});
#elif G_RESTART_MODE == 2
      foreach() {
        g.x[] = 0.;
        g.y[] = 0.;
      }
      boundary ({g.x, g.y});
#endif
    }

#endif // !REPLAY_FLUX

    // Reset ALL stracers to zero at EVERY multigrid level.  reset() zeroes
    // owned cells; boundary() communicates zeroed leaf ghost values.  But
    // coarse-level ghost cells (owned by a neighbouring MPI rank) are NOT
    // updated by boundary() — they retain stale checkpoint values.  When
    // h_relax runs at a coarse V-cycle level it reads a[1] which can be a
    // coarse ghost cell, producing a spurious non-zero restricted residual
    // (b[] >> 1) and driving a gradual multigrid divergence over ~5 periods.
    // restriction() propagates leaf=0 to coarse own cells and then runs
    // halo_restriction which communicates those zeroed coarse values to
    // neighbouring ranks' ghost slots, closing the ghost-cell gap.
    //
    // [PROJECT FIXED, 2026-09-15] The zeroing is now conditional. Note what
    // the diagnosis above actually says: the problem is STALE COARSE-LEVEL
    // values, and restriction() is what fixes it, by recomputing coarse cells
    // from the leaves and communicating them to neighbours' ghost slots. The
    // leaf values are not the problem -- and on a CONTINUATION they are the
    // entire experiment. Measured (probe job 6411170): c2 and oxy round-trip
    // through restore() bit-exact and were then destroyed here, which made
    // every chained mixing/kLa run silently measure from a blank tracer
    // field. restore() and solid() were both ruled out by direct measurement
    // before this line was touched. See diary.md 2026-09-15 (8).
    if (!params.restart_continue)
      reset (stracers, 0.);
    boundary (stracers);
    restriction (stracers);
  } else {
    // ── Fresh start ────────────────────────────────────────────────────────
    // Parametric bag geometry (dimensionless semi-axes)
    double a_nd = params.geometry_a / L_bio;
    double b_nd = params.geometry_b / L_bio;  // == Ly

    // [PROJECT CHANGED] Upstream: `double y_init = 0.0;` — fixed liquid
    // level (always half-full, matching Kim's own 0.5 fill fraction). Made
    // parametric for the sweep/optimization pipeline. At fill_level=0.5 this
    // reduces to y_fill = b_nd*(2*0.5-1) = 0, i.e. y_init=0.0 exactly,
    // matching upstream for Kim's own condition.
    // Fill level: liquid occupies fill_level fraction of bag height, measured from bottom
    double y_fill = b_nd * (2.*params.fill_level - 1.);
    fraction(f, y_fill - y);

    // [PROJECT CHANGED] Upstream's `init` event only ever builds a plain
    // rectangle here: `solid(cs, fs, intersection(-(y-0.5*Ly), -(-y-0.5*Ly)))`
    // — note this bounds ONLY y; the x-direction walls are upstream's plain
    // domain box edges (u.n[left]/u.n[right] boundary conditions below), not
    // embedded at all. This fork generalized the shape to a parametric
    // superellipse so bag corner-rounding could be a sweep variable. VERIFIED
    // (diary.md, 2026-07-28) that this generalization does NOT change the
    // effective geometry/BCs at Kim's own parameters: `params.geometry_n=8`
    // hits the `>= 8` branch below, which builds the exact same sharp
    // rectangle as upstream's intersection() call (a_nd - fabs(x), not a
    // rounded pow() formula), and because a_nd (=1, by construction) exceeds
    // the domain box's own half-width (0.5), the x-constraint in that
    // rectangle is never actually binding inside the box either way — so the
    // x-walls are, in practice, upstream's plain box edges here too. No
    // discrepancy vs. upstream survives this check.
    // Superellipse solid: |x/a|^n + |y/b|^n = 1
    // n >= 8 → perfect rectangle (avoids pow() singularities at sharp corners)
    #if EMBED
    if (params.geometry_n >= 8.)
      solid(cs, fs, intersection(a_nd - fabs(x), b_nd - fabs(y)));
    else
      solid(cs, fs, 1. - pow(fabs(x/a_nd), params.geometry_n)
                       - pow(fabs(y/b_nd), params.geometry_n));
    #endif
  }

  // [PROJECT ADDED, 2026-09-15] Kept from the bisection that found the
  // stracers reset (diary.md 2026-09-15 (8)): post_restore vs post_init
  // brackets everything event init does to the restored state, which is
  // where the tracer loss actually lived. Cheap, and it is the pair that
  // would have caught it years earlier.
  if (restart_file)
    write_restart_diagnostic ("post_init");
}


#if REPLAY_FLUX
// iforce.h (through tension.h) clamps f to [0,1] in its event acceleration, every step, after the VOF step
// and the scalar diffusion and before the next VOF step. VOF round-off leaves f ~ -1e-17 or 1+1e-17 in
// some cells, and the oxygen reset below tests f == 0 exactly, so a replay that skips this clamp resets
// a different set of gas cells (measured: the pure-gas cell count drifts by 3 per step, and oxy with it).
// It must run before the other user events of the next step, so it is registered first among them.
event clamp_f (i++)
{
  foreach()
    f[] = clamp (f[], 0., 1.);
}
#endif
