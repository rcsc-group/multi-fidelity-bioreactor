// init.h -- event init: build the bag and the initial fields, or restore them from a checkpoint.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.
//
// A restart is the hard part. Basilisk's dump/restore does not keep everything a fresh start sets up, and
// each missing piece fails silently (wrong numbers, no crash). The steps below put each piece back. The
// reasons, with measurements, are in docs_site/developer/solver-pitfalls.md; the full investigation
// notes are in docs/history/solver-comment-archive.md (internal).

// The bag wall as an embedded boundary: superellipse |x/a|^n + |y/b|^n = 1. For n >= 8 we use a sharp
// rectangle (no pow() at the corners); at Kim et al.'s values (n = 8, a = L_bio) this is their geometry.
// solid() writes leaf cells only (see restore_geometry()).
static void bag_solid (void)
{
  double a_nd = params.geometry_a / L_bio;
  double b_nd = params.geometry_b / L_bio;
  if (params.geometry_n >= 8.)
    solid (cs, fs, intersection(a_nd - fabs(x), b_nd - fabs(y)));
  else
    solid (cs, fs, 1. - pow(fabs(x/a_nd), params.geometry_n)
                      - pow(fabs(y/b_nd), params.geometry_n));
}

// restore() drops the refine/prolongation/restriction attributes of every field: they are set by the
// `defaults (i = 0)` / `metric (i = 0)` events of henry_oxy2.h, centered.h, embed.h and vof.h, and a
// restart never sees i = 0. Without the embed-aware versions, the multigrid spreads NaN from solid cells
// (stracers) and the pressure near the wall drifts until the run fails (u, p, g: SIGFPE about one period
// after the restart). Pitfall "attributes-after-restore". Must run BEFORE any refine().
static void reapply_field_attributes (void)
{
#if TREE
  for (scalar s in stracers) {
    s.refine = refine_embed_linear;
    set_prolongation (s, refine_embed_linear);
    set_restriction (s, restriction_volume_average);
  }
#endif
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
  cs.refine = embed_fraction_refine;
  cs.prolongation = fraction_refine;
  foreach_dimension()
    fs.x.prolongation = embed_face_fraction_refine_x;
#endif // TREE && EMBED
#if TREE
  f.refine = fraction_refine;
  set_prolongation (f, fraction_refine);
#endif // TREE
}

#if CROSS_LEVEL_WARMSTART
// Warm start from a checkpoint written at a COARSER level (-DCROSS_LEVEL_WARMSTART=1, test builds only).
// restore() gives back the tree as dumped; refine() makes it uniform at this run's level.
//   - Use refine(level < L), not adapt_wavelet(): wavelet details are exactly 0 in flat regions, so no
//     tolerance makes adapt_wavelet refine everywhere.
//   - Rebuild cs/fs with solid() BEFORE refine(): refine_embed_linear picks its interpolation from the
//     coarse fs, and after restore() about half of the open faces have fs = 0, so half the cells get
//     piecewise-constant values (a checkerboard; +30% on shear, +200% on dissipation). With this order an
//     L6 -> L7 warm start matches the L7 cold start within 0.05%. Pitfall "fs-after-restore".
//     -DCROSS_LEVEL_PREREFINE_SOLID=0 gives the old, wrong order (for A/B tests).
//   - Check the result: each refined leaf gives 4 children per level, so n_after/n_before must be an
//     exact power of 4 (1 for a same-level restart). Cells outside the bag stay coarse, so the count is
//     not NN*NN.
static void cross_level_refine (void)
{
#ifndef CROSS_LEVEL_PREREFINE_SOLID
#define CROSS_LEVEL_PREREFINE_SOLID 1
#endif
#if EMBED && CROSS_LEVEL_PREREFINE_SOLID
  bag_solid ();
#endif
  assert_fs_sane ("pre-refine");

  long n_before = grid->n;
  refine (level < params.fidelity);
  long n_after = grid->n;

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
// Replay with the scalars on a finer grid than the recorded flow (params.replay_scalar_level). Same order
// as cross_level_refine(): solid() on the coarse grid first, then refine().
static void replay_refine (void)
{
  if (depth() != params.fidelity) {
    fprintf (stderr, "ERROR: replay: restored grid has depth %d, params.fidelity says %d\n",
             depth(), params.fidelity);
    exit (1);
  }
#if EMBED
  bag_solid ();
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
  // the fs = 1 coarse-level convention of a fresh run, also on the new levels (see restore_checkpoint())
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
#endif

// fs (a face field) is not in the dump, and what restore() leaves in it is partly zero. Rebuild cs/fs
// from the geometry (the bag does not move: rocking is a body force), then restrict them, because solid()
// writes leaf cells only and the coarse levels feed the multigrid every step. This is Basilisk's own
// pattern for solid() outside i = 0 (src/test/neumann3D.c). Pitfall "fs-after-restore".
static void restore_geometry (void)
{
#if EMBED
  bag_solid ();
#if TREE && !REPLAY_FLUX      // a replay keeps the coarse levels it set up in restore_checkpoint()
  restriction ({cs, fs});
#endif
  assert_fs_sane ("post-solid");   // a fully open face must not have fs = 0
#endif
}

#if !REPLAY_FLUX
// A restart may change omega_b or theta_max, which changes the velocity unit U_bio. Convert the stored
// u, uf, p, pf to the new unit: u scales with su = U_bio_prev/U_bio, pressure with su^2. (Tests show the
// effect is below the cycle-to-cycle spread even for a 14% change; it is kept because it is the correct
// conversion. -DLEGACY_SU=1 uses the old omega_b-only ratio, for A/B tests.)
// g (pressure gradient + acceleration) has terms with different scalings, so no single factor is right.
// It only feeds the first predictor half-step, and a cold start runs that step with g = 0, so we set
// g = 0, as a cold start does. (Tested: scaling g by su^2, or keeping it, gives the same result.)
static void rescale_to_new_units (void)
{
  if (params.omega_b_prev > 0.) {
    double H_bio_here  = 2.*L_bio*Ly;
    double T_per_prev  = 2.*pi/params.omega_b_prev;
    double Th_max_prev = params.theta_max_prev[0] * pi / 180.;
    double V_bio_prev  = L_bio/4*(H_bio_here + 0.5*L_bio*tan(Th_max_prev));
    double U_bio_prev  = V_bio_prev/(H_bio_here*0.5)/T_per_prev;
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
      pf[]  *= su * su;
    }
    // uf is used directly (advection, VOF) before it is rebuilt from u, so it must be converted too.
    foreach_face()
      uf.x[] *= su;
    boundary ({u.x, u.y, p, pf, uf.x, uf.y});
    foreach() {
      g.x[] = 0.;
      g.y[] = 0.;
    }
    boundary ({g.x, g.y});
  }
}
#endif // !REPLAY_FLUX

event init (t = 0)
{
  if (restart_file) {
    // ---- Restart from a checkpoint (timing was set in set_run_window()) ----
    // p and pf are `nodump` by default, and restore() skips nodump fields on the READING side, whatever
    // the writer did. So un-set the flag around restore(), or p comes back as exactly zero.
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
    // A fresh coupled run sets fs = 1 on all coarse levels once (embed.h, event metric) and never
    // restricts fs after solid(); the scalar multigrid reads those coarse levels. Reproduce that here,
    // and skip the restriction after solid() (restore_geometry()), or the oxygen solve differs (1e-3
    // after one step).
    foreach_face()
      fs.x[] = 1.;
    restriction ({fs});
#endif

    reapply_field_attributes ();
#if CROSS_LEVEL_WARMSTART
    cross_level_refine ();
#endif
#if REPLAY_FLUX
    if (params.replay_scalar_level > params.fidelity)
      replay_refine ();
#endif
    restore_geometry ();

#if REPLAY_FLUX
    // The coupled run applies fractions_cleanup(cs, fs) every step (event adapt, centered.h), so the cs/fs
    // it advected f with are the cleaned ones. One call reaches the same fixed point. (The uf = 0 on
    // closed faces that event adapt also does is already in the recorded uf.)
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
      // Scalars interpolated into cells that the fine geometry calls solid (cs = 0) leak into cut cells
      // and make the oxygen solve diverge in the first step. A fresh run has 0 there.
      foreach()
        if (cs[] == 0.)
          for (scalar s in stracers)
            s[] = 0.;
      boundary (stracers);
      restriction (stracers);
    }
    replay_open ();
#endif

#if !REPLAY_FLUX
    rescale_to_new_units ();
#endif

    // Scalars. A new-condition restart starts them from zero; a continuation (restart_continue = 1) keeps
    // them: they are the experiment. In both cases restriction() refreshes the coarse levels and the MPI
    // ghost cells, which otherwise keep stale checkpoint values and make the multigrid diverge slowly.
    if (!params.restart_continue)
      reset (stracers, 0.);
    boundary (stracers);
    restriction (stracers);
  } else {
    // ---- Fresh start ----
    // Liquid fills params.fill_level of the bag height (0.5 -> interface at y = 0, as in Kim et al.).
    double b_nd = params.geometry_b / L_bio;  // == Ly
    double y_fill = b_nd * (2.*params.fill_level - 1.);
    fraction(f, y_fill - y);
    #if EMBED
    bag_solid ();
    #endif
  }

  // Compare with restart_diagnostic_post_restore.txt: the difference is what event init changed.
  if (restart_file)
    write_restart_diagnostic ("post_init");
}


#if REPLAY_FLUX
// The coupled run clamps f to [0,1] every step (iforce.h, event acceleration). Round-off leaves f at
// -1e-17 or 1+1e-17, and the oxygen reset tests f == 0 exactly, so the replay must clamp too. Registered
// here, first among the user events, so it runs where the coupled run's clamp runs.
event clamp_f (i++)
{
  foreach()
    f[] = clamp (f[], 0., 1.);
}
#endif
