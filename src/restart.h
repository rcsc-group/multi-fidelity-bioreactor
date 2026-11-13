// restart.h -- Restart checks: assert_fs_sane() and write_restart_diagnostic().
// Part of BioReactor.c: included at one fixed place; the place sets the event order.


// [PROJECT ADDED, 2026-09-07] Restart-correctness diagnostic (Test C,
// diary.md 2026-09-07): field-level summary statistics (Basilisk's
// statsf(), already used elsewhere in this file, MPI-safe), written to
// restart_diagnostic.txt. Called BOTH right before dump() (this run's own
// pre-dump state, in event dump_checkpoint below) and immediately after
// restore() succeeds on any restart (this run's post-restore state, in
// event init below) -- BEFORE any of the restart-specific fixups (fs
// recompute, omega_b rescaling, stracer reset) so it reflects exactly
// what restore() itself reconstructed, unmodified by anything else.
// Comparing the two files across a dump-then-restore pair (source run's
// file vs the restarted run's file) checks whether restore() faithfully
// reconstructs dump()'s content, independent of what the ramp/forcing
// does afterward -- the question the missing-dAk/dt-term and Delta_theta
// investigations couldn't answer on their own.
// [PROJECT ADDED, 2026-09-14] Invariant guard for the cross-level bug of
// 2026-09-14 (diary.md): fs does not survive dump/restore -- Basilisk's
// dump_list() skips face fields outright (output.h:1037) -- and what comes
// back has ~50% of FULLY OPEN faces at hard zero. Since refine_embed_linear
// gates its interpolation branches on fs truthiness, that silently degrades
// half the domain to piecewise-constant injection during a cross-level
// refine(), costing +30-42% on tau and +126-220% on dissipation with no
// crash, no warning, and a physically plausible-looking result.
//
// The invariant: a face whose BOTH neighbouring cells are entirely fluid
// cannot have zero open area. Violations mean fs is stale/partial and must
// not be used. Cheap (one pass), so it runs on every restart rather than
// behind a debug flag -- this class of bug is silent, which is exactly why
// it needs an always-on check.
static void assert_fs_sane (const char *tag)
{
#if EMBED
  int nopen = 0, nzero = 0;
  foreach_face(reduction(+:nopen) reduction(+:nzero))
    if (cs[] > 0.99 && cs[-1] > 0.99) {
      nopen++;
      if (fs.x[] == 0.) nzero++;
    }
  if (nzero > 0) {
    if (pid() == 0)
      fprintf (stderr, "ERROR [%s]: fs is not usable -- %d of %d fully-open "
               "faces (both neighbours cs>0.99) have fs==0. Face fields do "
               "not survive dump/restore; recompute cs/fs with solid() "
               "BEFORE any refine() that interpolates u/p/g. "
               "See diary.md 2026-09-14.\n", tag, nzero, nopen);
    exit (1);
  }
  if (pid() == 0)
    fprintf (ferr, "fs sanity [%s]: OK (%d fully-open faces, none zero)\n",
             tag, nopen);
#endif
}

#if REPLAY_FLUX
// Replay has no u, p, g: only the fields it advances.
static void write_restart_diagnostic (const char *tag)
{
  stats sfv  = statsf (f);
  stats sc2  = statsf (c2);
  stats soxy = statsf (oxy);
#if _MPI
  if (pid() == 0)
#endif
  {
    char fname[64];
    sprintf (fname, "restart_diagnostic_%s.txt", tag);
    FILE *fp = fopen (fname, "w");
    if (fp) {
      fprintf (fp, "tag %s\nt %.17g\n", tag, t);
      fprintf (fp, "f_min %.17g\nf_max %.17g\nf_sum %.17g\nf_stddev %.17g\n",
               sfv.min, sfv.max, sfv.sum, sfv.stddev);
      fprintf (fp, "c2_min %.17g\nc2_max %.17g\nc2_sum %.17g\nc2_stddev %.17g\n",
               sc2.min, sc2.max, sc2.sum, sc2.stddev);
      fprintf (fp, "oxy_min %.17g\noxy_max %.17g\noxy_sum %.17g\noxy_stddev %.17g\n",
               soxy.min, soxy.max, soxy.sum, soxy.stddev);
      fclose (fp);
    }
  }
}
#else
static void write_restart_diagnostic (const char *tag)
{
  stats su   = statsf (u.x);
  stats sv   = statsf (u.y);
  stats sp   = statsf (p);
  stats sfv  = statsf (f);
  // [PROJECT ADDED, 2026-09-07] Extended past u.x/u.y/p/f (which all
  // matched exactly in the first pass) to every OTHER field the su-
  // rescaling code (event init's restart branch, a few hundred lines
  // below) itself flags as needing special restart handling: pf (half-
  // step pressure), g.x/g.y (BCG pressure-gradient/acceleration term),
  // uf.x/uf.y (face-centered velocity). The kick test (diary.md
  // 2026-09-07) ruled out the continuous PDE itself as the source of the
  // escape-to-a-different-branch behavior -- a direct, tiny, non-
  // checkpoint velocity perturbation causes no drift at all -- so
  // whatever's wrong must be a field the FIRST diagnostic pass didn't
  // check.
  stats spf  = statsf (pf);
  stats sgx  = statsf (g.x);
  stats sgy  = statsf (g.y);
  // [PROJECT ADDED, 2026-09-15] The passive tracer and oxygen were NOT
  // covered here, and that gap is exactly why their loss across a restart
  // went unnoticed: every field this diagnostic DID track round-tripped to
  // full double precision, so the restart looked clean. Probe 6410668
  // measured c2 mass 0.0534 -> 0 and oxy 0.0803 -> 0 across a seam while
  // u/p/f/g were bit-exact. A guard must cover the fields the experiment
  // depends on, not only the ones a previous investigation happened to
  // suspect. Fig 9/10 (mixing) and Fig 11/12 (kLa) all hinge on these.
  stats sc2  = statsf (c2);
  stats soxy = statsf (oxy);
#if _MPI
  if (pid() == 0)
#endif
  {
    char fname[64];
    sprintf (fname, "restart_diagnostic_%s.txt", tag);
    FILE *fp = fopen (fname, "w");
    if (fp) {
      fprintf (fp, "tag %s\n", tag);
      fprintf (fp, "t %.17g\n", t);
      fprintf (fp, "ux_min %.17g\nux_max %.17g\nux_sum %.17g\nux_stddev %.17g\n",
               su.min, su.max, su.sum, su.stddev);
      fprintf (fp, "uy_min %.17g\nuy_max %.17g\nuy_sum %.17g\nuy_stddev %.17g\n",
               sv.min, sv.max, sv.sum, sv.stddev);
      fprintf (fp, "p_min %.17g\np_max %.17g\np_sum %.17g\np_stddev %.17g\n",
               sp.min, sp.max, sp.sum, sp.stddev);
      fprintf (fp, "f_min %.17g\nf_max %.17g\nf_sum %.17g\nf_stddev %.17g\n",
               sfv.min, sfv.max, sfv.sum, sfv.stddev);
      fprintf (fp, "pf_min %.17g\npf_max %.17g\npf_sum %.17g\npf_stddev %.17g\n",
               spf.min, spf.max, spf.sum, spf.stddev);
      fprintf (fp, "gx_min %.17g\ngx_max %.17g\ngx_sum %.17g\ngx_stddev %.17g\n",
               sgx.min, sgx.max, sgx.sum, sgx.stddev);
      fprintf (fp, "gy_min %.17g\ngy_max %.17g\ngy_sum %.17g\ngy_stddev %.17g\n",
               sgy.min, sgy.max, sgy.sum, sgy.stddev);
      fprintf (fp, "c2_min %.17g\nc2_max %.17g\nc2_sum %.17g\nc2_stddev %.17g\n",
               sc2.min, sc2.max, sc2.sum, sc2.stddev);
      fprintf (fp, "oxy_min %.17g\noxy_max %.17g\noxy_sum %.17g\noxy_stddev %.17g\n",
               soxy.min, soxy.max, soxy.sum, soxy.stddev);
      fclose (fp);
    }
  }
}
#endif // REPLAY_FLUX
