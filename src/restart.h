// restart.h -- restart checks: assert_fs_sane() and write_restart_diagnostic().
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

// Stop if the embedded face fractions fs are not usable. A face between two cells that are both fully
// fluid must be open (fs > 0). After restore() about half of such faces are 0 (face fields are not in the
// dump), and refine_embed_linear, which tests fs, then silently degrades to piecewise-constant values:
// +30% on shear stress, +200% on dissipation, no crash. The check costs one pass, so it runs on every
// restart. Pitfall "fs-after-restore" in docs_site/developer/solver-pitfalls.md.
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
               "See docs_site/developer/solver-pitfalls.md.\n", tag, nzero, nopen);
    exit (1);
  }
  if (pid() == 0)
    fprintf (ferr, "fs sanity [%s]: OK (%d fully-open faces, none zero)\n",
             tag, nopen);
#endif
}

// restart_diagnostic_<tag>.txt: min, max, sum and standard deviation of every field a restart must carry.
// Written before the dump ("pre_dump", checkpoint.h), right after restore() ("post_restore") and at the end
// of event init ("post_init"). pre_dump of one run vs post_restore of the next tests restore() alone;
// post_restore vs post_init shows what event init changed. Cover every field the experiment depends on:
// the scalars c2 and oxy were once lost across restarts while u, p and f round-tripped exactly.
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
  stats spf  = statsf (pf);
  stats sgx  = statsf (g.x);
  stats sgy  = statsf (g.y);
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
