/**
# Probe: does anything change the coarse levels of cs (= cm, embed.h) during the run?

Diary 2026-11-23 (profiling, pre-registration "cm restriction"). Appended to a COPY of src/BioReactor.c by
tools/profile/cm_probe/run.sh; never part of src/. At step 1 (after the first fractions_cleanup() of event
adapt, centered.h) the probe stores cs in cs0 on every cell of the tree (foreach_cell, as in
src/test/balance5.c). Each later step until PROBE_STEPS it prints, split by local/remote and leaf/coarse
cells: the max |cs - cs0| and the number of subnormal or non-finite cs values. One line per step on stderr,
prefix "cm_probe". */

#ifndef PROBE_STEPS
# define PROBE_STEPS 300
#endif

scalar cs0[];

event cm_probe (i++) {
  if (i >= 1 && i <= PROBE_STEPS) {
    if (i == 1)
      foreach_cell()
        cs0[] = cs[];
    else {
      // d[0] local leaf, d[1] local coarse, d[2] remote leaf, d[3] remote coarse (ghost/halo cells)
      double d[4] = {0., 0., 0., 0.};
      double nbad[4] = {0., 0., 0., 0.};
      foreach_cell() {
        int k = (is_local(cell) ? 0 : 2) + (is_leaf(cell) ? 0 : 1);
        double e = fabs (cs[] - cs0[]);
        if (e > d[k] || e != e)
          d[k] = e;
        if (!isfinite (cs[]) || (cs[] != 0. && fabs (cs[]) < 2.2250738585072014e-308))
          nbad[k] += 1.;
      }
#if _MPI
      mpi_all_reduce_array (d, MPI_DOUBLE, MPI_MAX, 4);
      mpi_all_reduce_array (nbad, MPI_DOUBLE, MPI_SUM, 4);
#endif
      if (pid() == 0)
        fprintf (stderr, "cm_probe %d %.3e %.3e %.3e %.3e %g %g %g %g\n", i,
                 d[0], d[1], d[2], d[3], nbad[0], nbad[1], nbad[2], nbad[3]);
    }
  }
}
