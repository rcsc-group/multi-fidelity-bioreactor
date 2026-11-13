// scalars.h -- release of the tracer and the oxygen (events tracer, oxygen).
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

// Bulk-liquid sums of tracer c2 over cells with f >= 0.5: volume b0, first moment b1, second moment b2 of
// the LIQUID concentration cl = phi1/f. (c2 is the total amount per cell, liquid + gas, so c2*f would
// weight by f^2 and jump when the free surface crosses a grid row.) Columns b0 b1 b2 of tr_oxy.dat, and
// mix_release.dat at the release step.
static void mix_bulk_sums (double * b0, double * b1, double * b2)
{
  double s0 = 0., s1 = 0., s2 = 0.;
  foreach (reduction(+:s0) reduction(+:s1) reduction(+:s2))
    if (f[] >= 0.5) {
      double cl = c2[]/(f[]*c2.alpha + (1. - f[]))*c2.alpha;   // liquid concentration phi1/f
      s0 += dv()*f[];
      s1 += dv()*cl*f[];
      s2 += dv()*sq(cl)*f[];
    }
  *b0 = s0, *b1 = s1, *b2 = s2;
}

// ---------------------------------------------------------------------------------------------------
// event tracer: at t_mix, set the initial tracer fields (concentration 1 in a region of the liquid).
//   c2  top half of the liquid (the production tracer: mixing time)
//   c1  left half;  c3  circle at the centre;  c  horizontal line, same area as the circle
// (Kim et al., Figs. 5-6. The bottom half needs no tracer: the variance is the same under c -> 1 - c.)
// c, c1, c3 only with EXTRA_TRACERS=1; each costs two advected fields and a diffusion solve per step.
// ---------------------------------------------------------------------------------------------------
#if TRACER
#if REPLAY_FLUX
// A replay releases by time, not by event time: restore() marks an event at the restored t as done
// without running it, and a replay that starts at the release (n_mix_cycles = 0) has exactly that t.
static int tracer_done = 0;
event tracer (i++) {
  if (!tracer_done && t >= t_mix - REC_EPS*(1. + fabs (t))) {
  tracer_done = 1;
#else
event tracer(t = t_mix){
#endif
  // A continuation (restart_continue = 1) keeps the tracer of the previous segment: no new release.
  // An `if`, not an early `return`: a bare return in an event stops the time loop (pitfall "bare-return").
  if (!params.restart_continue) {

  double h_tr;
  h_tr = (M_PI*R_tr*R_tr);  // circle area = line area

#if EXTRA_TRACERS
  {
    // The line is at least one cell thick (at L6 the exact band is 0.23 cells and would set nothing).
    // Where that floor applies, the line and the circle do not have the same area.
    foreach() {
      if ((f[] > 0.5) && (cs[] == 1)) {
        double half_h = fmax (0.5 * h_tr, Delta);
        if (x <= x_tr)
          c1[] = 1.0;                             // left half
        if (sq(x - x_tr) + sq(y - y_tr) <= sq(R_tr))
          c3[] = 1.0;                             // circle at the centre
        if (fabs (y - y_tr) <= half_h)
          c[]  = 1.0;                             // line across the middle
      }
    }
    boundary ({c, c1, c3});
  }
#endif

  #if VERTICAL_MIXUP
  {
    // Top half of the liquid column, for any fill level: y_fill = Ly(2 fill - 1), bottom -Ly.
    double y_liq_mid = Ly * (params.fill_level - 1.0);
    foreach(){
      // SEED_KIM (default): full liquid cells only (f == 1), as Kim et al. do. SEED_KIM=0 also seeds
      // interface cells (f > 0.5), which puts a liquid concentration 1/f > 1 there.
#ifndef SEED_KIM
#define SEED_KIM 1
#endif
#if SEED_KIM
      if ((f[] == 1) && (cs[]==1) && (y >= y_liq_mid))
#else
      if ((f[] > 0.5) && (cs[]==1) && (y >= y_liq_mid))
#endif
        c2[] = 1.0;
    }
    // Refresh MPI ghost cells after writing c2 directly: a stale halo gives a huge gradient (alpha = 1e30)
    // in the next diffusion solve, which fails under MPI floating-point traps.
    boundary ({c2});
  }
  #endif
  }  // end if (!params.restart_continue)
  if (REPLAY_FLUX || !params.restart_continue) {
    double b0, b1, b2;
    mix_bulk_sums (&b0, &b1, &b2);
    if (pid() == 0) {
      FILE * fr = fopen ("mix_release.dat", "w");
      fprintf (fr, "i t b0 b1 b2\n%d %.17g %.17g %.17g %.17g\n", i, t, b0, b1, b2);
      fclose (fr);
    }
  }

#if REPLAY_FLUX
  }
#endif
}
#endif

// ---------------------------------------------------------------------------------------------------
// event oxygen: from t_mix on, EVERY step, set oxy = 1 in pure gas cells: the air is a reservoir at
// saturation. (A replay must repeat this every step too: it changes the scalar between advections.)
// ---------------------------------------------------------------------------------------------------
#if OXYGEN
#if REPLAY_FLUX
event oxygen (i++) {
  if (t >= t_mix - REC_EPS*(1. + fabs (t))) {
#else
event oxygen (t=t_mix; i++){
#endif

#if OXYGEN_AIR
  foreach(){
    if ((f[] == 0) && (cs[]==1))
      oxy[] = 1.;
  }
  boundary ({oxy});   // refresh MPI ghost cells, as for c2 above
#endif
#if REPLAY_FLUX
  }
#endif
}
#endif
