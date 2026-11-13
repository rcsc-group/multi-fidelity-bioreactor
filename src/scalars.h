// scalars.h -- Tracer and oxygen release (events tracer, oxygen).
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

// Bulk-liquid tracer sums (cells with f >= 0.5), diary 2026-10-29. c2 is the TOTAL amount per cell (henry_oxy2:
// phi1 = a f alpha, phi2 = a (1-f), a = c/(f alpha + 1 - f)); the liquid concentration phi1/f is bounded only
// away from f -> 0. c2*f (the old c2_liq_sum) counts c f^2 and moves when the free surface crosses a grid row.
// tr_oxy.dat columns 13-15 (b0 b1 b2) and mix_release.dat (the same at the release step) carry these sums.
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

// ================================================================== //
//                           TRACER SETUP                             //
// ================================================================== //
#if TRACER
#if REPLAY_FLUX
// Release by time instead of by event time: restore() runs events(false) up to the restored t, which
// consumes an event whose time equals that t without running it, and a replay that releases at its
// start (n_mix_cycles = 0) has exactly that t. The replay starts at the release step, as the coupled run.
static int tracer_done = 0;
event tracer (i++) {
  if (!tracer_done && t >= t_mix - REC_EPS*(1. + fabs (t))) {
  tracer_done = 1;
#else
event tracer(t = t_mix){
#endif
  // On a CONTINUATION the tracer was injected in an earlier segment and has
  // been mixing ever since; re-injecting would overwrite that partially-mixed
  // field with a fresh segregated one and restart chi from 0 mid-experiment
  // (measured: sigma^2 jumping 0.0079 -> 0.2127, probe job 6410668).
  // An `if` wrapper, NOT an early `return` -- a valueless return in a Basilisk
  // event compiles to a nonzero int return and silently stops the time loop.
  if (!params.restart_continue) {

  double h_tr;
  h_tr = (M_PI*R_tr*R_tr);  // Area of circular tracer patch

  // circular shape tracer at the center of the liquid
  // fraction(c, -(sq(x-x_tr) + sq(y-y_tr) - sq(R_tr)) ); 

  // tracer released as a line (same area)
  // fraction(c, intersection( -(y-y_tr - 0.5*h_tr), -(-(y-y_tr + 0.5*h_tr)) ));

  // [PROJECT ADDED, 2026-09-16] Kim's Fig. 5 compares the degree of mixing for
  // FOUR initial tracer configurations, and Fig. 6 shows their fields at
  // chi = 0 and chi = 0.5. Upstream left the circle and the line commented out
  // and never initialised c1 or c3 at all, so only the top-half case existed
  // and neither figure could be drawn.
  //
  // One run gives all of them: the four soluble tracers are independent
  // scalars advected by the same flow, and tr_oxy.dat already logs the first
  // and second moments of each. Four runs would differ only in noise.
  //
  //   c2  top half   (upstream's VERTICAL_MIXUP; left exactly as it was)
  //   c1  left half
  //   c3  circle at the centre of the water
  //   c   line across the middle of the water, same area as the circle
  //
  // Kim's fourth panel in Fig. 6 is the BOTTOM half, which needs no tracer of
  // its own: sigma^2 is invariant under c -> 1 - c, so the bottom half's chi(t)
  // is identical to the top half's, and its field is 1 - c2 wherever there is
  // water. Spending a fifth advected scalar on it would buy a copy.
  //
  // Gated on EXTRA_TRACERS because advecting them is not free: each costs two
  // VOF-advected fields and a multigrid diffusion solve per timestep, and
  // dropping the three idle ones is where the mixing sweeps got 41% of their
  // speed (fc70305). Production keeps EXTRA_TRACERS=0 and one tracer; the
  // Fig. 5/6 build turns them on.
#if EXTRA_TRACERS
  {
    // The line spans the water and is thinned to the circle's area, but never
    // below one cell: at L6 the exact band is 0.23 cells thick and would
    // initialise NOTHING, giving a tracer that is identically zero and a chi
    // with the right shape and no meaning. Where the floor binds, the line
    // and the circle no longer have equal areas -- postprocess reads each
    // tracer's own sigma^2_max at injection, so chi stays correct, but the
    // two configurations are then not area-matched and should not be compared
    // as if they were.
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

  // Vertical mixing-top side
  #if VERTICAL_MIXUP
  {
    // Midpoint of the actual liquid column — fill-level-aware.
    // y_fill = Ly*(2*fill-1), y_bot = -Ly → midpoint = Ly*(fill-1).
    // The old hardcoded -Ly/4 was only correct for fill=0.5 and left
    // fill<=0.3 with zero tracer coverage.
    double y_liq_mid = Ly * (params.fill_level - 1.0);
    foreach(){
      // f[]>0.5 (majority-liquid) instead of f[]==1 (exact) so that near-interface
      // cells are included — at coarse resolutions f never reaches exactly 1.0 in
      // the top liquid layer, causing c2 to remain 0 throughout the run.
#ifndef SEED_KIM
#define SEED_KIM 1   // default since 2026-10-29 (diary: f > 0.5 put liquid concentration 1/f > 1)
#endif
#if SEED_KIM   // Kim's rule: full liquid cells only (his BioReactor.c event tracer); diary 2026-10-29 seed test
      if ((f[] == 1) && (cs[]==1) && (y >= y_liq_mid))
#else
      if ((f[] > 0.5) && (cs[]==1) && (y >= y_liq_mid))
#endif
        c2[] = 1.0;    // Upper half of liquid
    }
    // Synchronise MPI ghost cells after direct foreach() write.
    // Without this, ranks that border the upper liquid region see stale c2
    // values in neighbour halos at the next tracer_diffusion call, producing
    // a spurious large gradient (c_tracer_alpha=1e30 amplifies beta enormously)
    // that crashes h_relax under MPI FP trapping.  Same fix as boundary({oxy})
    // in the oxygen event.
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


// ================================================================== //
//                           OXYGEN SETUP                             //
// ================================================================== //
#if OXYGEN
// From t_mix on this runs EVERY step: gas cells (f == 0) are reset to oxy = 1, i.e. the gas is a reservoir
// held at saturation. The replay must repeat it each step (it modifies the scalar between advections).
#if REPLAY_FLUX
event oxygen (i++) {
  if (t >= t_mix - REC_EPS*(1. + fabs (t))) {
#else
event oxygen (t=t_mix; i++){
#endif

#if OXYGEN_AIR
  foreach(){
    if ((f[] == 0) && (cs[]==1))
      oxy[] = 1.;     // Oxygen in gas regions only
  }
  // Synchronise MPI ghost cells after direct foreach() write.
  // Without this, ranks that border gas regions see stale oxy=0 in
  // neighbour halos at the first tracer_diffusion call, producing a
  // spurious large gradient that crashes h_relax under MPI FP trapping.
  boundary ({oxy});
#endif
#if REPLAY_FLUX
  }
#endif
}
#endif
