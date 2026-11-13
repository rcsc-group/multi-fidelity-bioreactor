// replay_events.h -- record/replay hooks on the VOF step (event vof, event rec_checkpoint).
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

// Actions of one event name run in reverse order of registration, so our `event vof` runs BEFORE the VOF
// actions of henry_oxy2.h and vof.h: it sees uf and dt exactly as they will be used. (qcc can place an
// event relative to a library event only by registration order, hence the include at this place.)
#if RECORD_FLUX
// End of the step that lands on t_rec0: registered after the release events, so a release at t_rec0 is in
// the checkpoint.
event rec_checkpoint (t = t_rec0)
{
  rec_start();
  rec_state = 1;
}
#endif

// -DPREREG_DIAG=1: per-step sums from the release on, with the same code as on the Kim et al. solver side
// (kim_pipeline/prereg_src/prereg_MK.h), so that both solvers are sampled at the same steps
// (tr_oxy.dat samples only every t_out).
#ifndef PREREG_DIAG
#define PREREG_DIAG 0
#endif
#if PREREG_DIAG
static void pre_diag (double tt)
{
  // S1, S2, O: the f-weighted sums of Kim's normcal (c*f). C1, C2, OL: c is the TOTAL amount per cell
  // (henry_oxy2.h splits it as phi1 = a*f*alpha, phi2 = a*(1-f), a = c/(f*alpha + 1 - f)), so the tracer in
  // the liquid is sum phi1 dv; C2 = sum (phi1/f)^2 f dv is the liquid concentration variance moment.
  double s0 = 0., s1 = 0., s2 = 0., so = 0., c1 = 0., c2s = 0., ol = 0., b0 = 0., b1 = 0., b2 = 0.;
  foreach (reduction(+:s0) reduction(+:s1) reduction(+:s2) reduction(+:so)
           reduction(+:c1) reduction(+:c2s) reduction(+:ol)
           reduction(+:b0) reduction(+:b1) reduction(+:b2)) {
    s0 += dv()*f[];
    s1 += dv()*c2[]*f[];
    s2 += dv()*sq(c2[])*f[];
    so += dv()*oxy[]*f[];
    double ac = c2[]/(f[]*c2.alpha + (1. - f[]));
    double ao = oxy[]/(f[]*oxy.alpha + (1. - f[]));
    c1 += dv()*ac*f[]*c2.alpha;
    c2s += f[] >= 1e-6 ? dv()*sq(ac*c2.alpha)*f[] : 0.;   // phi1/f blows up as f -> 0
    ol += dv()*ao*f[]*oxy.alpha;
    if (f[] >= 0.5) {           // bulk liquid: concentration phi1/f is bounded here (B = chi without f -> 0)
      b0 += dv()*f[];
      b1 += dv()*ac*f[]*c2.alpha;
      b2 += dv()*sq(ac*c2.alpha)*f[];
    }
  }
  static FILE * fpd = NULL;
  if (!fpd && pid() == 0) {
    fpd = fopen ("prereg_diag.dat", "w");
    fprintf (fpd, "i t dt S0 S1 S2 O C1 C2 OL B0 B1 B2\n");
  }
  if (fpd) {
    fprintf (fpd, "%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
             iter, tt, dt, s0, s1, s2, so, c1, c2s, ol, b0, b1, b2);
    fflush (fpd);
  }
}
#endif

#if RECORD_FLUX || REPLAY_FLUX
event vof (i++)
{
#if RECORD_FLUX
  if (rec_state == 1) {
    if (t <= t_rec_end + REC_EPS*(1. + fabs (t))) {
#if REC_CHECKSUM
      rec_checksum();
#endif
      rec_write_step();
    }
    else {
      fclose (fp_rec);
      fp_rec = NULL;
      rec_state = 2;
    }
  }
#endif
#if REPLAY_FLUX
  replay_load_step();
#endif
#if PREREG_DIAG
  if (t >= t_mix - REC_EPS*(1. + fabs (t)))
    pre_diag (t);
#endif
}
#endif
