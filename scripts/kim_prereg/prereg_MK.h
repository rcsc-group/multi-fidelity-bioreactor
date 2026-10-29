/**
# Offline-vs-DNS test (pre-registered in diary.md, 2026-10-29). Not part of Kim's code.

Producer (codes_bio_run/BioReactor.c), -DPREREG_REC=1: from the release t_mix, for PREREG_NPER periods,
at the start of every step (event vof overload, before the scalar advection) it writes per rank
prereg_rec.<pid>.bin: per step int i, double t, double dt, long n, then for each local leaf
u.x, u.y, f, uf.x[], uf.x[1], uf.y[], uf.y[0,1] (PRE_NPC doubles). The record idiom is copied from
src/BioReactor.c of multi-fidelity-bioreactor (rec_gather / rec_scatter).

Offline (codes/TracerTest_test.c), -DPREREG_MODE=m (0 = Kim's code unchanged):
  1 exact : load u, f, uf of every record step
  2 S     : load u, f of every record step; uf rebuilt as Kim's line uf.x[] = fm.x[]*face_value(u.x,0)
  3 H     : load u, f, uf only at the first step of each of PREREG_HOLD snapshots per period; hold between
  4 W     : as 1, but only the record steps of the first period, looped until PREREG_NPER periods
PREREG_FRESET = 1: f is overwritten every step with the loaded (in mode 3: held) f, as Kim's Field_define
does; 0: f is loaded once at the first step, then advected by the solver. uf is put in place every step.
dt of each step = the producer's dt of that record step.

Both, -DPREREG_DIAG=1: prereg_diag.dat, one line per step at the same point (start of the step):
  i t dt S0 S1 S2 O   with S0 = sum dv f, S1 = sum dv c2 f, S2 = sum dv c2^2 f, O = sum dv oxy f.
In the offline run, t is the record's time, so the two files share one time axis.
*/

#ifndef PREREG_REC
#define PREREG_REC 0
#endif
#ifndef PREREG_DIAG
#define PREREG_DIAG 0
#endif
#ifndef PREREG_MODE
#define PREREG_MODE 0
#endif
#ifndef PREREG_NPER
#define PREREG_NPER 3
#endif
#ifndef PREREG_HOLD
#define PREREG_HOLD 100
#endif
#ifndef PREREG_SHOLD
#define PREREG_SHOLD 0
#endif
#ifndef PREREG_FRESET
#define PREREG_FRESET 1
#endif
#ifndef PREREG_DIR
#define PREREG_DIR "."
#endif
#define PRE_NPC 7
#define PRE_EPS 1e-9

static FILE * pre_fp = NULL;

static long pre_leaf_count (void)
{
  long n = 0;
  foreach_cell() {
    if (is_leaf(cell)) {
      if (is_local(cell))
        n++;
      continue;
    }
  }
  return n;
}

static void pre_name (char * name, size_t len)
{
  snprintf (name, len, "%s/prereg_rec.%d.bin", PREREG_DIR, pid());
}

#if PREREG_DIAG
static void pre_diag (double tt)
{
  // S1, S2, O: f-weighted sums as in his normcal (c*f). C1, C2, OL: c is the TOTAL amount per cell
  // (henry_oxy2_MK.h: phi1 = a*f*alpha, phi2 = a*(1-f), a = c/(f*alpha + 1 - f)); liquid part = sum phi1 dv.
  // Same code as src/BioReactor.c pre_diag.
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

/** ## Producer side */
#if !defined(PREREG_OFFLINE) && (PREREG_REC || PREREG_DIAG)
static int pre_state = 0;   // 0 before the window, 1 recording, 2 done

event vof (i++)
{
  if (PREREG_REC && pre_state == 0 && t >= t_mix - PRE_EPS*(1. + fabs (t))) {
    char name[600];
    pre_name (name, sizeof name);
    pre_fp = fopen (name, "wb");
    if (!pre_fp) {
      fprintf (stderr, "ERROR: prereg: cannot open %s\n", name);
      exit (1);
    }
    pre_state = 1;
  }
  if (pre_state == 1 && t > t_mix + PREREG_NPER*T_per_st + PRE_EPS*(1. + fabs (t))) {
    fclose (pre_fp);
    pre_fp = NULL;
    pre_state = 2;
  }
  if (pre_state == 1) {
    long n = pre_leaf_count();
    double * buf = (double *) malloc (n*PRE_NPC*sizeof(double));
    double * w = buf;
    foreach_cell() {
      if (is_leaf(cell)) {
        if (is_local(cell)) {
          *w++ = u.x[];
          *w++ = u.y[];
          *w++ = f[];
          *w++ = uf.x[];
          *w++ = uf.x[1];
          *w++ = uf.y[];
          *w++ = uf.y[0,1];
        }
        continue;
      }
    }
    int ii = i;
    double tt = t, hdt = dt;
    fwrite (&ii, sizeof(int), 1, pre_fp);
    fwrite (&tt, sizeof(double), 1, pre_fp);
    fwrite (&hdt, sizeof(double), 1, pre_fp);
    fwrite (&n, sizeof(long), 1, pre_fp);
    fwrite (buf, sizeof(double), n*PRE_NPC, pre_fp);
    fflush (pre_fp);
    free (buf);
  }
#if PREREG_DIAG
  if (t >= t_mix - PRE_EPS*(1. + fabs (t)) &&
      t <= t_mix + PREREG_NPER*T_per_st + PRE_EPS*(1. + fabs (t)))
    pre_diag (t);
#endif
}
#endif // producer

/** ## Offline side */
#if defined(PREREG_OFFLINE) && PREREG_MODE
static long pre_pos0 = 0, pre_nsteps = 0, pre_k = 0, pre_lap = 0;
static double pre_t0 = 0., pre_tlast = 0.;
static int pre_stop = 0;
static long pre_snap = -1;
scalar pre_fh[];
face vector pre_ufh[];

static int pre_head (int * hi, double * ht, double * hdt, long * hn)
{
  return fread (hi, sizeof(int), 1, pre_fp) == 1 && fread (ht, sizeof(double), 1, pre_fp) == 1 &&
         fread (hdt, sizeof(double), 1, pre_fp) == 1 && fread (hn, sizeof(long), 1, pre_fp) == 1;
}

static void pre_open (void)
{
  char name[600];
  pre_name (name, sizeof name);
  pre_fp = fopen (name, "rb");
  if (!pre_fp) {
    fprintf (stderr, "ERROR: prereg: cannot open %s\n", name);
    exit (1);
  }
  long ncell = pre_leaf_count(), hn;
  int hi;
  double ht, hdt;
  while (pre_head (&hi, &ht, &hdt, &hn)) {
    if (hn != ncell) {
      fprintf (stderr, "ERROR: prereg: record has %ld cells per step on rank %d, the grid has %ld\n",
               hn, pid(), ncell);
      exit (1);
    }
    if (pre_nsteps == 0)
      pre_t0 = ht;
    pre_tlast = ht;
    pre_nsteps++;
    fseek (pre_fp, (long) (hn*PRE_NPC*sizeof(double)), SEEK_CUR);
  }
  pre_lap = pre_nsteps;
#if PREREG_MODE == 4
  // one period: the steps whose time is before t0 + T
  pre_lap = 0;
  fseek (pre_fp, 0, SEEK_SET);
  while (pre_head (&hi, &ht, &hdt, &hn)) {
    if (ht < pre_t0 + T_per_st - PRE_EPS*(1. + fabs (ht)))
      pre_lap++;
    fseek (pre_fp, (long) (hn*PRE_NPC*sizeof(double)), SEEK_CUR);
  }
#endif
  fseek (pre_fp, 0, SEEK_SET);
  pre_pos0 = 0;
  if (pid() == 0)
    fprintf (stderr, "prereg: mode %d, %ld steps in record, pass of %ld, t0 %.17g, last %.17g, freset %d\n",
             PREREG_MODE, pre_nsteps, pre_lap, pre_t0, pre_tlast, PREREG_FRESET);
}


// Kim's step rule (timestep_MK.h without its growth limiter): CFL_MK * min over faces of Delta*fm/|uf|, capped
// at DT_MK (= DT_MAX). A record step longer than that is taken as n equal substeps with the same fields, so
// the offline clock still lands on every record time.
static double pre_cfl_dt (void)
{
  double dtmax = DT_MK/CFL_MK;
  foreach_face (reduction(min:dtmax))
    if (pre_ufh.x[] != 0. && fm.x[] > 0.) {
      double d = Delta*fm.x[]/fabs (pre_ufh.x[]);
      if (d < dtmax)
        dtmax = d;
    }
  return dtmax*CFL_MK;
}

// Kim's seeding (his event initial_conditions, = the producer's events tracer and oxygen), done on the
// record's first f before the first advection, which is where the producer seeds (its first diagnostic row
// already holds the tracer).
static void pre_seed (void)
{
#if TRACER
#if HORIZONTAL_MIXL
  fraction (c, 0 - x);
  foreach()
    if (f[] < 1 || (cs[] == 0))
      c[] = 0;
#endif
#if HORIZONTAL_MIXR
  fraction (c1, 0 + x);
  foreach()
    if ((cs[] == 0) || (f[] < 1))
      c1[] = 0;
#endif
#if VERTICAL_MIXUP
  foreach()
    if ((f[] == 1) && (cs[] == 1) && (y >= -Ly*0.5*0.5))
      c2[] = 1.0;
#endif
#if VERTICAL_MIXDOWN
  foreach()
    if ((f[] == 1) && (cs[] == 1) && (y <= -Ly*0.5*0.5))
      c3[] = 1.0;
#endif
#endif
#if OXYGEN && OXYGEN_AIR
  foreach()
    oxy[] = ((f[] == 0) && (cs[] == 1)) ? 1. : 0.;
#endif
  boundary ({c, c1, c2, c3, oxy});
}

static int pre_sub = 0, pre_nsub = 1, pre_pass = 0;
static double pre_hdt = 0.;

// Overload of event vof: runs before the solver's own vof event of the same step.
event vof (i++)
{
  if (i == 0) {
    // Kim's geometry, built in his event restore_all / Field_define (disabled in these modes): same lines
    solid (cs, fs, intersection (-(y - 0.5*Ly), -(-y - 0.5*Ly)));
    fractions_cleanup (cs, fs);
    boundary ({cs});
    boundary ((scalar *){fs});
    pre_open();
  }
  if (pre_sub == 0) {
    if (pre_k >= pre_lap) {      // only reached in mode 4: loop to the first record
      fseek (pre_fp, pre_pos0, SEEK_SET);
      pre_k = 0;
      pre_pass++;
    }
    int hi;
    double ht, hdt;
    long hn;
    if (!pre_head (&hi, &ht, &hdt, &hn)) {
      fprintf (stderr, "ERROR: prereg: record ends early at step %ld\n", pre_k);
      exit (1);
    }
    double * buf = (double *) malloc (hn*PRE_NPC*sizeof(double));
    if (fread (buf, sizeof(double), hn*PRE_NPC, pre_fp) != (size_t) (hn*PRE_NPC)) {
      fprintf (stderr, "ERROR: prereg: truncated record at step %ld\n", pre_k);
      exit (1);
    }
    // time check: the offline clock starts at 0, the record at pre_t0 (first pass only)
    if (pre_pass == 0 && fabs (t - (ht - pre_t0)) > PRE_EPS*(1. + fabs (ht))) {
      fprintf (stderr, "ERROR: prereg: offline t=%.17g but record t-t0=%.17g at i=%d\n", t, ht - pre_t0, i);
      exit (1);
    }
#if PREREG_MODE == 3 || (PREREG_MODE == 2 && PREREG_SHOLD)   // PREREG_SHOLD: rebuild from held snapshots (= his)
    long snap = (long) floor ((ht - pre_t0)/(T_per_st/PREREG_HOLD) + PRE_EPS);
    int load = (snap != pre_snap);
    pre_snap = snap;
#else
    int load = 1;
#endif
    // The loaded (or held) fields go to pre_fh / pre_ufh; every step then puts them in place, as Kim's
    // event Field_define puts u_save, ff_save, uf_save in place every step.
    if (load) {
      const double * r = buf;
      foreach_cell() {
        if (is_leaf(cell)) {
          if (is_local(cell)) {
            u.x[] = r[0];
            u.y[] = r[1];
            pre_fh[] = r[2];
            pre_ufh.x[]    = r[3];
            pre_ufh.x[1]   = r[4];
            pre_ufh.y[]    = r[5];
            pre_ufh.y[0,1] = r[6];
            r += PRE_NPC;
          }
          continue;
        }
      }
      boundary ({u.x, u.y, pre_fh});
#if PREREG_MODE == 2
      foreach_face()
        pre_ufh.x[] = fm.x[]*face_value(u.x,0);   // Kim's line (TracerTest_test.c, event restore_all)
#endif
      boundary ((scalar *){pre_ufh});
    }
    free (buf);
    pre_hdt = hdt;
    double dts = pre_cfl_dt();
    // only the rebuilt field (mode 2) can break the CFL limit of the producer's dt; the true flux keeps one step
    pre_nsub = (PREREG_MODE == 2 && hdt > dts*(1. + PRE_EPS)) ? (int) ceil (hdt/dts - PRE_EPS) : 1;
  }
  if (PREREG_FRESET || i == 0) {
    foreach()
      f[] = pre_fh[];
    boundary ({f});
  }
  foreach_face()
    uf.x[] = pre_ufh.x[];
  boundary ((scalar *){uf});
  if (i == 0)
    pre_seed();
#if PREREG_DIAG
  if (pre_sub == 0)
    pre_diag (pre_t0 + t);
#endif
  dt = pre_hdt/pre_nsub;
  tnext = t + dt;
  if (++pre_sub == pre_nsub) {
    pre_sub = 0;
    pre_k++;
#if PREREG_MODE != 4
    if (pre_k == pre_lap)
      pre_stop = 1;
#else
    if (t + dt >= PREREG_NPER*T_per_st - PRE_EPS*(1. + fabs (t)))
      pre_stop = 1;
#endif
  }
  static long pre_nsub_total = 0, pre_nsub_max = 0;
  if (pre_sub == 0) {
    pre_nsub_total += pre_nsub;
    if (pre_nsub > pre_nsub_max)
      pre_nsub_max = pre_nsub;
  }
  if (pre_stop && pid() == 0)
    fprintf (stderr, "prereg: done, %ld record steps, %ld solver steps, max substeps %ld\n",
             pre_k + pre_pass*pre_lap, pre_nsub_total, pre_nsub_max);
}

// Same idiom as src/BioReactor.c event end_replay: the condition keeps the loop alive, the return 1 ends it
// after the last record step.
event pre_end (i++; i < 2000000000) {
  if (pre_stop)
    return 1;
  return 0;
}
#endif // PREREG_MODE

/** Offline, Kim's own workflow (PREREG_MODE 0) with the diagnostic only: the offline clock starts at the first
snapshot, t_ini (= T_INI_FILE_FIRST). Overload of event vof, so it is taken at the same point of the step. */
#if defined(PREREG_OFFLINE) && !PREREG_MODE && PREREG_DIAG
event vof (i++)
{
  pre_diag (t_ini + t);
}
#endif
