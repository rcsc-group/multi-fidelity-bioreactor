// replay.h -- Flux record/replay: file format, write (-DRECORD_FLUX) and read (-DREPLAY_FLUX).
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

//                  FLUX RECORD / REPLAY (see top of file)            //
// ================================================================== //
#if RECORD_FLUX || REPLAY_FLUX
#if REPLAY_FLUX && (VIDEOS || CROSS_LEVEL_WARMSTART)
#error "REPLAY_FLUX has no velocity field: VIDEOS and CROSS_LEVEL_WARMSTART are not available"
#endif

// Record format, one file per rank (<rec_dir>/flux_rec.<pid>.bin), read by scripts/replay_gate.py:
//   file header: char magic[8], int pid, int npe
//   per step   : int i, double t, double dt, long long n, then n*REC_NPC doubles: for each leaf cell and
//                each axis, the left and the right face value of uf (x: uf.x[], uf.x[1]; y: uf.y[], uf.y[0,1]).
// The right/top faces are stored too: on the domain edge they live in a ghost cell, which a leaf-only
// traversal would not reach.
#define REC_NPC (2*dimension)
// The cells a record is made of. Recording: the leaves of the flow grid. Replaying: the cells of the flow
// level params.fidelity, which are the leaves too unless the scalars live on a finer grid
// (params.replay_scalar_level > params.fidelity), where they are the parents of the leaves.
#if REPLAY_FLUX
#define REC_AT_CELL (level == params.fidelity)
#else
#define REC_AT_CELL is_leaf(cell)
#endif
#define REC_EPS 1e-9                        // relative tolerance when comparing times
// -DREC_KIM_FLUX=1 (with -DRECORD_FLUX; default 0): store ufk = fm*face_value(u) instead of uf. It emulates
// an offline solver that rebuilds the advecting face velocity from the stored cell-centred velocity (not
// divergence-free), to measure the error of that approach. Same record format; uf, u and the run are untouched.
// -DREC_KIM_MASK=1: also set ufk = 0 on faces with fs = 0, as Kim's Field_define does (diary 2026-10-29).
#ifndef REC_KIM_MASK
#define REC_KIM_MASK 0
#endif
#ifndef REC_KIM_FLUX
#define REC_KIM_FLUX 0
#endif
#if REC_KIM_FLUX && !RECORD_FLUX
#error "REC_KIM_FLUX needs RECORD_FLUX"
#endif
#if REC_KIM_FLUX && dimension != 2
#error "REC_KIM_FLUX is 2D only"
#endif
#if REC_KIM_FLUX
face vector ufk[];                          // declared once at file scope (see the note at tau_field)
#endif
static const char rec_magic[8] = {'F','L','X','R','E','C','0','1'};
static double t_rec0 = -1., t_rec_end = 1e30;   // record window; set in main()
static FILE * fp_rec = NULL;

// Leaf cells owned by this rank. foreach_cell() is serial, so the gather/scatter order is the traversal
// order and identical in the recording and the replaying run (same restored tree).
static long long rec_leaf_count (void)
{
  long long n = 0;
  foreach_cell() {
    if (REC_AT_CELL) {
      if (is_local(cell))
        n++;
      continue;
    }
  }
  return n;
}

static void rec_file_name (char * name, size_t len)
{
  snprintf (name, len, "%s/flux_rec.%d.bin", params.rec_dir, pid());
}

static void rec_gather (double * buf)
{
  double * w = buf;
  foreach_cell() {
    if (REC_AT_CELL) {
      if (is_local(cell))
        foreach_dimension() {
#if REC_KIM_FLUX
          *w++ = ufk.x[];
          *w++ = ufk.x[1];
#else
          *w++ = uf.x[];
          *w++ = uf.x[1];
#endif
        }
      continue;
    }
  }
}

static void rec_scatter (const double * buf)
{
  const double * r = buf;
  foreach_cell() {
    if (REC_AT_CELL) {
      if (is_local(cell))
        foreach_dimension() {
          uf.x[]  = *r++;
          uf.x[1] = *r++;
        }
      continue;
    }
  }
}

// -DREC_CHECKSUM=1: weighted sums of f, uf, c2 and oxy at the same point of the step (before the VOF
// advection, after uf is in place) in both modes, into rec_checksum.txt, to find the first step at which a
// replay leaves the recorded run.
#ifndef REC_CHECKSUM
#define REC_CHECKSUM 0
#endif
#if REC_CHECKSUM
static void rec_checksum_to (const char * file)
{
  double cf = 0., cu = 0., cc2 = 0., co = 0., cog = 0., cogn = 0., nz = 0., ntiny = 0., nneg = 0., cf2 = 0.,
    mc2 = 0., mo = 0., mf = 0.;
  foreach (reduction(+:cf) reduction(+:cu) reduction(+:cc2) reduction(+:co) reduction(+:cog) reduction(+:cogn)
           reduction(+:nz) reduction(+:ntiny) reduction(+:nneg) reduction(+:cf2) reduction(+:mc2) reduction(+:mo)
           reduction(+:mf)) {
    // total (liquid + gas) mass of c2 and oxy: the transported variable c = phi1 + phi2 integrates to it
    mc2 += dv()*c2[];
    mo  += dv()*oxy[];
    mf  += dv()*f[];
    cf2 += f[]*f[];
    nz += (f[] == 0.);
    ntiny += (f[] != 0. && fabs (f[]) < 1e-9);
    nneg += (f[] < 0.);
    double w = 1. + 0.37*x + 0.71*y;
    cf  += f[]*w;
    cu  += (uf.x[] + 2.*uf.y[])*w;
    cc2 += c2[]*w;
    co  += oxy[]*w;
    if (f[] == 0.) {              // pure gas cells, and how many of them
      cog  += oxy[]*w;
      cogn += cs[];
    }
  }
  static FILE * fps[2] = {NULL, NULL};
  FILE ** pf_ = file[12] == '_' ? &fps[1] : &fps[0];   // "rec_checksum_norm.txt" vs "rec_checksum.txt"
  FILE * fp = *pf_;
  if (!fp && pid() == 0)
    fp = *pf_ = fopen (file, "w");
  if (fp) {
    fprintf (fp, "%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %g %g %g %.17g %.17g %.17g %.17g\n", iter, t, cf, cu, cc2, co, cog, cogn, nz, ntiny, nneg, cf2, mc2, mo, mf);
    fflush (fp);
  }
}
static void rec_checksum (void) { rec_checksum_to ("rec_checksum.txt"); }
// where the largest |oxy| is (serial only): position, f, cs, oxy
static void rec_where_max (void)
{
  double big = 0., bx = 0., by = 0., bf = 0., bcs = 0., bo = 0.;
  foreach (serial) {
    if (fabs (oxy[]) > big) {
      big = fabs (oxy[]); bx = x; by = y; bf = f[]; bcs = cs[]; bo = oxy[];
    }
  }
  fprintf (stderr, "where_max i=%d x=%g y=%g f=%g cs=%g oxy=%g\n", iter, bx, by, bf, bcs, bo);
}
#endif
#endif // RECORD_FLUX || REPLAY_FLUX

#if RECORD_FLUX
// Called from event rec_checkpoint at the START of the step that lands on t_rec0, after the release events
// of that step (user events run before the vof/diffusion/projection events of the same step). The dump is
// the state that step starts from, and that step is the first recorded one.
// restore() ends by marking every event of the restored iteration as done, so a restart resumes at the
// NEXT iteration. To resume AT this step the dump is written with iter - 1 in its header (same t).
// p and pf are forced into the dump as in event dump_checkpoint, so that restore() of it is the ordinary
// restart.
static void rec_start (void)
{
  p.nodump = pf.nodump = false;
  iter--;
  dump (file = "replay_start.dump");
  iter++;
  p.nodump = pf.nodump = true;
  if (pid() == 0)
    fprintf (ferr, "record: replay_start.dump written at i=%d t=%.17g\n", iter, t);
  char name[600];
  rec_file_name (name, sizeof name);
  fp_rec = fopen (name, "wb");
  if (!fp_rec) {
    perror (name);
    exit (1);
  }
  int rank = pid(), np = npe();
  fwrite (rec_magic, 1, 8, fp_rec);
  fwrite (&rank, sizeof(int), 1, fp_rec);
  fwrite (&np, sizeof(int), 1, fp_rec);
}

static void rec_write_step (void)
{
  static int first = 1;
  if (first && pid() == 0)
    fprintf (ferr, "record: first step i=%d t=%.17g dt=%.17g\n", iter, t, dt);
  first = 0;
  long long n = rec_leaf_count();
  double * buf = (double *) malloc (n*REC_NPC*sizeof(double));
#if REC_KIM_FLUX
  // foreach_face() also visits the face on the upper domain edge (it lives in a ghost cell), computed from
  // the ghost value of u; boundary() then fills the faces that sit in the halos of the other ranks and of
  // coarser neighbours, so that ufk.x[1] is defined the same way uf.x[1] is.
  // noauto (hence the explicit 2D order x, y): without it qcc calls boundary() on u first, which refreshes
  // the ghost and halo cells of u and changes the coupled run (measured: i = 200 on, L6, 2 ranks). With it,
  // u is only read, as the solver holds it; the faces that reach a ghost or halo cell of u (~0.5 % of them
  // at the first recorded step) use the value u has there.
  foreach_face (x, y, noauto)
    ufk.x[] = (REC_KIM_MASK && !fs.x[]) ? 0. : fm.x[]*face_value (u.x, 0);   // mask: Kim TracerTest_test.c:930
  boundary ((scalar *){ufk});
#endif
  rec_gather (buf);
  int hi = iter;
  double ht = t, hdt = dt;
  fwrite (&hi, sizeof(int), 1, fp_rec);
  fwrite (&ht, sizeof(double), 1, fp_rec);
  fwrite (&hdt, sizeof(double), 1, fp_rec);
  fwrite (&n, sizeof(long long), 1, fp_rec);
  fwrite (buf, sizeof(double), n*REC_NPC, fp_rec);
  fflush (fp_rec);
  free (buf);
}

static int rec_state = 0;   // 0 before the window, 1 recording, 2 window closed

// Empty events: they only make dtnext() land a step exactly on both ends of the window, so that the
// checkpoint is at t_rec0 and the last recorded step ends at t_rec_end (n_rec_periods periods later; the
// step that starts at t_rec_end is not run, event dump_checkpoint stops the run first).
event rec_land_start (t = t_rec0) {}
event rec_land_end   (t = t_rec_end) {}
#endif // RECORD_FLUX

#if REPLAY_FLUX
static long rec_nuse = 0;       // records of the file that belong to the replay
static long rec_lap = 0;        // steps of one pass over the record (one fewer than rec_nuse when looping)
static long rec_k = 0;          // next step of the current pass
static long rec_pos0 = 0;       // file offset of the first step
static int  rec_pass = 0;       // number of completed passes
static int  replay_stop = 0;    // the last step of a non-looping record has been read

static int rec_read_head (int * hi, double * ht, double * hdt, long long * hn)
{
  return fread (hi, sizeof(int), 1, fp_rec) == 1 && fread (ht, sizeof(double), 1, fp_rec) == 1 &&
         fread (hdt, sizeof(double), 1, fp_rec) == 1 && fread (hn, sizeof(long long), 1, fp_rec) == 1;
}

// Called from event init once the checkpoint is restored: opens the record, indexes it and checks that it
// starts where the checkpoint is (the dump header carries iter - 1, see rec_start). It also sets the time of
// the first step, which the restored run does not know (centered.h does this with event("stability")).
// The records of one pass add up to exactly the window length when its last step landed on the window end.
static void replay_open (void)
{
  char name[600];
  rec_file_name (name, sizeof name);
  fp_rec = fopen (name, "rb");
  if (!fp_rec) {
    fprintf (stderr, "ERROR: replay: cannot open %s\n", name);
    exit (1);
  }
  char magic[8];
  int rank = -1, np = -1;
  if (fread (magic, 1, 8, fp_rec) != 8 || memcmp (magic, rec_magic, 8) ||
      fread (&rank, sizeof(int), 1, fp_rec) != 1 || fread (&np, sizeof(int), 1, fp_rec) != 1 ||
      rank != pid() || np != npe()) {
    fprintf (stderr, "ERROR: replay: %s is not a flux record for rank %d of %d\n", name, pid(), npe());
    exit (1);
  }
  rec_pos0 = ftell (fp_rec);
  long long ncell = rec_leaf_count(), hn;
  int hi, i_first = -1;
  double ht, hdt, t_first = 0.;
  long nall = 0;
  rec_nuse = 0;
  while (rec_read_head (&hi, &ht, &hdt, &hn)) {
    if (hn != ncell) {
      fprintf (stderr, "ERROR: replay: record has %lld leaf cells per step, the restored grid has %lld\n",
               hn, ncell);
      exit (1);
    }
    if (nall == 0)
      i_first = hi, t_first = ht;
    if (params.n_rec_periods <= 0 ||
        ht < t_first + params.n_rec_periods*T_per_st - REC_EPS*(1. + fabs (ht)))
      rec_nuse++;
    nall++;
    fseek (fp_rec, (long) (hn*REC_NPC*sizeof(double)), SEEK_CUR);
  }
  rec_lap = rec_nuse;
  if (rec_lap < 1) {
    fprintf (stderr, "ERROR: replay: %ld usable records in %s is too few\n", rec_nuse, name);
    exit (1);
  }
  if (i_first != iter + 1 || fabs (t_first - t) > REC_EPS*(1. + fabs (t))) {
    fprintf (stderr, "ERROR: replay: record starts at i=%d t=%.17g but the checkpoint resumes at i=%d t=%.17g\n",
             i_first, t_first, iter + 1, t);
    exit (1);
  }
  tnext = t_first;
  fseek (fp_rec, rec_pos0, SEEK_SET);
  rec_k = 0;
  if (pid() == 0)
    fprintf (ferr, "replay: %ld records (%ld in file), pass of %ld steps, loop %d, from i=%d t=%.17g\n",
             rec_nuse, nall, rec_lap, params.replay_loop, i_first, t_first);
}

// Scalars on a finer grid than the flow: the record holds uf at the flow level; the faces of the finer
// levels come from Basilisk's own solenoidal face refinement, refine_face_solenoidal() of
// grid/tree-common.h (Popinet 2009, sec. 3.3: linear interpolation plus a local projection that keeps the
// divergence of the children equal to that of the parent), applied level by level. That function is written
// for the moment a cell is refined: it writes the faces shared with a neighbour only when the neighbour is
// not refined yet (is_refined()). The tree is already refined, so the cells of the level being
// prolonged are marked as leaves for the duration of the call, which is the state refine_cell() sees.
static void rec_prolong (void)
{
  for (int l = params.fidelity; l < params.replay_scalar_level; l++) {
    foreach_cell() {
      if (level == l) {
        cell.flags |= leaf;
        continue;
      }
    }
    foreach_level(l)
      refine_face_solenoidal (point, uf.x);
    foreach_cell() {
      if (level == l) {
        cell.flags &= ~leaf;
        continue;
      }
    }
  }
  // The solenoidal refinement knows nothing of the embedded boundary, so it puts flux on faces the wall
  // closes (fs = 0). vof.h divides by fm (un = uf*dt/(Delta*fm + SEPS)), so such a flux gives un ~ 1e30 and
  // the scalars diverge in the first step (measured). The coupled run clears them every step in event adapt
  // of centered.h, and so does this: the same line. It is the only place the prolonged flux is changed
  // after the canonical routine, and it is what makes the divergence at wall cells nonzero (reported).
  foreach_face()
    if (uf.x[] && !fs.x[])
      uf.x[] = 0.;
  // params.replay_project_flux = 1: make the prolonged flux divergence-free on the fine grid with Basilisk's
  // own project() (poisson.h), called as centered.h's event projection calls it (alpha = fm, here with
  // dt = 1: only the divergence-free flux is wanted, and project() scales its tolerance by dt^2). The
  // tolerance is TOLERANCE/sq(dt) on the divergence residual, so it is set to replay_project_tol for the
  // call and put back for the scalar solves. fm = 0 on closed faces, so they stay zero.
  if (params.replay_project_flux) {
    scalar pp[];
    foreach()
      pp[] = 0.;
    double tol0 = TOLERANCE;
    TOLERANCE = params.replay_project_tol > 0. ? params.replay_project_tol : 1e-10;
    mgstats mp = project (uf, pp, fm, 1., 4);
    TOLERANCE = tol0;
    static FILE * fpm = NULL;
    if (!fpm && pid() == 0)
      fpm = fopen ("replay_project.txt", "w");
    if (fpm) {
      fprintf (fpm, "%d %d %g %g\n", iter, mp.i, mp.resb, mp.resa);
      fflush (fpm);
    }
  }
}

// params.replay_diag = 1: per step, into replay_diag.txt, max |div uf|*Delta/max|uf| on the leaf grid, over
// the cells away from the wall (cs = 1 and the four faces open) and over all cells.
static void replay_diag (void)
{
  double dfar = 0., dall = 0., umax = 0.;
  foreach (reduction(max:dfar) reduction(max:dall) reduction(max:umax)) {
    double d = 0.;
    foreach_dimension() {
      d += uf.x[1] - uf.x[];
      if (fabs (uf.x[]) > umax)
        umax = fabs (uf.x[]);
    }
    d = fabs (d);
    if (d > dall)
      dall = d;
    if (cs[] == 1. && fs.x[] == 1. && fs.x[1] == 1. && fs.y[] == 1. && fs.y[0,1] == 1. && d > dfar)
      dfar = d;
  }
  static FILE * fp = NULL;
  if (!fp && pid() == 0)
    fp = fopen ("replay_diag.txt", "w");
  if (fp) {
    fprintf (fp, "%d %.17g %.17g %.17g %.17g\n", iter, t, dfar, dall, umax);
    fflush (fp);
  }
}

// Called from event vof at the start of every step: takes the place of stability + the Navier-Stokes
// solve. Sets dt and the time of the next step, and loads the flux that step advects with.
// With scalars on a grid finer than the flow by k levels, the recorded dt violates the CFL limit of the VOF
// step (it was set for the flow grid: measured, CFL 0.5 -> 0.91 at one level, and the scalar solves blow
// up). Each recorded step is therefore taken as n = 2^k substeps of dt/n with the same flux, the flux being
// constant over a recorded step in the coupled run too. Substeps are iterations, so `i` advances n times
// per record and the i check is dropped (the time check stays).
static int rec_sub = 0;         // substeps of the current record already taken
static double rec_hdt = 0.;     // dt of the current record

static void replay_load_step (void)
{
  int nsub = 1;
  if (params.replay_scalar_level > params.fidelity)
    nsub = 1 << (params.replay_scalar_level - params.fidelity);
  if (rec_sub == 0) {
    if (rec_k >= rec_lap) {     // end of the pass: time goes on, the flow restarts from its first step
      fseek (fp_rec, rec_pos0, SEEK_SET);
      rec_k = 0;
      rec_pass++;
    }
    int hi;
    double ht;
    long long hn;
    if (!rec_read_head (&hi, &ht, &rec_hdt, &hn)) {
      fprintf (stderr, "ERROR: replay: record ends early at step %ld of the pass\n", rec_k);
      exit (1);
    }
    if (rec_pass == 0 && ((nsub == 1 && hi != iter) || fabs (ht - t) > REC_EPS*(1. + fabs (t)))) {
      fprintf (stderr, "ERROR: replay: out of step with the record: i=%d t=%.17g, record has i=%d t=%.17g\n",
               iter, t, hi, ht);
      exit (1);
    }
    double * buf = (double *) malloc (hn*REC_NPC*sizeof(double));
    if (fread (buf, sizeof(double), hn*REC_NPC, fp_rec) != (size_t) (hn*REC_NPC)) {
      fprintf (stderr, "ERROR: replay: truncated record at step %ld\n", rec_k);
      exit (1);
    }
    // params.replay_hold = n > 0 emulates an offline solver that stores n snapshots per period: the flux is
    // loaded at the first step of each snapshot interval (time from the first record) and held over the rest.
    static double rec_t_first = 0.;
    static long rec_snap = -1;
    if (rec_pass == 0 && rec_k == 0)
      rec_t_first = ht;
    long snap = params.replay_hold > 0 ?
      (long) floor ((ht - rec_t_first)/(T_per_st/params.replay_hold) + REC_EPS) : -2;
    if (snap != rec_snap || params.replay_hold <= 0)
      rec_scatter (buf);
    rec_snap = snap;
    free (buf);
    if (params.replay_scalar_level > params.fidelity)
      rec_prolong();
    if (params.replay_diag)
      replay_diag();
  }
#if REC_CHECKSUM
  rec_checksum();
#endif
#if REC_CHECKSUM
  rec_where_max();
#endif
  dt = rec_hdt/nsub;
  tnext = t + dt;               // dt already lands on every event time of the recorded run
  if (++rec_sub == nsub) {
    rec_sub = 0;
    rec_k++;
    if (!params.replay_loop && rec_k == rec_lap)
      replay_stop = 1;
  }
}

// The Basilisk time events `t += dt_out` are dropped when t steps over their time, and a looped replay
// does not land on the lattice of output times. So replay outputs at the first step whose time is at or
// after the next lattice point. In the first pass the steps do land on it, so the times agree with the
// recorded run.
static int out_due (double * t_next, double dt_out)
{
  if (*t_next < 0.)
    *t_next = dt_out*ceil (t/dt_out - REC_EPS);
  if (t < *t_next - REC_EPS*(1. + fabs (t)))
    return 0;
  *t_next = dt_out*(floor (t/dt_out + REC_EPS) + 1.);
  return 1;
}
#endif // REPLAY_FLUX
