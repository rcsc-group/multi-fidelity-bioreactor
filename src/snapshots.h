// snapshots.h -- Uniform-grid field snapshots and the uv field written with the checkpoint.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

#if !REPLAY_FLUX
// [PROJECT ADDED, 2026-09-16] Velocity field on a uniform grid, written once,
// at the same instant as the checkpoint. Kim's Fig. A.16(c) is the discrete
// L2 norm of the velocity error between each resolution and a finer
// reference, which needs the FIELDS at one common instant across levels --
// and nothing this code wrote was per-cell velocity: normf.dat carries rms
// scalars, frames_tau carries (f, tau, eps), and the Basilisk dump stores a
// tree the analysis scripts cannot read. Written unconditionally rather than
// under VIDEOS because it is one file per run (8 MB at L10) and a run whose
// velocity field cannot be compared to anything is a run that has to be
// repeated to answer a convergence question.
//
// The instant is the checkpoint's own, i.e. a period boundary (theta=0), not
// Kim's t/T_p = 29.77 (maximum rocking angle): t_dump_checkpoint is pinned to
// a zero-crossing by construction. For a convergence measure the phase only
// has to be the SAME across levels, which this guarantees and an arbitrary
// mid-cycle time would not.
// [PROJECT ADDED, 2026-09-16] Field snapshot on a uniform grid: velocity,
// vorticity, the interface, the embedded geometry, all four tracer
// configurations and oxygen. Kim's Figs. 3, 6 and 7 are FIELDS at instants
// that are not known in advance -- four phases of one rocking cycle, the
// instant chi reaches 0.5, the instants C* reaches 0.10/0.25/0.50 -- so the
// solver records a window on the video cadence and the analysis picks the
// frame it wants from the time series it already has. No threshold logic
// enters the solver, where getting it wrong costs a rerun.
//
// Nine planes, float32, in this order:
//   ux, uy, omega, f, cs, c (line), c1 (left half), c2 (top half), c3 (circle)
// followed by oxy. Readers must take the plane count from the header.
#define SNAP_NPLANES 10

static void write_snapshot (const char *path, double t_nd)
{
  int    n  = NN;
  double dx = L0 / n;
  scalar omega[];
  vorticity (u, omega);

  scalar *planes = {u.x, u.y, omega, f, cs, c, c1, c2, c3, oxy};
  int np = SNAP_NPLANES;
  float *buf = (float *) malloc ((size_t) np * n * n * sizeof (float));
  int k = 0;
  for (scalar s in planes) {
    for (int j = 0; j < n; j++) {
      double yj = Y0 + (j + 0.5) * dx;
      for (int i = 0; i < n; i++) {
        double xi = X0 + (i + 0.5) * dx;
        // interpolate() is collective: every rank calls it for every cell.
        buf[(size_t) k * n * n + j * n + i] = (float) interpolate (s, xi, yj);
      }
    }
    k++;
  }
#if _MPI
  if (pid() == 0)
#endif
  {
    FILE *fp = fopen (path, "wb");
    if (fp) {
      fwrite (&n,    sizeof(int),    1, fp);
      fwrite (&np,   sizeof(int),    1, fp);
      fwrite (&t_nd, sizeof(double), 1, fp);
      fwrite (buf,   sizeof(float), (size_t) np * n * n, fp);
      fclose (fp);
    } else
      fprintf (stderr, "write_snapshot: cannot open %s\n", path);
  }
  free (buf);
}

static int    _snap_idx  = 0;
static double _last_snap = -1e30;

event field_snapshots (i++)
{
  if (params.snap_end_cycle > 0.) {
    double c0 = params.snap_start_cycle * T_per_st;
    double c1_ = params.snap_end_cycle   * T_per_st;
    if (t >= c0 && t <= c1_ && t - _last_snap >= dt_video) {
      _last_snap = t;
#if _MPI
      if (pid() == 0)
#endif
      { if (_snap_idx == 0) system ("mkdir -p fields"); }
      char path[512];
      sprintf (path, "fields/snap_%06d.bin", _snap_idx++);
      write_snapshot (path, t);
    }
  }
}

static void write_uv_field (double t_nd)
{
  int    n  = NN;
  double dx = L0 / n;
  float *bux = (float *) malloc (n * n * sizeof (float));
  float *buy = (float *) malloc (n * n * sizeof (float));
  float *bf  = (float *) malloc (n * n * sizeof (float));
  float *bcs = (float *) malloc (n * n * sizeof (float));
  int idx = 0;
  for (int j = 0; j < n; j++) {
    double yj = Y0 + (j + 0.5) * dx;
    for (int i = 0; i < n; i++) {
      double xi = X0 + (i + 0.5) * dx;
      // interpolate() is collective -- every rank must call it for all cells.
      bux[idx] = (float) interpolate (u.x, xi, yj);
      buy[idx] = (float) interpolate (u.y, xi, yj);
      bf [idx] = (float) interpolate (f,   xi, yj);
      bcs[idx] = (float) interpolate (cs,  xi, yj);
      idx++;
    }
  }
#if _MPI
  if (pid() == 0)
#endif
  {
    FILE *fp = fopen ("uv_field.bin", "wb");
    if (fp) {
      fwrite (&n,    sizeof(int),    1, fp);
      fwrite (&t_nd, sizeof(double), 1, fp);
      fwrite (bux,   sizeof(float), n * n, fp);
      fwrite (buy,   sizeof(float), n * n, fp);
      fwrite (bf,    sizeof(float), n * n, fp);
      fwrite (bcs,   sizeof(float), n * n, fp);
      fclose (fp);
    } else
      fprintf (stderr, "write_uv_field: cannot open uv_field.bin\n");
  }
  free (bux); free (buy); free (bf); free (bcs);
}

#endif // !REPLAY_FLUX
