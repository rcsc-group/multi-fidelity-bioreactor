// snapshots.h -- fields on a uniform n x n grid (n = 2^level), as raw float32 files for Python.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

#if !REPLAY_FLUX
// fields/snap_NNNNNN.bin: written on the video cadence between params.snap_start_cycle and snap_end_cycle
// (off when snap_end_cycle = 0). The analysis picks the frames it needs (a phase, the instant chi = 0.5,
// ...) afterwards, so no threshold logic is needed in the solver.
// Format: int n, int nplanes, double t, then nplanes float planes of n*n values, in this order:
//   ux, uy, omega, f, cs, c (line), c1 (left half), c2 (top half), c3 (circle), oxy.
// Readers must take the plane count from the header.
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

// uv_field.bin: velocity at the checkpoint instant (a period boundary), for grid-convergence norms across
// levels (Kim et al., Fig. A.16c). Format: int n, double t, then float planes ux, uy, f, cs of n*n values.
// interpolate() is collective: every rank calls it for every point.
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
