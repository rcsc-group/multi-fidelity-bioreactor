// sample_grid.h -- a scalar on the uniform n x n grid of cell centres (X0 + (i + 0.5) L0/n), as float32.
// Part of BioReactor.c: included before the writers that use it (forcing.h, snapshots.h).
//
// Same values as interpolate (s, x, y) point by point (bilinear, interpolate_linear(), cartesian-common.h),
// at the same points, but with ONE MPI reduction per scalar instead of one per point: each rank fills the
// points it owns (foreach_point, as in output_field(), output.h), the others stay nodata, and the min over
// ranks picks the owner's value, as interpolate()'s reduction (min:val) does. Profiling 2026-11-23 (diary):
// point-wise interpolate() made 94 % of all MPI reductions of an L8 run (4 x NN^2 per period end).

static void sample_grid (scalar s, int n, float * out)
{
  double dx = L0 / n;
  size_t nn = (size_t) n * n;
  double * v = (double *) malloc (nn * sizeof (double));
  for (size_t k = 0; k < nn; k++)
    v[k] = nodata;
  for (int j = 0; j < n; j++) {
    double yj = Y0 + (j + 0.5) * dx;
    for (int i = 0; i < n; i++) {
      double xi = X0 + (i + 0.5) * dx;
      size_t k = (size_t) j * n + i;
      foreach_point (xi, yj) {
        double * alias = v; // so that qcc considers 'v' a local variable (as in output_field())
        alias[k] = interpolate_linear (point, s, xi, yj, 0.);
      }
    }
  }
#if _MPI
  mpi_all_reduce_array (v, MPI_DOUBLE, MPI_MIN, (int) nn);
#endif
  for (size_t k = 0; k < nn; k++)
    out[k] = (float) v[k];
  free (v);
}
