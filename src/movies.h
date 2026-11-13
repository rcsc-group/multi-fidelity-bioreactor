// movies.h -- video frames (-DVIDEOS=1): raw float32 fields on a uniform grid, rendered by
// scripts/render_videos.py. Part of BioReactor.c: included at one fixed place; the place sets the event order.
//
// Both events check the frame time every step (`t += dt_video` needs a compile-time constant, and
// dt_video comes from params). interpolate() is collective, so every rank interpolates; rank 0 writes.

#if VIDEOS
int _vframe = 0;
double _last_vdump = -1e30;

// frames/frame_NNNNNN.bin, from t_mix to t_end: the volume fraction f.
// Format: int n, double t, double Th (angle), double xh (horizontal offset), double w_bio (rad/s),
// double theta envelope (deg), then n*n floats of f, row-major from the bottom (j = 0 is y = Y0).
// render_videos.py load_frame() reads the same order.
event movies_output(i++)
{
  if (t < t_mix || t > t_end || t - _last_vdump < dt_video)
    return 0;   // never a bare `return;` in an event (pitfall "bare-return")
  _last_vdump = t;

#if _MPI
  if (pid() == 0)
#endif
  {
    if (_vframe == 0)
      system("mkdir -p frames");
  }

  int    n    = NN;
  double t_nd = t;
  double dx   = L0 / n;

  // Horizontal displacement (lab frame): X = sum_k A_k sin(k w_h t + phi_k)
  double xh_nd = 0.0;
  if (params.omega_h > 0.) {
    double w_h_nd = params.omega_h * T_bio;
    for (int k = 1; k <= params.n_harmonics; k++) {
      double wk = k * w_h_nd;
      xh_nd += (params.amplitude_h[k-1] / L_bio) * sin(wk*t + params.phi_horizontal[k-1]);
    }
  }

  float *buf = (float *)malloc(n * n * sizeof(float));
  int idx = 0;
  for (int j = 0; j < n; j++) {
    double yj = Y0 + (j + 0.5) * dx;
    for (int i = 0; i < n; i++) {
      double xi = X0 + (i + 0.5) * dx;
      buf[idx++] = (float)interpolate(f, xi, yj);
    }
  }

#if _MPI
  if (pid() == 0)
#endif
  {
    char fpath[512];
    sprintf(fpath, "frames/frame_%06d.bin", _vframe);
    FILE *fp = fopen(fpath, "wb");
    if (fp) {
      fwrite(&n,              sizeof(int),    1, fp);
      fwrite(&t_nd,           sizeof(double), 1, fp);
      fwrite(&Th,             sizeof(double), 1, fp);
      fwrite(&xh_nd,          sizeof(double), 1, fp);
      fwrite(&w_bio,          sizeof(double), 1, fp);
      fwrite(&_theta_env_deg, sizeof(double), 1, fp);
      fwrite(buf,    sizeof(float), n * n, fp);
      fclose(fp);
    } else {
      fprintf(stderr, "movies_output: cannot open %s\n", fpath);
    }
  }

  free(buf);
  _vframe++;
}
#endif

// frames_tau/frame_NNNNNN.bin, from t = 0 (shear exists before the release) to t_end: f, the SIGNED shear
// stress and the dissipation rate, same stencil and mask as event normcal.
// Format: int n, double t, double Th, double xh, then three n*n float planes: f, tau, ediss.
// (Older recordings: two planes, or |tau|; see docs_site/reference/output-files.md.)
#if VIDEOS
// Declared once here, not inside the event (a field declared in an event is allocated on every call).
scalar tau_field[];
scalar ediss_field[];
int _vframe_tau = 0;
double _last_vdump_tau = -1e30;

event movies_output_tau(i++)
{
  if (t > t_end || t - _last_vdump_tau < dt_video)
    return 0;   // never a bare `return;` in an event (pitfall "bare-return")
  _last_vdump_tau = t;

#if _MPI
  if (pid() == 0)
#endif
  {
    if (_vframe_tau == 0)
      system("mkdir -p frames_tau");
  }

  foreach() {
    if (cs[] > 0. && f[] > 0.5) {   // liquid inside the bag; 0 elsewhere
      double du_dy = (u.x[0,1] - u.x[0,-1]) / (2.*Delta);
      double dv_dx = (u.y[1]   - u.y[-1])   / (2.*Delta);
      double du_dx = (u.x[1]   - u.x[-1])   / (2.*Delta);
      double dv_dy = (u.y[0,1] - u.y[0,-1]) / (2.*Delta);
      tau_field[]   = mu(f[]) * (du_dy + dv_dx);
      ediss_field[] = mu(f[]) * (2.*du_dx*du_dx + 2.*dv_dy*dv_dy + (du_dy+dv_dx)*(du_dy+dv_dx));
    } else {
      tau_field[]   = 0.;
      ediss_field[] = 0.;
    }
  }

  int    n    = NN;
  double t_nd = t;
  double dx   = L0 / n;

  double xh_nd = 0.0;
  if (params.omega_h > 0.) {
    double w_h_nd = params.omega_h * T_bio;
    for (int k = 1; k <= params.n_harmonics; k++) {
      double wk = k * w_h_nd;
      xh_nd += (params.amplitude_h[k-1] / L_bio) * sin(wk*t + params.phi_horizontal[k-1]);
    }
  }

  float *buf_f     = (float *)malloc(n * n * sizeof(float));
  float *buf_tau   = (float *)malloc(n * n * sizeof(float));
  float *buf_ediss = (float *)malloc(n * n * sizeof(float));
  int idx = 0;
  for (int j = 0; j < n; j++) {
    double yj = Y0 + (j + 0.5) * dx;
    for (int i = 0; i < n; i++) {
      double xi = X0 + (i + 0.5) * dx;
      buf_f[idx]     = (float)interpolate(f, xi, yj);
      buf_tau[idx]   = (float)interpolate(tau_field, xi, yj);
      buf_ediss[idx] = (float)interpolate(ediss_field, xi, yj);
      idx++;
    }
  }

#if _MPI
  if (pid() == 0)
#endif
  {
    char fpath[512];
    sprintf(fpath, "frames_tau/frame_%06d.bin", _vframe_tau);
    FILE *fp = fopen(fpath, "wb");
    if (fp) {
      fwrite(&n,         sizeof(int),    1, fp);
      fwrite(&t_nd,      sizeof(double), 1, fp);
      fwrite(&Th,        sizeof(double), 1, fp);
      fwrite(&xh_nd,     sizeof(double), 1, fp);
      fwrite(buf_f,      sizeof(float), n * n, fp);
      fwrite(buf_tau,    sizeof(float), n * n, fp);
      fwrite(buf_ediss,  sizeof(float), n * n, fp);
      fclose(fp);
    } else {
      fprintf(stderr, "movies_output_tau: cannot open %s\n", fpath);
    }
  }

  free(buf_f);
  free(buf_tau);
  free(buf_ediss);
  _vframe_tau++;
}
#endif
