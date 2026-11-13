// movies.h -- Video frames (-DVIDEOS): events movies_output, movies_output_tau.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

// Make videos — field data exported as binary frames; rendered by scripts/render_videos.py
#if VIDEOS
int _vframe = 0;
double _last_vdump = -1e30;

// [PROJECT CHANGED, 2026-07-30] Basilisk's `t += dt_video` event-repeat
// syntax requires dt_video to be a compile-time constant (qcc error:
// "events can only use a single condition" once dt_video became a
// runtime value driven by params.frames_per_period) -- switched to a
// plain per-timestep check so the cadence can be runtime-configurable.
event movies_output(i++)
{
  if (t < t_mix || t > t_end || t - _last_vdump < dt_video)
    return;
  _last_vdump = t;

  // MPI-safe video output:
  // interpolate() is collective (MPI_Allreduce per point) — ALL ranks must call it.
  // File I/O is guarded to rank 0 only so no two ranks open the same file.
  // _vframe increments on all ranks so the counter stays in sync.
#if _MPI
  if (pid() == 0)
#endif
  {
    if (_vframe == 0)
      system("mkdir -p frames");
  }

  int    n    = NN;  // NN = 1<<fidelity, set in main()
  double t_nd = t;
  double dx   = L0 / n;

  // Horizontal displacement (lab frame): X_lab = sum_k A_k * sin(k*w_h*t + phi_k)
  double xh_nd = 0.0;
  if (params.omega_h > 0.) {
    double w_h_nd = params.omega_h * T_bio;
    for (int k = 1; k <= params.n_harmonics; k++) {
      double wk = k * w_h_nd;
      xh_nd += (params.amplitude_h[k-1] / L_bio) * sin(wk*t + params.phi_horizontal[k-1]);
    }
  }

  // VOF field f interpolated onto n×n uniform grid, row-major (j=0 → y=Y0 = bottom).
  // Run on ALL ranks — interpolate() uses MPI_Allreduce internally and requires
  // every rank to participate for each (xi, yj) query.
  float *buf = (float *)malloc(n * n * sizeof(float));
  int idx = 0;
  for (int j = 0; j < n; j++) {
    double yj = Y0 + (j + 0.5) * dx;
    for (int i = 0; i < n; i++) {
      double xi = X0 + (i + 0.5) * dx;
      buf[idx++] = (float)interpolate(f, xi, yj);
    }
  }

  // Only rank 0 writes; all ranks already have the correct buf values.
#if _MPI
  if (pid() == 0)
#endif
  {
    char fpath[512];
    sprintf(fpath, "frames/frame_%06d.bin", _vframe);
    FILE *fp = fopen(fpath, "wb");
    if (fp) {
      // [PROJECT ADDED, 2026-09-07] w_bio (dimensional rad/s)/_theta_env_deg
      // appended after xh_nd for the checkpoint-chaining visualization
      // (diary.md 2026-09-07) -- see render_videos.py's load_frame() for
      // the matching read order. w_bio is set once per segment in main()
      // and never interpolated across a restart (unlike the amplitude
      // envelope) -- deliberately dumped raw so a frequency-changing
      // checkpoint transition shows its true instantaneous jump, not a
      // smoothed value.
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

// [PROJECT ADDED, 2026-08-05] Shear-stress field video export, requested to
// visualize where tau_100_max (the pointwise-max KPI) occurs and how it
// moves over time. No upstream equivalent (Kim et al.'s driver computes no
// shear stress at all -- see the CAUTION above event normcal).
//
// Separate from movies_output (VOF field) rather than folding into it,
// because the gating condition must differ: movies_output only starts at
// t_mix (oxygen/tracer injection), but shear stress is present from t=0 and
// is unrelated to injection timing -- gating this on t_mix would produce zero
// frames for any run with n_mix_cycles large enough that t_mix > t_end (as
// happened with l10_kim_seg2, the run this was written to visualize).
//
// tau_field[] is declared here rather than reusing a scalar declared inside
// an event body: the CAUTION comment above event normcal already documents
// that declaring a new scalar (there, tau_liq[]) inside an event leaks a
// Basilisk scalar on every call and corrupts the global scalar list,
// segfaulting at fidelity >= 7. Declaring it once at file scope (like c[],
// oxy[] above) avoids that hazard entirely -- it's allocated once, not
// once per event call.
#if VIDEOS
scalar tau_field[];
// [PROJECT ADDED, 2026-08-08] ediss_field: same rationale as tau_field --
// declared once at file scope, not inside the event, to avoid the
// documented scalar-leak/segfault hazard. Needed for Fig. 8b/c
// reproduction (normalized histograms of tau and epsilon across the
// liquid, at the instant each domain-mean peaks) -- shear_stress.dat only
// carries the scalar mean, not the full spatial distribution.
scalar ediss_field[];
int _vframe_tau = 0;
double _last_vdump_tau = -1e30;

event movies_output_tau(i++)
{
  if (t > t_end || t - _last_vdump_tau < dt_video)
    return;
  _last_vdump_tau = t;

#if _MPI
  if (pid() == 0)
#endif
  {
    if (_vframe_tau == 0)
      system("mkdir -p frames_tau");
  }

  // Same stencils as event normcal's tau/ediss computation, masked to
  // liquid (f[] > 0.5) -- kept identical on purpose so the video fields
  // match the tau_100_max/tau_mean_max/ediss_mean KPIs in shear_stress.dat.
  //
  // [PROJECT CHANGED, 2026-08-09] tau_field now stores the SIGNED shear
  // stress, not fabs() -- Kim et al.'s Fig. 8b histogram bins the signed
  // quantity (their Fig_tau_Ediss.pdf shows it centered near zero and
  // ranging -15e-3 to +15e-3 Pa), and a fabs()'d field can't be un-signed
  // in post -- the reverse (taking abs() of a signed field at render/
  // analysis time, if a magnitude heatmap is wanted) is always possible.
  // NOTE: this changes the on-disk semantics of frames_tau/*.bin's second
  // buffer for any run recorded after this commit; older recordings
  // (l10_kim_tau_video, l10_kim_fig8, etc.) still store fabs()'d values.
  foreach() {
    // same bag mask as the averaging pass (2026-09-15 fix): outside the
    // embedded bag there is no fluid, so these fields must be zero there.
    if (cs[] > 0. && f[] > 0.5) {
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

  // Three buffers per frame: f (for the bag/liquid mask), tau_field, and
  // ediss_field (2026-08-08 addition, for Fig. 8b/c histograms). All
  // interpolated on ALL ranks (interpolate() is collective); only rank 0
  // writes. NOTE: this is a FORMAT CHANGE from the original 2-buffer
  // (f, tau) frames_tau/*.bin -- older recordings (e.g. l10_kim_tau_video,
  // l10_kim_strict_mask_test, l10_kim_signed_test) are still 2-buffer and
  // must be read with the old loader; anything recorded after this commit
  // is 3-buffer.
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


