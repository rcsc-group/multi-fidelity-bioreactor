// forcing.h -- the rocking motion (event acceleration), steady-streaming vorticity, droplet removal.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.
//
// The bag does not move in the computation. We solve in the frame of the bag, so the rocking appears as
// body forces: gravity at the tilt angle Th(t), plus the Coriolis, centrifugal and Euler (angular
// acceleration) forces about the rocking axis, which is L_piv below the centre.

#if ACCELERATION
#if STREAMING_VORTICITY
// ---------------------------------------------------------------------------------------------------
// Steady-streaming vorticity (Kim et al.): average the WATER velocity over one rocking period, take the
// curl of that average, then average |curl| over the water. streaming.dat: one row per period.
// streaming_field.bin: the last period's mean-flow fields on a uniform grid (for streamline plots).
// ---------------------------------------------------------------------------------------------------
event streaming_vorticity (i++)
{
  // Liquid-weighted time average ubar = INT(u f dt) / INT(f dt). A plain INT(u dt)/T mixes water and air
  // in the cells the interface sweeps through and gives a large false mean shear there. (Check: the
  // mean-flow vorticity cannot exceed the time average of |vorticity|.) u*dt accumulates exactly for any
  // time step; the window restarts each period.
  foreach() {
    u_acc.x[] += u.x[]*f[]*dt;
    u_acc.y[] += u.y[]*f[]*dt;
    f_acc[]   += f[]*dt;
  }
  streaming_acc_dt += dt;

  if (streaming_acc_dt >= T_per_st) {
    vector ubar[];
    scalar omega_bar[];
    foreach() {
      double w = f_acc[];
      ubar.x[] = w > 0. ? u_acc.x[]/w : 0.;
      ubar.y[] = w > 0. ? u_acc.y[]/w : 0.;
    }
    boundary ((scalar *){ubar});
    vorticity (ubar, omega_bar);

    // Average over cells that were water for (almost) the whole period and are not cut by the wall
    // (cs = 1, as in Kim et al.'s postprocessing): cut cells next to the wall have very large gradients.
    double vsum = 0., varea = 0., vmax = 0.;
    foreach (reduction(+:vsum) reduction(+:varea) reduction(max:vmax)) {
      if (cs[] > 1. - 1e-10 && f_acc[] > 0.99*streaming_acc_dt) {
        double dA = cs[]*(Delta*Delta);
        double a  = fabs(omega_bar[]);
        vsum  += a*dA;
        varea += dA;
        if (a > vmax) vmax = a;
      }
    }
    if (pid() == 0 && fp_streaming) {
      fprintf (fp_streaming, "%i %g %g %g %g %g \n", i, t,
               varea > 0. ? vsum/varea : 0., vmax, varea, streaming_acc_dt);
      fflush (fp_streaming);
    }

    // streaming_field.bin (overwritten every period; the last one is the most settled):
    // int n, int nplanes, double t, double window, then float planes ubar.x, ubar.y, omega_bar, f_acc on
    // an n x n grid. f_acc / window = fraction of the period the cell was water.
    {
      int    n_s  = NN;
      double dx_s = L0 / n_s;
      scalar *sp = {ubar.x, ubar.y, omega_bar, f_acc};
      int np_s = 4;
      float *bs = (float *) malloc ((size_t) np_s * n_s * n_s * sizeof (float));
      int ks = 0;
      for (scalar s in sp) {
        for (int j = 0; j < n_s; j++) {
          double yj = Y0 + (j + 0.5) * dx_s;
          for (int ii = 0; ii < n_s; ii++) {
            double xi = X0 + (ii + 0.5) * dx_s;
            bs[(size_t) ks * n_s * n_s + j * n_s + ii] =
              (float) interpolate (s, xi, yj);
          }
        }
        ks++;
      }
#if _MPI
      if (pid() == 0)
#endif
      {
        FILE *fs_ = fopen ("streaming_field.bin", "wb");
        if (fs_) {
          fwrite (&n_s,  sizeof(int),    1, fs_);
          fwrite (&np_s, sizeof(int),    1, fs_);
          fwrite (&t,    sizeof(double), 1, fs_);
          fwrite (&streaming_acc_dt, sizeof(double), 1, fs_);
          fwrite (bs, sizeof(float), (size_t) np_s * n_s * n_s, fs_);
          fclose (fs_);
        }
      }
      free (bs);
    }
    foreach() {
      u_acc.x[] = 0.;
      u_acc.y[] = 0.;
      f_acc[]   = 0.;
    }
    streaming_acc_dt = 0.;
  }
}
#endif

// ---------------------------------------------------------------------------------------------------
// event acceleration: body force of the rocking motion, every step.
//
// Angle: Th(t) = sum_k Ak sin(k w (t + t_phase_offset) + phk), k = 1..n_harmonics (1 = Kim et al.).
// Ramp: amplitude Ak and phase phk go from their *_prev values (0 on a fresh start) to the target values
// with a smooth step alpha = 3x^2 - 2x^3, x = (t - t_ramp_start)/(N_RAMP_CYCLES T). Th_d and Th_2d are
// the exact first and second time derivatives, including the alpha' and alpha'' terms of the ramp. (Leaving
// those out makes the forcing inconsistent with the angle during the ramp, and a restart then settles in
// a different state.) The frequency is not ramped.
// Optional horizontal shaking (omega_h > 0), off by default.
// ---------------------------------------------------------------------------------------------------
event acceleration(i++)
{
  double elapsed  = t - t_ramp_start;
  double ramp_dur = N_RAMP_CYCLES * T_per_st;
  int    in_ramp  = elapsed < ramp_dur;
  double x_ss     = in_ramp ? elapsed / ramp_dur : 1.0;
  double alpha    = 3.*x_ss*x_ss - 2.*x_ss*x_ss*x_ss;   // 0 -> 1
  double dalpha_dt   = in_ramp ? 6.*x_ss*(1.-x_ss)/ramp_dur : 0.0;
  double d2alpha_dt2 = in_ramp ? 6.*(1.-2.*x_ss)/(ramp_dur*ramp_dur) : 0.0;

  Th = 0;  Th_d = 0;  Th_2d = 0;
  for (int k = 1; k <= params.n_harmonics; k++) {
    double wk   = k * w_bio_st;
    double dAk  = params.theta_max[k-1]   - params.theta_max_prev[k-1];    // deg
    double dphk = params.phi_angular[k-1] - params.phi_angular_prev[k-1];  // rad
    double Ak   = (params.theta_max_prev[k-1]   + alpha*dAk)  * pi / 180.;
    double phk  =  params.phi_angular_prev[k-1] + alpha*dphk;
    double Ak_d   = dalpha_dt   * dAk * pi / 180.;
    double Ak_2d  = d2alpha_dt2 * dAk * pi / 180.;
    double phk_d  = dalpha_dt   * dphk;
    double phk_2d = d2alpha_dt2 * dphk;
    if (k == 1)
      _theta_env_deg = Ak * 180. / pi;

    double theta_arg = wk*(t + t_phase_offset) + phk;
    double theta_d    = wk + phk_d;
    double s = sin(theta_arg), c = cos(theta_arg);

    Th    += Ak * s;
    Th_d  += Ak_d*s + Ak*theta_d*c;
    Th_2d += Ak_2d*s + 2.*Ak_d*theta_d*c + Ak*phk_2d*c - Ak*theta_d*theta_d*s;
  }

  double x_acc = 0.;   // horizontal shaking
  {
    double omega_h_now = (1.-alpha)*params.omega_h_prev + alpha*params.omega_h;
    if (omega_h_now > 0.) {
      double w_h_st = omega_h_now * T_bio;
      for (int k = 1; k <= params.n_harmonics; k++) {
        double wk_h = k * w_h_st;
        double Ah   = (1.-alpha)*params.amplitude_h_prev[k-1]
                    +     alpha *params.amplitude_h[k-1];
        double phh  = (1.-alpha)*params.phi_horizontal_prev[k-1]
                    +     alpha *params.phi_horizontal[k-1];
        x_acc += (Ah / L_bio) * sq(wk_h) * sin(wk_h*(t + t_phase_offset) + phh);
      }
    }
  }

  // a = gravity + Coriolis + centrifugal + Euler (+ horizontal shaking, x only)
  face vector av = a;
  foreach_face(x)
    av.x[] = -sin(Th)/(Fr*Fr) + 2*Th_d*(u.y[] + u.y[-1,0])*0.5
    + Th_d*Th_d*(x+L_piv*sin(Th)) + Th_2d*(y+L_piv*cos(Th)) + x_acc;
  foreach_face(y)
    av.y[] = -cos(Th)/(Fr*Fr) - 2*Th_d*(u.x[] + u.x[0,-1])*0.5
    + Th_d*Th_d*(y+L_piv*cos(Th)) - Th_2d*(x+L_piv*sin(Th));
  a = av;
}
#endif

// Remove small droplets and bubbles (params.remove_drop = 1; Kim et al.'s REMOVE_DROP option).
// An `if`, not a bare `return`: a bare return in an event stops the run (pitfall "bare-return").
event remove_drop(i++) {
  if (params.remove_drop) {
    remove_droplets(f, remove_minsize, remove_threshold, false);  // droplets
    remove_droplets(f, remove_minsize, remove_threshold, true);   // bubbles
  }
}
