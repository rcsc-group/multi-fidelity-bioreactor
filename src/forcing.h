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
// Angle Th(t), its rate Th_d and acceleration Th_2d from forcing_now() (forcing_state.h): the forcing ramps
// from its source (rest on a fresh start; the dump's label on a restart) to params.json over N_RAMP_CYCLES,
// in amplitude, phase and frequency. Th_d and Th_2d are the exact time derivatives, ramp terms included
// (pitfall "ramp-derivatives"). Optional horizontal shaking (omega_h > 0), off by default.
// ---------------------------------------------------------------------------------------------------
event acceleration(i++)
{
  ForcingNow F = forcing_now (t);
  _forcing_last = F;  _forcing_last_t = t;   // for the checkpoint label (forcing_state.h)
  Th = F.Th;  Th_d = F.Th_d;  Th_2d = F.Th_2d;
  if (_restart_kind > 0 && !_restart_forcing_written) {   // first step after a restart, t = t_checkpoint
    write_restart_forcing (F, t, _restart_kind);          // (not i == 0: restore() brings back i)
    _restart_forcing_written = 1;
  }
  _theta_env_deg = F.theta_max[0];
  double x_acc = F.x_acc;   // horizontal shaking

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
