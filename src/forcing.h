// forcing.h -- Rocking motion: events acceleration, streaming_vorticity, remove_drop, velocity_kick.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.


// ================================================================== //
//                           ACCELERATION                             //
// ================================================================== //
#if ACCELERATION
#if STREAMING_VORTICITY
// Kim's steady-streaming vorticity: average u over ONE rocking period, then
// take the curl of that average, then spatially average |curl| over the WATER.
//
// Accumulating u*dt every step (rather than sampling) makes the window average
// exact for any variable timestep. The window is reset each period, so a
// restart loses at most one window rather than corrupting the series.
//
// Masked with `cs[] > 0. && f[] > 0.5` and area-weighted by cs[]*Delta*Delta,
// the same bag mask the tau/EDR diagnostics use -- `f` alone fills the whole
// lower half-DOMAIN including cells outside the embedded bag (diary.md
// 2026-09-15 (2)). Normalising by that liquid area is what makes this a water
// average, as Kim defines it, rather than a whole-fluid-volume average.
event streaming_vorticity (i++)
{
  // LIQUID-CONDITIONED time average: accumulate u*f dt and f dt, so
  // ubar = INT(u f dt)/INT(f dt) is the mean velocity OF THE WATER.
  //
  // [PROJECT FIXED, 2026-09-15] A plain INT(u dt)/T average is contaminated by
  // the free surface: the interface sloshes through a band of cells, so those
  // cells hold water for part of the period and air for the rest, and the
  // unconditioned average blends the two into an enormous spurious mean shear.
  // Measured (job 6427163, L6): it gave <|xi_bar|> = 10.9 non-dim with a local
  // max of 241, ABOVE the time-averaged |instantaneous vorticity| (~12.4 once
  // correctly normalised by liquid area) -- impossible for a genuine mean-flow
  // vorticity, since omega is linear in u so omega(ubar) = mean of omega(u)
  // and |mean| <= mean|.|. That bound is the check that caught it.
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

    // Average only over cells that were water for essentially the WHOLE
    // period. Cells the interface swept through have a noisy ubar and, more
    // importantly, noisy neighbours -- and the vorticity stencil reads
    // neighbours. Kim's quantity is defined in the water, not at the surface.
    double vsum = 0., varea = 0., vmax = 0.;
    foreach (reduction(+:vsum) reduction(+:varea) reduction(max:vmax)) {
      // cs == 1, not cs > 0: Kim's own postprocessing masks on `sol_2D == 1`
      // (bio_stress.m:412), excluding embed CUT cells. Measured here: the
      // f-conditioned average alone left vor_max pinned at 241 non-dim, so the
      // extreme was never at the interface -- it is at the bag wall, where cut
      // cells carry huge velocity gradients. Same mask as the Kim-exact
      // tau/EDR columns.
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

    // [PROJECT ADDED, 2026-09-16] The steady-streaming FIELD, not just its
    // spatial mean. Kim's Fig. 4 is the mean-flow vorticity field and the
    // volume fraction field, each overlaid with streamlines of the mean flow,
    // and streamlines need ubar itself -- streaming.dat's single number cannot
    // produce them. Rewritten every window rather than accumulated: the file
    // that matters is the last (most settled) one, and at L10 each is 16 MB.
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
          // The window length is what f_acc is normalised by; without it a
          // reader cannot tell a cell that was water throughout from one the
          // interface swept.
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

event acceleration(i++)
{
  // [PROJECT CHANGED] Upstream's ramp is a plain LINEAR interpolation over a
  // fixed 30 PHYSICAL SECONDS, applied only to amplitude (phase/frequency
  // unramped): `Th_max2 = (Th_max/t_change_st)*t; Th = Th_max2*sin(w_bio_st*t);`
  // (see upstream BioReactor.c lines ~397-407). Replaced with a smooth-step
  // (3x²-2x³) interpolation of BOTH amplitude and phase, over N_RAMP_CYCLES
  // cycles instead of a fixed time, so that: (a) checkpoint-restart segments
  // that warm-start a NEW condition from an old one can interpolate between
  // two nonzero states (upstream's ramp only makes sense starting from
  // Th_max2=0, i.e. a true cold start) and (b) the ramp duration scales with
  // omega_b rather than being a fixed physical duration. See N_RAMP_CYCLES's
  // definition above for the still-open question of whether this changes
  // quasi-steady-state results relative to upstream's 30s linear ramp.
  // Smooth-step parameter interpolation: alpha: 0→1 over N_RAMP_CYCLES.
  // For fresh runs, *_prev fields are 0 → reproduces the original cold-start ramp.
  // For restarts, *_prev fields carry the previous segment's values → smooth transition
  // between two fully-forced steady states without ever underdriving the system.
  double elapsed  = t - t_ramp_start;
  double ramp_dur = N_RAMP_CYCLES * T_per_st;
  int    in_ramp  = elapsed < ramp_dur;
  double x_ss     = in_ramp ? elapsed / ramp_dur : 1.0;
  double alpha    = 3.*x_ss*x_ss - 2.*x_ss*x_ss*x_ss;   // smooth-step ∈ [0,1]
  // [PROJECT FIXED, 2026-09-07] d(alpha)/dt and d²(alpha)/dt², needed below.
  // Th_d/Th_2d previously assumed Ak(t)/phk(t) constant in time (correct
  // only once the ramp has finished, alpha'=alpha''=0) -- but during any
  // restart with theta_max_prev != theta_max, Ak(t) = Ak_prev+alpha(t)*dAk
  // genuinely varies over the ramp, and the old formulas silently dropped
  // the resulting Ak'(t)*sin(...) term entirely. This meant the velocity/
  // acceleration fed to the embedded-boundary condition was inconsistent
  // with the angle Th(t) actually being imposed, for the full N_RAMP_CYCLES
  // window, regardless of how small the parameter change was -- diary.md
  // 2026-09-07: an infinitesimal Delta_theta_max=0.1deg restart was found
  // to settle into a persistently different (14-23% offset) quasi-steady
  // tau_mean than an independent cold start at the IDENTICAL target
  // condition, while the trivial Delta_theta_max=0 identity-restart
  // control settled instantly and exactly -- consistent with a spurious
  // "kick" from this missing term tipping the flow into a different
  // attractor, present whenever alpha'(t) != 0 (i.e. whenever anything is
  // actually being ramped), independent of the ramp's amplitude.
  // d/dx[3x²-2x³] = 6x-6x², chain-ruled by dx_ss/dt = 1/ramp_dur (0 once
  // clamped past the ramp window, matching alpha itself being clamped).
  double dalpha_dt   = in_ramp ? 6.*x_ss*(1.-x_ss)/ramp_dur : 0.0;
  double d2alpha_dt2 = in_ramp ? 6.*(1.-2.*x_ss)/(ramp_dur*ramp_dur) : 0.0;

  // [PROJECT CHANGED] Upstream applies a single harmonic unconditionally:
  // `Th=Th_max*sin(w_bio_st*t); Th_d=w_bio_st*Th_max*cos(...); Th_2d=-w_bio_st²*
  // Th_max*sin(...)`. Generalized to a sum over params.n_harmonics so the
  // sweep/optimization pipeline can explore multi-harmonic rocking motions;
  // at n_harmonics=1 this reduces to exactly upstream's formula (verified
  // algebraically, diary.md 2026-07-28 — Ak=theta_max[0]*pi/180, phk=0).
  // Multi-harmonic angular forcing with smooth-step interpolation of amplitude and phase.
  // For each harmonic k: Ak and phk are interpolated from _prev → current over N_RAMP_CYCLES.
  // wk (frequency) is NOT ramped -- see N_RAMP_CYCLES's definition above;
  // that remains a real, separate, unfixed discontinuity at omega_b-
  // changing restarts, orthogonal to the Ak(t)/phk(t) fix below.
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

  // [PROJECT ADDED] Horizontal translational forcing has no upstream
  // equivalent at all — Kim et al.'s bioreactor only rocks (angular motion).
  // Added so the sweep/optimization pipeline can also explore horizontal
  // shaking/translation as a control parameter. Zero by default
  // (omega_h=0), which recovers upstream's pure-rocking motion exactly.
  // Multi-harmonic horizontal forcing with smooth-step interpolation.
  double x_acc = 0.;
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

  face vector av = a;
  // av.x: gravity + Coriolis + centrifugal + Euler (azimuthal) + horizontal translation
  // av.y: gravity + Coriolis + centrifugal + Euler (azimuthal)
  foreach_face(x)
    av.x[] = -sin(Th)/(Fr*Fr) + 2*Th_d*(u.y[] + u.y[-1,0])*0.5
    + Th_d*Th_d*(x+L_piv*sin(Th)) + Th_2d*(y+L_piv*cos(Th)) + x_acc;
  foreach_face(y)
    av.y[] = -cos(Th)/(Fr*Fr) - 2*Th_d*(u.x[] + u.x[0,-1])*0.5
    + Th_d*Th_d*(y+L_piv*cos(Th)) - Th_2d*(x+L_piv*sin(Th));
  a = av;
}
#endif

// [PROJECT ADDED] Reinstated upstream event (see REMOVE_DROP note near the
// tag.h include above). Runtime-gated on params.remove_drop instead of a
// compile-time #if, so a single binary can run both arms of the on/off
// matrix. Upstream ran this unconditionally once enabled, every step
// (i++); logic unchanged from upstream's remove_drop event.
event remove_drop(i++) {
  // `if` wrapper, not a bare `return;`: a valueless return in an event is a nonzero int return that
  // stops the run (CLAUDE.md). Same behaviour as before for remove_drop = 0 and 1.
  if (params.remove_drop) {
    remove_droplets(f, remove_minsize, remove_threshold, false);  // remove droplets
    remove_droplets(f, remove_minsize, remove_threshold, true);   // remove bubbles
  }
}

// [PROJECT ADDED, 2026-09-07] One-shot velocity-field "kick" test (diary.md
// 2026-09-07): direct test of whether the cold-start-beats-warm-start
// escape-to-a-different-branch behavior is a genuine property of the
// continuous dynamical system (any tiny perturbation escapes to the same
// branch, regardless of how it's introduced) or something specific to the
// checkpoint dump/restore code path (already ruled out as a slow-ramp
// artifact -- a 40-cycle quasi-static ramp landed on the same branch as a
// 3-cycle one). Enabled only when compiled with -DKICK_TEST=1; fires ONCE
// at t=KICK_T (nondim) in a single, otherwise-uninterrupted continuous run
// -- no dump/restore involved at all, so this isolates the continuous PDE's
// own dynamics from anything checkpoint-related.
//
// [PROJECT AMENDED, 2026-09-07] KICK_MODE. The original kick (now mode 0)
// was u *= (1+KICK_EPS): a UNIFORM rescale. That is a bad probe of
// stability, and its null result was over-read. A uniform rescale keeps the
// field divergence-free, preserves every symmetry of the state, and is very
// nearly tangent to the one-parameter family of periodic orbits -- i.e. it
// is close to the single direction guaranteed NOT to excite a transverse
// mode. Its decay says almost nothing.
//
// Mode 1 is the honest probe: u *= (1 + KICK_EPS*noise()) per cell, a
// spatially-random, zero-mean, KICK_EPS-relative perturbation that breaks
// the symmetries and projects onto every mode, including the ones a uniform
// rescale cannot reach. It is not divergence-free; the projection step at
// the next timestep handles that, which is precisely what a generic
// perturbation does.
//
// The question this settles: the 2x2 (diary.md 2026-09-07) showed the
// escape is driven by the restored STATE, not the ramp code path -- a
// restart whose state is ~0.3% off the target's own converged state
// diverges, while a bit-exact one and a 1e-7 forcing change do not. Either
// the reference state is linearly stable with a finite basin that 0.3%
// escapes (then a 0.3% random kick on a pure cold start must also escape,
// and nothing about this is restart-specific), or it is not (then the
// restart path is still doing something extra).
#ifndef KICK_TEST
#define KICK_TEST 0
#endif
#if KICK_TEST
#ifndef KICK_T
#define KICK_T 12.15
#endif
#ifndef KICK_EPS
#define KICK_EPS 1e-4
#endif
#ifndef KICK_MODE
#define KICK_MODE 0            // 0 = uniform rescale, 1 = per-cell random
#endif
event velocity_kick (t = KICK_T) {
  if (pid() == 0)
    fprintf (ferr, "kick: mode %d, relative amplitude %g, at t=%.4g\n",
             KICK_MODE, KICK_EPS, t);
  foreach() {
#if KICK_MODE == 1
    u.x[] *= (1. + KICK_EPS*noise());
    u.y[] *= (1. + KICK_EPS*noise());
#else
    u.x[] *= (1. + KICK_EPS);
    u.y[] *= (1. + KICK_EPS);
#endif
  }
  boundary ({u.x, u.y});
}
#endif
