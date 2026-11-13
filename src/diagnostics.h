// diagnostics.h -- Time series: events logstats, pressure_diagnostics, normcal (normf.dat, tr_oxy.dat, ...).
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

//  Log performance and runtime
#if REPLAY_FLUX
static double t_next_log = -1.;
event logstats (i++) {
  if (out_due (&t_next_log, 0.1)) {
#else
event logstats (t+=0.1; t <= t_end) {
#endif

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    if (pid() == 0){
      fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
      fflush(fp_stats);
    }
#if REPLAY_FLUX
  }
#endif
}

// Pressure Poisson residuals — health-check build only (DIAGNOSTICS=1)
#if DIAGNOSTICS && !REPLAY_FLUX
event pressure_diagnostics (t+=0.1; t<=t_end) {
  static FILE *fp_pdiag = NULL;
  if (!fp_pdiag) {
    fp_pdiag = fopen("pressure_diag.dat", "w");
    fprintf(fp_pdiag, "i t mgp_resa mgu_resa mgp_i mgu_i\n");
  }
  if (pid() == 0) {
    fprintf(fp_pdiag, "%d %g %g %g %d %d\n", i, t, mgp.resa, mgu.resa, mgp.i, mgu.i);
    fflush(fp_pdiag);
  }
}
#endif

// Compute simulation statistics
#if NORMCAL
event normcal (t+=t_out; t<=t_end){
    //timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
#if REC_CHECKSUM && RECORD_FLUX
    rec_checksum_to ("rec_checksum_norm.txt");
#endif

    scalar ux_liq[],uy_liq[],ux_liq_abs[],omega[],omega_liq[],oxy_liq[],f_liq[],posY[],c_liq[],c1_liq[],c2_liq[],c3_liq[];
    double omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max;
    double ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max;
    double f_liq_sum,f_liq_interf,posY_max,posY_min,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2;
    
    vorticity (u, omega); // vorticity

    // only liquid velocity
    foreach(){
      ux_liq[]  = u.x[]*f[];
      uy_liq[]  = u.y[]*f[];
      // BUG (upstream rcsc-group/BioReactor): oxy[]*f[] mixes gas-phase oxygen into
      // the integral at interface cells, inflating oxy_liq_sum above its physical bound.
      // FIX: extract only the liquid-phase contribution via Henry's law activity.
      // oxy[] = a*(f*alpha + (1-f))  where a = gas-equivalent activity.
      // Liquid-phase concentration (normalised by saturation = alpha) = f*a = f*oxy/(f*alpha+(1-f)).
      // Result: oxy_liq=0 in pure gas, oxy_liq=1 in pure liquid at saturation.
      oxy_liq[] = f[]*oxy[]/(f[]*c_oxy_alpha + (1.-f[]) + 1e-10);
      omega_liq[] = omega[]*f[];
      // BUG (upstream): (1-cs[])*f[] is zero everywhere inside the bag (cs=1 in fluid cells),
      // so statsf2(f_liq).sum ≈ 0.002 (only bag-wall cut cells) instead of ~0.28 (liquid volume).
      // Dividing oxy_liq_sum by this tiny number pushed C* >> 1 immediately.
      // FIX: use f[] directly — statsf2(f).sum = true liquid volume via embed-aware dv().
      f_liq[]   = f[];
      c_liq[]   = c[]*f[];
      c1_liq[]  = c1[]*f[];
      c2_liq[]  = c2[]*f[];
      c3_liq[]  = c3[]*f[];
    }

    position (f, posY, {0,1,0});  // (0,1,0) indicates the unit vector in the y-direction

    omega_liq_avg = normf(omega_liq).avg;
    omega_liq_rms = normf(omega_liq).rms;
    omega_liq_vol = normf(omega_liq).volume;
    omega_liq_max = normf(omega_liq).max;
    ux_liq_avg    = normf(ux_liq).avg;
    ux_liq_rms    = normf(ux_liq).rms;
    ux_liq_vol    = normf(ux_liq).volume;
    ux_liq_max    = normf(ux_liq).max;
    uy_liq_avg    = normf(uy_liq).avg;
    uy_liq_rms    = normf(uy_liq).rms;
    uy_liq_vol    = normf(uy_liq).volume;
    uy_liq_max    = normf(uy_liq).max;
    
    f_liq_sum     = statsf2(f_liq).sum;
    f_liq_interf  = interface_area(f);
    // [PROJECT ADDED, 2026-09-16] SIGNED spatial means of the liquid velocity.
    // Basilisk's normf(v).avg is sum(|v| dv)/volume -- an absolute value -- so
    // ux_liq_avg rectifies the signal: it oscillates at TWICE the rocking
    // frequency and never goes negative. Kim's Fig. 2 plots the signed means
    // against the rocking angle, and a rectified series cannot show the phase
    // relationship that figure exists to show.
    double ux_sig = 0., uy_sig = 0., wsum = 0.;
    foreach (reduction(+:ux_sig) reduction(+:uy_sig) reduction(+:wsum)) {
      double dA = dv();
      ux_sig += u.x[]*f[]*dA;
      uy_sig += u.y[]*f[]*dA;
      wsum   += dA;
    }
    double ux_liq_savg = wsum > 0. ? ux_sig/wsum : 0.;
    double uy_liq_savg = wsum > 0. ? uy_sig/wsum : 0.;

    // [PROJECT ADDED, 2026-09-16] The SPANWISE component, for Kim's Fig. A.18
    // (2D vs 3D). The solver compiles and runs unchanged under -grid=octree --
    // the embedded bag is a function of x and y only, so it extrudes into a
    // slab and the liquid volume comes out identical to 2D -- but nothing
    // logged u_z, so the one quantity the comparison turns on was the one
    // quantity missing. Zero by construction in 2D, where these two columns
    // are written as zeros to keep the file's shape the same at both
    // dimensionalities.
    double uz_liq_savg = 0., uz_liq_rms = 0.;
#if dimension > 2
    {
      double uz_sig = 0., uz_sq = 0.;
      foreach (reduction(+:uz_sig) reduction(+:uz_sq)) {
        double dA = dv();
        uz_sig += u.z[]*f[]*dA;
        uz_sq  += sq(u.z[]*f[])*dA;
      }
      if (wsum > 0.) {
        uz_liq_savg = uz_sig/wsum;
        uz_liq_rms  = sqrt (uz_sq/wsum);
      }
    }
#endif

    posY_max      = statsf(posY).max;

    // [PROJECT ADDED, 2026-09-16] Interface height at the two ENDS of the bag.
    // Kim's Fig. 14 is "the maximum surface elevation at the left end", and
    // Fig. 15 is the frequency spectrum of the surface elevation -- neither is
    // the global posY_max, which is wherever the crest happens to be and
    // wanders across the bag over a cycle. The probe is the outermost 2% of
    // the bag width; narrower than that and a coarse grid has no interface
    // cell in it on some timesteps, leaving the column empty at random.
    double x_lo = X0 + 0.02*L0, x_hi = X0 + L0 - 0.02*L0;
    double y_left = -HUGE, y_right = -HUGE;
    foreach (reduction(max:y_left) reduction(max:y_right)) {
      if (posY[] != nodata && cs[] > 0.) {
        if (x <= x_lo && posY[] > y_left)  y_left  = posY[];
        if (x >= x_hi && posY[] > y_right) y_right = posY[];
      }
    }
    if (y_left  == -HUGE) y_left  = nodata;
    if (y_right == -HUGE) y_right = nodata;
    posY_min      = statsf(posY).min;

    oxy_liq_sum   = statsf2(oxy_liq).sum;
    oxy_liq_sum2  = statsf2(oxy_liq).sum2;
    c_liq_sum     = statsf2(c_liq).sum;
    c_liq_sum2    = statsf2(c_liq).sum2;
    c1_liq_sum    = statsf2(c1_liq).sum;
    c1_liq_sum2   = statsf2(c1_liq).sum2;
    c2_liq_sum    = statsf2(c2_liq).sum;
    c2_liq_sum2   = statsf2(c2_liq).sum2;
    c3_liq_sum    = statsf2(c3_liq).sum;
    c3_liq_sum2   = statsf2(c3_liq).sum2;

    // [PROJECT ADDED] The entire shear-stress KPI pipeline (this block through
    // the tau_95/98/100/mean percentile walk and fp_tau write below) has no
    // upstream equivalent whatsoever — Kim et al.'s driver computes no shear
    // stress at all (it is presumably computed in their own postprocessing,
    // outside this file). Added (commits bc1c74a, 6d08f53) to reproduce their
    // Fig. 13a tau_liq_max/tau_liq_mean KPIs from this project's own sweeps.
    //
    // CAUTION (diary.md, 2026-07-28): the comment below claims this stencil
    // "mirrors vorticity() in basilisk/src/utils.h" — checked directly against
    // basilisk/src/utils.h:286-292, and that claim is only half true. The real
    // vorticity() weights by face metric factors fm.x/fm.y/cm and divides by
    // 2*(cm[]+SEPS)*Delta, which correct the finite-difference stencil near
    // embedded cut-cells (i.e. near walls — exactly where peak shear stress is
    // largest). This implementation uses plain 2*Delta with no metric
    // correction at all. This is a real, confirmed discrepancy from the
    // pattern it claims to follow. NOT YET FIXED: has not been shown to
    // explain the ~5x mean-VELOCITY mismatch found in the same investigation
    // (velocity is a primitive field, this bug can only affect a derivative
    // quantity), so it is at most a partial explanation, and the user's
    // explicit instruction is to make no changes beyond what's strictly
    // numerically necessary — not yet applied pending further diagnosis.
    //
    // ── 98th-percentile shear stress τ₉₈(t) ──────────────────────────────────
    // tau = |μ(∂u/∂y + ∂v/∂x)| in the liquid domain (f[] > 0.5).
    // Masked to bulk liquid to exclude interface artefacts.
    //
    // Implementation: two-pass 200-bin histogram — no dynamic memory, MPI-safe.
    //   Pass 1: global tau_max via statsf2 (Basilisk internally does MPI_Allreduce).
    //   Pass 2: local bin counts → MPI_Allreduce(SUM) → global percentile walk.
    //
    // Stencil mirrors vorticity() in basilisk/src/utils.h: uses face-centred
    // velocity values u.x[0,±1] and u.y[±1] for 2nd-order centred differences.
    // Pass 1: global max and mean via inline reduction (no scalar allocation).
    // Declaring scalar tau_liq[] inside an event leaks a Basilisk scalar on every
    // call and corrupts the global scalar list, causing segfaults at fidelity ≥7.
    // [PROJECT ADDED, 2026-08-07] tau_*_strict mirrors Kim et al.'s own
    // postprocessing (dev/postprocessing/bio_stress.m, shared directly):
    // `in_liq_field2 = find(abs(al_2D_proj) > 1-1e-10)` -- pure liquid only,
    // excluding ANY interface/mixed cell, vs. our existing f[]>0.5 which
    // includes the whole liquid-side half of the interface transition band.
    // Computed in the SAME pass as the existing (f[]>0.5) reduction to avoid
    // a second loop -- both masks share the identical tau formula, just a
    // stricter threshold gates the _strict accumulators. 1e-6 instead of
    // Kim's 1e-10: our f[] is a VOF fraction from a finite-precision solve,
    // not his post-hoc alpha field, so exact machine-precision equality to 1
    // is not expected; 1e-6 is still far tighter than the interface's
    // O(1) width in f, so it excludes the same cells his 1e-10 would.
    // [PROJECT ADDED, 2026-08-07] tau_max_signed mirrors Kim et al.'s own
    // formula LITERALLY: dev/postprocessing/bio_stress.m:345,
    // `tau_field = mu_field.*(duxdy_field(:,3) + duydx_field(:,3))` -- no
    // fabs() anywhere in his script -- vs. our existing fabs(du_dy+dv_dx).
    // max(signed) <= max(|signed|) always, and shear reverses sign every
    // half rocking-cycle, so this could plausibly explain us reading
    // persistently higher than Kim by a wide margin. Init to a large
    // negative sentinel, not 0, since the signed quantity can legitimately
    // be negative throughout a liquid region.
    // [PROJECT ADDED, 2026-08-08] Energy dissipation rate (EDR), for Fig. 8a
    // reproduction: epsilon = mu*[2*(du_x/dx)^2 + 2*(du_y/dy)^2 +
    // (du_x/dy+du_y/dx)^2] -- matches Kim et al.'s own formula exactly
    // (Main.tex, Sec. "Shear stress and energy dissipation rate", and their
    // shared bio_stress.m:346 -- Ediss_field = mu_field.*(2*duxdx.^2 +
    // 2*duydy.^2 + (duxdy+duydx).^2)). Only the domain-mean is needed here
    // (Fig. 8a plots spatially-averaged tau and epsilon vs. time, not a
    // max); reuses the exact same f[]>0.5 mask as tau for consistency.
    double tau_max_val = 0., tau_sum = 0., tau_vol = 0.;
    double tau_max_strict = 0., tau_sum_strict = 0., tau_vol_strict = 0.;
    // [PROJECT ADDED, 2026-09-15] Kim-EXACT mask, read from his own
    // postprocessing (rcsc-group/BioReactor3D dev/postprocessing/bio_stress.m,
    // private repo): line 412 `in_solid = find(sol_2D == 1)` and line 546
    // `in_liq_field2 = find(abs(al_2D_proj)>1-1e-10)` -- FULLY-fluid AND
    // FULLY-liquid cells only, excluding every interface cell and every embed
    // cut cell. He also takes max()/mean() of the SIGNED tau (line 583-585:
    // `tau_field_liq2_dim = tau_field_liq2*U_bio/L_bio`, no abs()), where our
    // tau_100 is max|tau|. And he reports an Ediss MAX (line 600) that this
    // project never logged at all -- Kim's Fig 13 has four series per panel.
    double tau_kim_max = -HUGE, tau_kim_sum = 0., kim_vol = 0.;
    double ediss_kim_max = 0., ediss_kim_sum = 0.;
    double tau_max_signed = -1e30;
    double ediss_sum = 0.;
    // [PROJECT ADDED, 2026-08-09] tau_sum_signed: Kim et al.'s Fig. 8a
    // plots <tau_w'>, the domain-mean of the SIGNED shear stress (their
    // Fig_tau_Ediss.pdf shows it oscillating through zero, -1.8e-3 to
    // +1.8e-3 Pa) -- not mean(|tau|). tau_mean_val below has always been
    // mean(|tau|) (tau_sum accumulates the fabs()'d `tau`, not
    // `tau_signed`); that's a genuinely different statistic, not just a
    // scale/units difference, and the one Fig. 8a needs to actually match
    // Kim's plotted quantity. Caught by the user pointing out Fig 8a's
    // reproduction didn't match Kim's real figure's shape (sign-crossing)
    // once actually rendered and viewed side by side, not just captioned.
    double tau_sum_signed = 0.;
    foreach(reduction(max:tau_max_val) reduction(+:tau_sum) reduction(+:tau_vol)
            reduction(max:tau_max_strict) reduction(+:tau_sum_strict) reduction(+:tau_vol_strict)
            reduction(max:tau_max_signed) reduction(+:ediss_sum) reduction(+:tau_sum_signed)
            reduction(max:tau_kim_max) reduction(+:tau_kim_sum) reduction(+:kim_vol)
            reduction(max:ediss_kim_max) reduction(+:ediss_kim_sum)) {
      // [PROJECT FIX, 2026-09-15, diary.md] Mask on cs[] too, and weight by
      // it. `f` is initialised with fraction(f, y_init - y), which fills the
      // lower half of the WHOLE DOMAIN -- including everything outside the
      // embedded bag, where there is no fluid and u is identically zero.
      // Masking on `f[] > 0.5` alone therefore swept those cells into the
      // spatial averages: they add nothing to the numerator but inflate the
      // area in the denominator. Measured at L10 (runs/57f68830): cells with
      // f>0.5 cover 0.5000 of the domain while the liquid actually inside the
      // bag covers 0.1425 -- a 3.51x inflation, so every reported spatial mean
      // of tau and EDR was ~3x too small. That is the entire "~3x below Kim"
      // discrepancy: with the bag mask applied, tau_liq_mean becomes 1.685e-3
      // Pa vs Kim's 1.803e-3 (0.93x) and Ediss_liq_mean 0.2773 vs 0.2866
      // W/m3 (0.97x). Spatial MAXIMA were never affected (u=0 outside the bag
      // => tau=0 there), which is why peak tau always looked reasonable while
      // the means did not, and why the deficit was resolution-independent.
      // cs[] is the embed fluid fraction, so cs[]*Delta*Delta is the true
      // fluid area of a cut cell.
      if (cs[] > 0. && f[] > 0.5) {
        double dA = cs[]*(Delta*Delta);
        double du_dy = (u.x[0,1] - u.x[0,-1]) / (2.*Delta);
        double dv_dx = (u.y[1]   - u.y[-1])   / (2.*Delta);
        double du_dx = (u.x[1]   - u.x[-1])   / (2.*Delta);
        double dv_dy = (u.y[0,1] - u.y[0,-1]) / (2.*Delta);
        double tau_signed = mu(f[]) * (du_dy + dv_dx);
        double tau   = fabs(tau_signed);
        double ediss = mu(f[]) * (2.*du_dx*du_dx + 2.*dv_dy*dv_dy + (du_dy+dv_dx)*(du_dy+dv_dx));
        if (tau > tau_max_val) tau_max_val = tau;
        if (tau_signed > tau_max_signed) tau_max_signed = tau_signed;
        tau_sum += tau * dA;
        tau_sum_signed += tau_signed * dA;
        tau_vol += dA;
        ediss_sum += ediss * dA;
        if (f[] > 1. - 1e-6) {
          if (tau > tau_max_strict) tau_max_strict = tau;
          tau_sum_strict += tau * dA;
          tau_vol_strict += dA;
        }
        // Kim's own mask: cs==1 (no cut cells) AND f==1 (no interface cells),
        // on the SIGNED tau. Area-weighted, which is identical to his
        // arithmetic mean() on his uniform mesh and correct on a graded one.
        if (cs[] > 1. - 1e-10 && f[] > 1. - 1e-10) {
          if (tau_signed > tau_kim_max) tau_kim_max = tau_signed;
          if (ediss > ediss_kim_max)    ediss_kim_max = ediss;
          tau_kim_sum   += tau_signed * dA;
          ediss_kim_sum += ediss * dA;
          kim_vol       += dA;
        }
      }
    }
    double tau_mean_val    = (tau_vol > 0.)        ? tau_sum / tau_vol               : 0.;
    double tau_mean_strict = (tau_vol_strict > 0.) ? tau_sum_strict / tau_vol_strict : 0.;
    double tau_kim_mean   = (kim_vol > 0.) ? tau_kim_sum   / kim_vol : 0.;
    double ediss_kim_mean = (kim_vol > 0.) ? ediss_kim_sum / kim_vol : 0.;
    if (kim_vol <= 0.) { tau_kim_max = 0.; ediss_kim_max = 0.; }
    double ediss_mean_val  = (tau_vol > 0.)        ? ediss_sum / tau_vol             : 0.;
    double tau_mean_signed = (tau_vol > 0.)        ? tau_sum_signed / tau_vol        : 0.;
    if (tau_max_val < 1e-14) tau_max_val = 1e-14;  // guard /0
    if (tau_max_strict < 1e-14) tau_max_strict = 1e-14;

    // Pass 2: 200-bin histogram — recomputes tau inline, no scalar needed
    // [FIX, 2026-08-20, diary.md] `bins[b]++` used to be a manual array
    // increment with no protection -- not a pattern Basilisk's foreach()/
    // OpenMP auto-reduction recognizes (unlike the scalar max/sum
    // accumulators above, which DO get an automatic reduction clause), so
    // under -fopenmp with >1 thread it was a genuine, unprotected data
    // race: concurrent non-atomic read-modify-write to the same bins[b]
    // silently loses increments. MPI builds never hit this (qcc disables
    // OpenMP under -D_MPI=1, confirmed by its own "OpenMP cannot be used
    // with MPI (yet)" warning, so each MPI rank is single-threaded here)
    // -- only OpenMP-only builds (the DEFAULT for chain.py sweeps that
    // don't explicitly set mpi:true) were affected. Confirmed empirically,
    // not just by inspection: two identical OpenMP runs (same params, same
    // binary) gave tau_95 values 38% apart and tau_98 values 6.3% apart,
    // while tau_100/tau_mean/ediss_mean (computed via the auto-reduced
    // accumulators) differed by <1% (ordinary floating-point reduction-
    // order noise, not a bug). Both `#pragma omp atomic` and Basilisk's own
    // OMP(omp atomic) pragma-insertion macro (grid/config.h) fail to
    // compile at this exact position -- qcc's stencil-analysis AST walk
    // chokes on a pragma inside this foreach() body (confirmed: same
    // "expected expression before '}' token" error from both). Worked
    // around with per-thread-local bins (flat 1D array, indexed by
    // omp_get_thread_num() -- a plain function call, not a pragma, so qcc
    // parses it fine) merged into the final histogram after the loop.
    // Re-verified after this fix: the same reproducibility test now gives
    // tau_95/98 differing by ~0.87% between identical runs, matching
    // tau_100/mean's ordinary noise level -- the race is gone.
    #define TAU_BINS 2000
    int _tau_nthreads = 1;
    #if _OPENMP
      _tau_nthreads = omp_get_max_threads();
    #endif
    long *_tau_tbins = (long *) calloc(_tau_nthreads * TAU_BINS, sizeof(long));
    foreach() {
      // same bag mask as the averaging pass above (2026-09-15 fix)
      if (cs[] > 0. && f[] > 0.5) {
        double du_dy = (u.x[0,1] - u.x[0,-1]) / (2.*Delta);
        double dv_dx = (u.y[1]   - u.y[-1])   / (2.*Delta);
        double tau   = mu(f[]) * fabs(du_dy + dv_dx);
        int b = (int)(tau / tau_max_val * (TAU_BINS - 1));
        if (b < 0)         b = 0;
        if (b >= TAU_BINS) b = TAU_BINS - 1;
        int _tid = 0;
        #if _OPENMP
          _tid = omp_get_thread_num();
        #endif
        _tau_tbins[_tid * TAU_BINS + b]++;
      }
    }
    long bins[TAU_BINS];
    for (int k = 0; k < TAU_BINS; k++) {
      bins[k] = 0;
      for (int th = 0; th < _tau_nthreads; th++)
        bins[k] += _tau_tbins[th * TAU_BINS + k];
    }
    free(_tau_tbins);
    #if _MPI
    {
      long gbins[TAU_BINS];
      MPI_Allreduce(bins, gbins, TAU_BINS, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      for (int k = 0; k < TAU_BINS; k++) bins[k] = gbins[k];
    }
    #endif

    // Walk bins once to find 95th and 98th percentiles; 100th == tau_max_val
    long total = 0, cumul = 0;
    for (int k = 0; k < TAU_BINS; k++) total += bins[k];
    double tau_95_val = tau_max_val;
    double tau_98_val = tau_max_val;
    int    found_95   = 0;
    // Bin k holds tau in [k, k+1)*tau_max_val/(TAU_BINS-1) (see the index
    // map above). The percentile is the linear interpolation of the cumulative
    // count inside the crossing bin: lower edge + fraction*width (2026-10-10;
    // the old rule returned the upper edge, biasing tau_95 high by up to a bin).
    double _tau_w = tau_max_val / (TAU_BINS - 1.0);
    for (int k = 0; k < TAU_BINS; k++) {
      long prev = cumul;
      cumul += bins[k];
      if (!found_95 && cumul >= (long)(0.95 * (double)total)) {
        double fr = bins[k] > 0 ? (0.95*(double)total - (double)prev) / (double)bins[k] : 0.;
        if (fr < 0.) fr = 0.; if (fr > 1.) fr = 1.;
        tau_95_val = _tau_w * (k + fr);
        found_95   = 1;
      }
      if (cumul >= (long)(0.98 * (double)total)) {
        double fr = bins[k] > 0 ? (0.98*(double)total - (double)prev) / (double)bins[k] : 0.;
        if (fr < 0.) fr = 0.; if (fr > 1.) fr = 1.;
        tau_98_val = _tau_w * (k + fr);
        break;
      }
    }
    #undef TAU_BINS

    // [PROJECT ADDED, 2026-11-04] edr_95: 95th percentile of the local EDR
    // ediss = mu*(2 du/dx^2 + 2 dv/dy^2 + (du/dy+dv/dx)^2), on the SAME mask as
    // tau_95 (cs>0, f>0.5), with the same stencil as the ediss_mean pass above.
    // Same two-pass 2000-bin histogram as tau_95 (per-thread bins, no pragma
    // inside foreach). Also writes the per-step spatial max, so that
    // edr_mean(t) <= edr_95(t) <= edr_max(t) can be checked sample by sample.
    // Collective: every rank runs it, outside the pid() == 0 block.
    double edr_max_val = 0.;
    foreach(reduction(max:edr_max_val)) {
      if (cs[] > 0. && f[] > 0.5) {
        double du_dx = (u.x[1]   - u.x[-1])   / (2.*Delta);
        double dv_dy = (u.y[0,1] - u.y[0,-1]) / (2.*Delta);
        double du_dy = (u.x[0,1] - u.x[0,-1]) / (2.*Delta);
        double dv_dx = (u.y[1]   - u.y[-1])   / (2.*Delta);
        double edr_cell = mu(f[]) * (2.*du_dx*du_dx + 2.*dv_dy*dv_dy + (du_dy+dv_dx)*(du_dy+dv_dx));
        if (edr_cell > edr_max_val) edr_max_val = edr_cell;
      }
    }
    if (edr_max_val < 1e-14) edr_max_val = 1e-14;  // guard /0

    #define EDR_BINS 2000
    int _edr_nthreads = 1;
    #if _OPENMP
      _edr_nthreads = omp_get_max_threads();
    #endif
    long *_edr_tbins = (long *) calloc(_edr_nthreads * EDR_BINS, sizeof(long));
    foreach() {
      if (cs[] > 0. && f[] > 0.5) {
        double du_dx = (u.x[1]   - u.x[-1])   / (2.*Delta);
        double dv_dy = (u.y[0,1] - u.y[0,-1]) / (2.*Delta);
        double du_dy = (u.x[0,1] - u.x[0,-1]) / (2.*Delta);
        double dv_dx = (u.y[1]   - u.y[-1])   / (2.*Delta);
        double edr_cell = mu(f[]) * (2.*du_dx*du_dx + 2.*dv_dy*dv_dy + (du_dy+dv_dx)*(du_dy+dv_dx));
        int eb = (int)(edr_cell / edr_max_val * (EDR_BINS - 1));
        if (eb < 0)         eb = 0;
        if (eb >= EDR_BINS) eb = EDR_BINS - 1;
        int _etid = 0;
        #if _OPENMP
          _etid = omp_get_thread_num();
        #endif
        _edr_tbins[_etid * EDR_BINS + eb]++;
      }
    }
    long ebins[EDR_BINS];
    for (int k = 0; k < EDR_BINS; k++) {
      ebins[k] = 0;
      for (int th = 0; th < _edr_nthreads; th++)
        ebins[k] += _edr_tbins[th * EDR_BINS + k];
    }
    free(_edr_tbins);
    #if _MPI
    {
      long gebins[EDR_BINS];
      MPI_Allreduce(ebins, gebins, EDR_BINS, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      for (int k = 0; k < EDR_BINS; k++) ebins[k] = gebins[k];
    }
    #endif
    long etotal = 0, ecumul = 0;
    for (int k = 0; k < EDR_BINS; k++) etotal += ebins[k];
    double edr_95_val = edr_max_val;
    int    e_found    = 0;
    double _edr_w = edr_max_val / (EDR_BINS - 1.0);
    for (int k = 0; k < EDR_BINS; k++) {
      long eprev = ecumul;
      ecumul += ebins[k];
      if (!e_found && ecumul >= (long)(0.95 * (double)etotal)) {
        double fr = ebins[k] > 0 ? (0.95*(double)etotal - (double)eprev) / (double)ebins[k] : 0.;
        if (fr < 0.) fr = 0.; if (fr > 1.) fr = 1.;
        edr_95_val = _edr_w * (k + fr);
        e_found    = 1;
        break;
      }
    }
    #undef EDR_BINS

   double mb0, mb1, mb2;   // collective: every rank, outside the pid() == 0 block
   mix_bulk_sums (&mb0, &mb1, &mb2);
   // i, timestep, no of cells, real time elapsed, cpu time
   if (pid() == 0){

      fprintf(fp_norm, "%i %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",i,t,omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max,ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max,ux_liq_savg,uy_liq_savg,uz_liq_savg,uz_liq_rms);
      fflush(fp_norm);

      fprintf(fp_stats2, "%i %g %g %g %g %g %g %g \n",i,t,f_liq_sum,f_liq_interf,posY_max,posY_min,y_left,y_right);
      fflush(fp_stats2);

      fprintf(fp_stats3, "%i %g %g %g %g %g %g %g %g %g %g %g %.17g %.17g %.17g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2,mb0,mb1,mb2);
      //fprintf(fp_stats3, "%i %g %g %g %g %g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2);
      fflush(fp_stats3);

      fprintf(fp_tau, "%i %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n", i, t, tau_95_val, tau_98_val, tau_max_val, tau_mean_val, tau_max_strict, tau_mean_strict, tau_max_signed, ediss_mean_val, tau_mean_signed, tau_kim_max, tau_kim_mean, ediss_kim_max, ediss_kim_mean);
      fflush(fp_tau);

      fprintf(fp_edr, "%i %g %g %g \n", i, t, edr_95_val, edr_max_val);
      fflush(fp_edr);
   }
}
#endif

#if REPLAY_FLUX
// Replay version of event normcal: the interface, tracer and oxygen columns only. The formulas and the
// output format of vol_frac_interf.dat and tr_oxy.dat are copied from the event above and must stay in
// step with it; normf.dat and shear_stress.dat (velocity statistics) get their header and no rows.
static double t_next_norm = -1.;
event normcal (i++) {
  if (out_due (&t_next_norm, t_out)) {
#if REC_CHECKSUM
    rec_checksum_to ("rec_checksum_norm.txt");
#endif
    scalar oxy_liq[],f_liq[],posY[],c_liq[],c1_liq[],c2_liq[],c3_liq[];
    foreach(){
      oxy_liq[] = f[]*oxy[]/(f[]*c_oxy_alpha + (1.-f[]) + 1e-10);
      f_liq[]   = f[];
      c_liq[]   = c[]*f[];
      c1_liq[]  = c1[]*f[];
      c2_liq[]  = c2[]*f[];
      c3_liq[]  = c3[]*f[];
    }
    position (f, posY, {0,1,0});
    double f_liq_sum = statsf2(f_liq).sum;
    double f_liq_interf = interface_area(f);
    double posY_max = statsf(posY).max;
    double x_lo = X0 + 0.02*L0, x_hi = X0 + L0 - 0.02*L0;
    double y_left = -HUGE, y_right = -HUGE;
    foreach (reduction(max:y_left) reduction(max:y_right)) {
      if (posY[] != nodata && cs[] > 0.) {
        if (x <= x_lo && posY[] > y_left)  y_left  = posY[];
        if (x >= x_hi && posY[] > y_right) y_right = posY[];
      }
    }
    if (y_left  == -HUGE) y_left  = nodata;
    if (y_right == -HUGE) y_right = nodata;
    double posY_min = statsf(posY).min;
    double oxy_liq_sum   = statsf2(oxy_liq).sum;
    double oxy_liq_sum2  = statsf2(oxy_liq).sum2;
    double c_liq_sum     = statsf2(c_liq).sum;
    double c_liq_sum2    = statsf2(c_liq).sum2;
    double c1_liq_sum    = statsf2(c1_liq).sum;
    double c1_liq_sum2   = statsf2(c1_liq).sum2;
    double c2_liq_sum    = statsf2(c2_liq).sum;
    double c2_liq_sum2   = statsf2(c2_liq).sum2;
    double c3_liq_sum    = statsf2(c3_liq).sum;
    double c3_liq_sum2   = statsf2(c3_liq).sum2;
    double mb0, mb1, mb2;   // collective: every rank, outside the pid() == 0 block
    mix_bulk_sums (&mb0, &mb1, &mb2);
    if (pid() == 0) {
      fprintf(fp_stats2, "%i %g %g %g %g %g %g %g \n",i,t,f_liq_sum,f_liq_interf,posY_max,posY_min,y_left,y_right);
      fflush(fp_stats2);
      fprintf(fp_stats3, "%i %g %g %g %g %g %g %g %g %g %g %g %.17g %.17g %.17g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2,mb0,mb1,mb2);
      fflush(fp_stats3);
    }
  }
}

#endif // REPLAY_FLUX
