// diagnostics.h -- time series: events logstats, pressure_diagnostics, normcal.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.
// Output columns: docs_site/reference/output-files.md.

// logstats.dat: step, time, cells, wall clock (every 0.1 time units)
#if REPLAY_FLUX
static double t_next_log = -1.;
event logstats (i++) {
  if (out_due (&t_next_log, 0.1)) {
#else
event logstats (t+=0.1; t <= t_end) {
#endif

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
    if (pid() == 0){
      fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
      fflush(fp_stats);
    }
#if REPLAY_FLUX
  }
#endif
}

// pressure_diag.dat: multigrid residuals of the pressure and viscous solves (-DDIAGNOSTICS=1 only)
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

// ---------------------------------------------------------------------------------------------------
// Percentiles from a histogram. Each OpenMP thread fills its own bins (tbins[thread*nbins + bin]):
// a shared bins[b]++ inside foreach() is a data race (qcc adds reductions for scalars, not for arrays,
// and an omp pragma inside foreach() does not compile). merge_bins() sums the threads, then the MPI ranks.
// ---------------------------------------------------------------------------------------------------
static void merge_bins (long * tbins, int nthreads, long * bins, int nbins)
{
  for (int k = 0; k < nbins; k++) {
    bins[k] = 0;
    for (int th = 0; th < nthreads; th++)
      bins[k] += tbins[th * nbins + k];
  }
#if _MPI
  long gbins[nbins];
  MPI_Allreduce(bins, gbins, nbins, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  for (int k = 0; k < nbins; k++) bins[k] = gbins[k];
#endif
}

// Bin k holds values in [k, k+1)*width. The q-quantile is the lower edge of the bin where the cumulative
// count reaches q*total, plus the linear fraction inside that bin. Returns `fallback` if no bin reaches it.
static double bin_percentile (const long * bins, int nbins, double q, double width, double fallback)
{
  long total = 0, cumul = 0;
  for (int k = 0; k < nbins; k++) total += bins[k];
  for (int k = 0; k < nbins; k++) {
    long prev = cumul;
    cumul += bins[k];
    if (cumul >= (long)(q * (double)total)) {
      double fr = bins[k] > 0 ? (q*(double)total - (double)prev) / (double)bins[k] : 0.;
      if (fr < 0.) fr = 0.; if (fr > 1.) fr = 1.;
      return width * (k + fr);
    }
  }
  return fallback;
}

// ---------------------------------------------------------------------------------------------------
// event normcal: every t_out, one row in each of normf.dat, vol_frac_interf.dat, tr_oxy.dat,
// shear_stress.dat, edr_pct.dat.
//
// The liquid region is cs > 0 (inside the bag) AND f > 0.5 (liquid side). Both masks are needed: f is
// initialised over the whole domain, also outside the bag, where u = 0; without cs those cells inflate
// the area of every spatial mean (by 3.5x at L10). Cut cells are weighted by their fluid area cs*Delta^2.
//
// Velocity gradients use plain centred differences (2 Delta). Note: Basilisk's vorticity() also applies
// the embed metric (fm, cm) in cut cells; this stencil does not, so values in cut cells next to the wall
// are less accurate.
// ---------------------------------------------------------------------------------------------------
#if NORMCAL
event normcal (t+=t_out; t<=t_end){

#if REC_CHECKSUM && RECORD_FLUX
    rec_checksum_to ("rec_checksum_norm.txt");
#endif

    scalar ux_liq[],uy_liq[],ux_liq_abs[],omega[],omega_liq[],oxy_liq[],f_liq[],posY[],c_liq[],c1_liq[],c2_liq[],c3_liq[];
    double omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max;
    double ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max;
    double f_liq_sum,f_liq_interf,posY_max,posY_min,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2;

    vorticity (u, omega);

    foreach(){
      ux_liq[]  = u.x[]*f[];
      uy_liq[]  = u.y[]*f[];
      // Dissolved oxygen, liquid phase only. oxy holds a*(f*alpha + 1 - f), with a the gas-equivalent
      // activity; the liquid part, normalised by saturation, is f*a. (oxy*f would add gas-phase oxygen
      // in interface cells.) 0 in gas, 1 in liquid at saturation.
      oxy_liq[] = f[]*oxy[]/(f[]*c_oxy_alpha + (1.-f[]) + 1e-10);
      omega_liq[] = omega[]*f[];
      f_liq[]   = f[];        // liquid volume: statsf2(f).sum (Kim et al. used (1-cs)*f, which is ~0 in the bag)
      c_liq[]   = c[]*f[];
      c1_liq[]  = c1[]*f[];
      c2_liq[]  = c2[]*f[];
      c3_liq[]  = c3[]*f[];
    }

    position (f, posY, {0,1,0});  // interface height

    // normf().avg is the mean of |v| (rectified); the signed means are computed below.
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

    // Signed mean liquid velocity (Kim et al., Fig. 2 plots it against the rocking angle).
    double ux_sig = 0., uy_sig = 0., wsum = 0.;
    foreach (reduction(+:ux_sig) reduction(+:uy_sig) reduction(+:wsum)) {
      double dA = dv();
      ux_sig += u.x[]*f[]*dA;
      uy_sig += u.y[]*f[]*dA;
      wsum   += dA;
    }
    double ux_liq_savg = wsum > 0. ? ux_sig/wsum : 0.;
    double uy_liq_savg = wsum > 0. ? uy_sig/wsum : 0.;

    // Spanwise velocity, for 3D runs (-grid=octree; Kim et al., Fig. A.18). Zero columns in 2D.
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

    // Interface height at the two ends of the bag (outermost 2% of the width; narrower and a coarse grid
    // sometimes has no interface cell there). Kim et al., Figs. 14-15.
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

    // ---- Shear stress tau = mu (du/dy + dv/dx) and dissipation rate
    //      ediss = mu (2 (du/dx)^2 + 2 (dv/dy)^2 + (du/dy + dv/dx)^2), as in Kim et al. ----
    // Pass 1: maxima and area-weighted means, with three masks:
    //   default  cs > 0, f > 0.5 (as above), on |tau|
    //   _strict  also f > 1 - 1e-6: pure liquid, no interface cells
    //   _kim     cs = 1 and f = 1 (no cut cells, no interface cells), on the SIGNED tau: the mask and
    //            sign convention of Kim et al.'s own postprocessing
    // tau_mean_signed is the mean of the signed tau (Kim et al., Fig. 8a: it oscillates through zero).
    // (No scalar field is declared here for tau: reductions only.)
    double tau_max_val = 0., tau_sum = 0., tau_vol = 0.;
    double tau_max_strict = 0., tau_sum_strict = 0., tau_vol_strict = 0.;
    double tau_kim_max = -HUGE, tau_kim_sum = 0., kim_vol = 0.;
    double ediss_kim_max = 0., ediss_kim_sum = 0.;
    double tau_max_signed = -1e30;
    double ediss_sum = 0.;
    double tau_sum_signed = 0.;
    foreach(reduction(max:tau_max_val) reduction(+:tau_sum) reduction(+:tau_vol)
            reduction(max:tau_max_strict) reduction(+:tau_sum_strict) reduction(+:tau_vol_strict)
            reduction(max:tau_max_signed) reduction(+:ediss_sum) reduction(+:tau_sum_signed)
            reduction(max:tau_kim_max) reduction(+:tau_kim_sum) reduction(+:kim_vol)
            reduction(max:ediss_kim_max) reduction(+:ediss_kim_sum)) {
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

    // Pass 2: tau_95, tau_98 from a 2000-bin histogram of |tau| on [0, tau_max]
    #define TAU_BINS 2000
    int _tau_nthreads = 1;
    #if _OPENMP
      _tau_nthreads = omp_get_max_threads();
    #endif
    long *_tau_tbins = (long *) calloc(_tau_nthreads * TAU_BINS, sizeof(long));
    foreach() {
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
    merge_bins (_tau_tbins, _tau_nthreads, bins, TAU_BINS);
    free(_tau_tbins);
    double _tau_w = tau_max_val / (TAU_BINS - 1.0);
    double tau_95_val = bin_percentile (bins, TAU_BINS, 0.95, _tau_w, tau_max_val);
    double tau_98_val = bin_percentile (bins, TAU_BINS, 0.98, _tau_w, tau_max_val);
    #undef TAU_BINS

    // edr_95: 95th percentile of the local dissipation rate, same mask and method as tau_95.
    // edr_max is written too, so that edr_mean <= edr_95 <= edr_max can be checked per sample.
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
    merge_bins (_edr_tbins, _edr_nthreads, ebins, EDR_BINS);
    free(_edr_tbins);
    double edr_95_val = bin_percentile (ebins, EDR_BINS, 0.95, edr_max_val / (EDR_BINS - 1.0), edr_max_val);
    #undef EDR_BINS

   double mb0, mb1, mb2;   // collective: every rank, outside the pid() == 0 block
   mix_bulk_sums (&mb0, &mb1, &mb2);
   if (pid() == 0){

      fprintf(fp_norm, "%i %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n",i,t,omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max,ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max,ux_liq_savg,uy_liq_savg,uz_liq_savg,uz_liq_rms);
      fflush(fp_norm);

      fprintf(fp_stats2, "%i %g %g %g %g %g %g %g \n",i,t,f_liq_sum,f_liq_interf,posY_max,posY_min,y_left,y_right);
      fflush(fp_stats2);

      fprintf(fp_stats3, "%i %g %g %g %g %g %g %g %g %g %g %g %.17g %.17g %.17g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2,mb0,mb1,mb2);
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
