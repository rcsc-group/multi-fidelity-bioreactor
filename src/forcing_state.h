// forcing_state.h -- the rocking forcing at any time, the checkpoint label, and the source of a restart's ramp.
// Part of BioReactor.c: functions only (no events), included once after the global declarations.
//
// After every start the forcing ramps from a SOURCE forcing to the TARGET forcing (params.json) over
// N_RAMP_CYCLES target periods, with the smooth step alpha = 3x^2 - 2x^3. Amplitude, phase AND frequency ramp.
//   fresh start: source = target frequency, zero amplitude (the bag starts at rest);
//   restart:     source = the forcing that made the restored state, read from the dump's LABEL
//                (checkpoint.phase). A continuation (same forcing) therefore does not ramp at all.
// params.json holds the target only. A restart with no label, or with a label AND *_prev keys in params.json,
// stops with an error: one source of truth, nothing silent. Pitfall "restart-ramp" in
// docs_site/developer/solver-pitfalls.md; design in docs/plans/2026-11-13-restart-label-plan.md.
//
// The angle uses an integrated phase, psi(t) = 2 pi cycles0 + (integral of w dt), so a frequency change has
// no chirp: the rate of the angle is the instantaneous w, and the angle and its rate are continuous at a restart.

typedef struct {
  double omega_b;                 // rad/s
  double theta_max[N_MAX];        // deg
  double phi_angular[N_MAX];      // rad
  double omega_h;                 // rad/s
  double amplitude_h[N_MAX];      // m
  double phi_horizontal[N_MAX];   // rad
  double cycles, cycles_h;        // forcing phase at the start of the ramp, in periods (unit-free)
  double U_bio;                   // m/s: velocity unit of the restored fields (0 on a fresh start)
} ForcingSource;

static ForcingSource fsrc;        // set in set_run_window()

typedef struct {
  double Th, Th_d, Th_2d;         // angle (rad) and its first two time derivatives, run units
  double x_acc;                   // horizontal-shaking acceleration, run units
  double alpha;                   // ramp fraction, 0 -> 1
  double psi, psi_h;              // phase of the fundamental (rad), rocking and horizontal shaking
  double omega_b, omega_h;        // instantaneous frequencies (rad/s)
  double theta_max[N_MAX], phi_angular[N_MAX], amplitude_h[N_MAX], phi_horizontal[N_MAX];
} ForcingNow;

// The functions below have no branches (fmin/fmax and 0/1 weights instead of if, ?: and continue): qcc's
// dimension check runs them with undefined values and splits its run at every undefined condition, which used
// up its call budget and failed the build with exit code 2 and no message (2026-11-13). Pitfall "qcc-dimension-check".

// Integral of the smooth step over [0, tau] for a ramp of length D: D (x^3 - x^4/2), then D/2 + (tau - D).
static double ramp_integral (double tau, double D)
{
  double x = fmin (tau/D, 1.);
  return D*(x*x*x - 0.5*x*x*x*x) + fmax (tau - D, 0.);
}

// The rocking phase psi(t) = 2 pi cycles0 + w0 tau + dw (integral of alpha) and its rate w0 + dw alpha.
// Light on purpose (no harmonic loop): t_at_psi() iterates on it.
static double psi_at (double t, double * psi_d)
{
  double D   = N_RAMP_CYCLES * T_per_st;
  double tau = t - t_ramp_start;
  double x   = fmin (tau/D, 1.);
  double w0  = fsrc.omega_b * T_bio, dw = w_bio_st - w0;
  *psi_d = w0 + dw*(3.*x*x - 2.*x*x*x);
  return 2.*pi*fsrc.cycles + w0*tau + dw*ramp_integral (tau, D);
}

// The forcing at time t (run units). Pure: the acceleration event, the label and the stop rule all use it.
// With source = target (fresh start, or a continuation) it gives the same numbers as the old formula
// Th = sum Ak sin(k w (t + offset) + phk).
static ForcingNow forcing_now (double t)
{
  ForcingNow F = {0};
  double D        = N_RAMP_CYCLES * T_per_st;
  double tau      = t - t_ramp_start;
  double in_ramp  = tau < D;                     // 1 during the ramp, 0 after
  double x        = fmin (tau/D, 1.);
  double alpha    = 3.*x*x - 2.*x*x*x;
  double alpha_d  = 6.*x*(1. - x)/D;               // 0 at x = 1
  double alpha_2d = in_ramp*6.*(1. - 2.*x)/(D*D);
  F.alpha = alpha;

  // rocking: phase psi = 2 pi cycles0 + w0 tau + dw (integral of alpha), rate w0 + dw alpha
  double w0 = fsrc.omega_b * T_bio, dw = w_bio_st - w0;
  F.psi = 2.*pi*fsrc.cycles + w0*tau + dw*ramp_integral (tau, D);
  double psi_d  = w0 + dw*alpha;
  double psi_2d = dw*alpha_d;
  F.omega_b = (1. - alpha)*fsrc.omega_b + alpha*w_bio;
  // Loops run to the constant N_MAX; on = 0 switches off the unused harmonics (a loop bounded by
  // params.n_harmonics has no bound in the dimension check, which does not read params.json).
  for (int k = 1; k <= N_MAX; k++) {
    double on     = k <= params.n_harmonics;
    double dAk    = params.theta_max[k-1]   - fsrc.theta_max[k-1];     // deg
    double dphk   = params.phi_angular[k-1] - fsrc.phi_angular[k-1];   // rad
    double Ak     = on*(fsrc.theta_max[k-1] + alpha*dAk) * pi / 180.;
    double phk    =  fsrc.phi_angular[k-1] + alpha*dphk;
    double Ak_d   = on*alpha_d  * dAk * pi / 180.;
    double Ak_2d  = on*alpha_2d * dAk * pi / 180.;
    double arg    = k*F.psi + phk;
    double arg_d  = k*psi_d  + alpha_d*dphk;
    double arg_2d = k*psi_2d + alpha_2d*dphk;
    double s = sin(arg), c = cos(arg);
    F.Th    += Ak*s;
    F.Th_d  += Ak_d*s + Ak*arg_d*c;
    F.Th_2d += Ak_2d*s + 2.*Ak_d*arg_d*c + Ak*arg_2d*c - Ak*arg_d*arg_d*s;
    F.theta_max[k-1]   = on*((1. - alpha)*fsrc.theta_max[k-1]   + alpha*params.theta_max[k-1]);
    F.phi_angular[k-1] = on*((1. - alpha)*fsrc.phi_angular[k-1] + alpha*params.phi_angular[k-1]);
  }

  // horizontal shaking (off unless omega_h > 0; with omega_h = 0 the (k w)^2 factor makes it exactly 0): same
  // integrated phase. The acceleration keeps its constant-frequency form A (k w)^2 sin; the ramp-rate terms
  // are left out (no run uses this yet).
  double wh0 = fsrc.omega_h * T_bio, dwh = params.omega_h * T_bio - wh0;
  F.psi_h = 2.*pi*fsrc.cycles_h + wh0*tau + dwh*ramp_integral (tau, D);
  double psi_h_d = wh0 + dwh*alpha;
  F.omega_h = (1. - alpha)*fsrc.omega_h + alpha*params.omega_h;
  for (int k = 1; k <= N_MAX; k++) {
    double on  = k <= params.n_harmonics;
    double Ah  = on*((1. - alpha)*fsrc.amplitude_h[k-1]    + alpha*params.amplitude_h[k-1]);
    double phh = on*((1. - alpha)*fsrc.phi_horizontal[k-1] + alpha*params.phi_horizontal[k-1]);
    F.x_acc += (Ah / L_bio) * sq(k*psi_h_d) * sin(k*F.psi_h + phh);
    F.amplitude_h[k-1] = Ah;
    F.phi_horizontal[k-1] = phh;
  }
  return F;
}

// First time t >= t_ramp_start at which the rocking phase reaches `target` (psi rises with t). Closed form
// after the ramp; inside it, Newton on psi_at() from the linear guess, clamped to the ramp. The loop has a
// fixed count and calls only psi_at(): qcc's dimension check runs this code with undefined values and
// follows every undefined branch, so a long loop over forcing_now() (an earlier bisection, 200 steps)
// used up its call budget and failed the build with exit code 2 and no message (2026-11-13).
static double t_at_psi (double target)
{
  double D = N_RAMP_CYCLES * T_per_st, rate;
  double psi_D = psi_at (t_ramp_start + D, &rate);
  if (target >= psi_D)
    return t_ramp_start + D + (target - psi_D)/w_bio_st;
  double psi0 = psi_at (t_ramp_start, &rate);
  double tau = (target - psi0)/rate;
  for (int it = 0; it < 20; it++) {
    tau = tau < 0. ? 0. : tau > D ? D : tau;
    double psi = psi_at (t_ramp_start + tau, &rate);
    tau -= (psi - target)/rate;
  }
  return t_ramp_start + (tau < 0. ? 0. : tau > D ? D : tau);
}

// ---- the label: checkpoint.phase, written next to every dump ----------------------------------------------
// One JSON line: the dump time t, the forcing phase in periods (cycles, cycles_h), the units of the dumped
// fields (U_bio, T_bio, SI), the forcing at that instant (SI and degrees; inside a ramp, the ramped values),
// the angle and its rate (Th, Th_d in rad/s, for the continuity guard), and t_phase_offset/T_per_st for
// binaries older than this label (they read only those two).
// The forcing the acceleration event computed last, and its time t_forcing (one step before a dump: see
// label_write).
static ForcingNow _forcing_last;
static double     _forcing_last_t = -1.;
static int        _restart_forcing_written = 0;   // restart_forcing.json written (forcing.h)

static void label_write (FILE * fp, double t_now)
{
  // The forcing cache (F) is one step old at a dump: the acceleration event runs after the dump event in a
  // step, and the last step's never runs (2026-11-13: one-dt phase lag after every restart; guard S12). The
  // phase is therefore taken AT t_now from psi_at() (light, exact); Th, Th_d, Th_2d stay at t_forcing and are
  // written with it (guard S13 extrapolates). No forcing_now() call here: qcc's dimension check runs this
  // event and a second forcing_now() path used up its call budget (exit 2, pitfall "qcc-dimension-check").
  ForcingNow F = _forcing_last;
  double rate;
  double cycles = psi_at (t_now, &rate)/(2.*pi);
  double cycles_h = (F.psi_h + F.omega_h*T_bio*(t_now - _forcing_last_t))/(2.*pi);   // exact after the ramp
  fprintf (fp, "{\"t\": %.17g, \"t_forcing\": %.17g, \"cycles\": %.17g, \"cycles_h\": %.17g, \"U_bio\": %.17g, \"T_bio\": %.17g, "
           "\"omega_b\": %.17g, \"omega_h\": %.17g, \"n_harmonics\": %d, ",
           t_now, _forcing_last_t, cycles, cycles_h, U_bio, T_bio, F.omega_b, F.omega_h, params.n_harmonics);
  fprintf (fp, "\"theta_max\": [%.17g, %.17g, %.17g], \"phi_angular\": [%.17g, %.17g, %.17g], "
           "\"amplitude_h\": [%.17g, %.17g, %.17g], \"phi_horizontal\": [%.17g, %.17g, %.17g], ",
           F.theta_max[0], F.theta_max[1], F.theta_max[2], F.phi_angular[0], F.phi_angular[1], F.phi_angular[2],
           F.amplitude_h[0], F.amplitude_h[1], F.amplitude_h[2],
           F.phi_horizontal[0], F.phi_horizontal[1], F.phi_horizontal[2]);
  fprintf (fp, "\"Th\": %.17g, \"Th_d\": %.17g, \"Th_2d\": %.17g, \"T_per_st\": %.17g, \"t_phase_offset\": %.17g}\n",
           F.Th, F.Th_d/T_bio, F.Th_2d/(T_bio*T_bio), T_per_st, cycles*T_per_st - t_now);
}

// One value of the label: the number after "key": (0 and *found = 0 if the key is absent). The label is our own
// one-line format (label_write), so a key lookup is enough. No token loop: qcc's dimension check runs this
// with an unknown file, and a loop over an unknown token count used up its call budget (2026-11-13).
static double label_value (const char * json, const char * key, int * found)
{
  char k[64];
  snprintf (k, sizeof(k), "\"%s\":", key);
  const char * c = strstr (json, k);
  *found = c != NULL;
  return c ? atof (c + strlen (k)) : 0.;
}

// A 3-vector of the label, "key": [a, b, c]; zeros if absent.
static void label_vector (const char * json, const char * key, double * v)
{
  char k[64];
  snprintf (k, sizeof(k), "\"%s\":", key);
  const char * c = strstr (json, k);
  v[0] = v[1] = v[2] = 0.;
  if (c)
    sscanf (c + strlen (k), " [%lf , %lf , %lf ]", &v[0], &v[1], &v[2]);
}

static double _label_t = -1.;    // the time the label was written at (label_read); must be the dump's time

// Read a label. Returns 0 if the file is missing, -1 if it holds no "t", 1 for a label with the forcing, 2 for
// an old label (only t, t_phase_offset, T_per_st: written before 2026-11-13). Fills s as far as it can.
static int label_read (const char * path, ForcingSource * s)
{
  FILE * fp = fopen (path, "r");
  if (!fp)
    return 0;
  char json[4096];
  memset (json, 0, sizeof(json));
  if (fgets (json, sizeof(json), fp) == NULL)
    json[0] = '\0';
  fclose (fp);
  ForcingSource z = {0};
  *s = z;
  int has_t, has_forcing, found;
  double t = label_value (json, "t", &has_t);
  if (!has_t)
    return -1;
  _label_t = t;
  s->omega_b  = label_value (json, "omega_b", &has_forcing);
  s->omega_h  = label_value (json, "omega_h", &found);
  s->cycles   = label_value (json, "cycles", &found);
  s->cycles_h = label_value (json, "cycles_h", &found);
  s->U_bio    = label_value (json, "U_bio", &found);
  label_vector (json, "theta_max", s->theta_max);
  label_vector (json, "phi_angular", s->phi_angular);
  label_vector (json, "amplitude_h", s->amplitude_h);
  label_vector (json, "phi_horizontal", s->phi_horizontal);
  if (has_forcing)
    return 1;
  double offset = label_value (json, "t_phase_offset", &found);
  double T_old  = label_value (json, "T_per_st", &found);
  s->cycles = T_old > 0. ? (t + offset)/T_old : 0.;   // old label: the phase only
  return 2;
}

// A refused restart sets this flag and main() returns 1 right after set_run_window(). Not exit(): qcc's
// dimension check runs main() in its interpreter, where fopen() is always NULL (ast/interpreter/overload.h),
// so the label is "missing" there; an exit() on that path fails the build with exit code 2 and no message.
static int _restart_refused = 0;
static int _restart_kind    = 0;   // 1: source from the label, 2: old label + *_prev

static void restart_refuse (const char * why)
{
  if (pid() == 0)
    fprintf (stderr, "ERROR: restart refused: %s\n"
             "  The dump's label (checkpoint.phase, written next to every dump) gives the forcing that made the\n"
             "  restored state; params.json gives the target only. Stage the label with the dump\n"
             "  (scripts/submit_continue.stage_phase). Pitfall \"restart-ramp\" in docs_site/developer/"
             "solver-pitfalls.md.\n", why);
  _restart_refused = 1;
}

// The source of a restart's ramp, from the label in the working directory (rules at the top of this file).
// Returns 1 if the source came from the label, 2 if it came from the old *_prev keys (old label), -1 if the
// restart is refused (message printed, _restart_refused set).
static int restart_source (void)
{
  ForcingSource s;
  int kind = label_read ("checkpoint.phase", &s);
  if (kind == 0)
    restart_refuse ("no checkpoint.phase (the dump's label) in the working directory");
  else if (kind < 0)
    restart_refuse ("checkpoint.phase in the working directory has no \"t\" (not a label)");
  else if (kind == 1 && params.has_prev_keys)
    restart_refuse ("params.json sets *_prev keys, but the label already gives the source forcing "
                    "(two sources); remove the *_prev keys");
  else if (kind == 2 && !params.has_prev_keys)
    restart_refuse ("checkpoint.phase is an old label (phase only, no forcing) and params.json has no "
                    "*_prev keys, so the source forcing is unknown");
  // The label must be the dump's own: a label left from another dump (another time) gives the wrong source.
  // params.t_checkpoint is the dump's exact time here (set_run_window reads it from the dump header); the
  // label prints its time with %.17g, so the two agree to round-off.
  else if (fabs (_label_t - params.t_checkpoint) > 1e-9*(fabs (_label_t) + fabs (params.t_checkpoint))) {
    static char why[256];
    snprintf (why, sizeof(why), "checkpoint.phase is the label of t=%.17g, but the dump is at t=%.17g "
              "(a stale label from another dump)", _label_t, params.t_checkpoint);
    restart_refuse (why);
  }
  if (_restart_refused)
    return -1;
  if (kind == 2) {
    // Old label: the source forcing comes from the *_prev keys, as before 2026-11-13. Allowed, but loud.
    if (pid() == 0)
      fprintf (stderr, "WARNING: old checkpoint label (no forcing); the ramp starts from the *_prev keys "
               "of params.json. Labels written from 2026-11-13 on carry the forcing.\n");
    s.omega_b = params.omega_b_prev;
    s.omega_h = params.omega_h_prev;
    for (int k = 0; k < N_MAX; k++) {
      s.theta_max[k]      = params.theta_max_prev[k];
      s.phi_angular[k]    = params.phi_angular_prev[k];
      s.amplitude_h[k]    = params.amplitude_h_prev[k];
      s.phi_horizontal[k] = params.phi_horizontal_prev[k];
    }
    s.cycles_h = 0.;
    s.U_bio = 0.;
    if (s.omega_b > 0.) {   // the old velocity unit, from the source rate and angle (set_scales() formula)
      double H = 2.*L_bio*Ly, Th0 = s.theta_max[0]*pi/180.;
      s.U_bio = L_bio/4*(H + 0.5*L_bio*tan(Th0))/(H*0.5)/(2.*pi/s.omega_b);
    }
  }
  fsrc = s;
  return kind;
}

// Fresh start: the bag starts at rest, so the amplitude ramps from zero at the target frequency.
static void fresh_source (void)
{
  ForcingSource s = {0};
  s.omega_b = w_bio;            // w_bio, not params.omega_b: bit-for-bit the target rate (no ramp in w)
  s.omega_h = params.omega_h;
  fsrc = s;
}

// restart_forcing.json: the forcing of THIS run at the restart time, to compare with the label (guard S13).
// Written by the acceleration event at its first step after a restart, from the forcing it has just computed:
// qcc's dimension check runs main() and the events and splits on every undefined condition, so each extra
// forcing_now() call costs much of its call budget (exit code 2, no message, 2026-11-13).
static void write_restart_forcing (ForcingNow F, double t_now, int kind)
{
  if (pid() != 0)
    return;
  FILE * fp = fopen ("restart_forcing.json", "w");
  if (!fp)
    return;
  fprintf (fp, "{\"t\": %.17g, \"T_bio\": %.17g, \"Th\": %.17g, \"Th_d\": %.17g, \"cycles\": %.17g, \"omega_b\": %.17g, "
           "\"omega_b_source\": %.17g, \"theta_max_source\": %.17g, \"U_bio_source\": %.17g, "
           "\"source\": \"%s\"}\n", t_now, T_bio, F.Th, F.Th_d/T_bio, F.psi/(2.*pi), F.omega_b, fsrc.omega_b,
           fsrc.theta_max[0], fsrc.U_bio, kind == 1 ? "label" : "old label + *_prev");
  fclose (fp);
}
