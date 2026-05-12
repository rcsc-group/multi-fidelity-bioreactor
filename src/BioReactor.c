// Primary author: Minki Kim
// Contributor: Radu Cimpeanu
// Contributor: Daniel Harris
// Date: 03/04/2025

// Core Basilisk modules for embedding geometries and centered Navier-Stokes solver
#include "embed.h"
#include "navier-stokes/centered.h"

// Define dynamic viscosity as a function of the volume fraction (f)
// Uses harmonic averaging for two-phase flow
#define mu(f) ( 1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)  )

// Modules for two-phase flow (e.g., air and water) with surface tension and volume conservation
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"

// Modules for oxygen transfer + tracers
#include "henry_oxy2.h"

// Custom headers foroxy specific configurations and visualization
#include "view3.h"
#include "utils2.h"

// Mathematical constants (e.g., M_PI from math.h)
#define _USE_MATH_DEFINES
#include <math.h>

// JSON parameter reader (jsmn-based); replaces argv[1..3] with params.json
#include "params_read.h"

// Flags to control the inclusion of features (set to 1 = enable, 0 = disable)
#define EMBED            1   // Enable embedded boundary for solid geometry
#define CONTACT          0   // Enable contact angle boundary condition
#define OXYGEN           1   // Enable oxygen concentration simulation
#define OXYGEN_CIRCLE    0   // Initial distribution (circle) of oxygen (if OXYGEN == 1)
#define OXYGEN_AIR       1   // Initial distribution (air side) of oxygen (if OXYGEN == 1)

// Mixing strategies for tracer release
#define TRACER           1   // Enable passive tracer simulation
#define HORIZONTAL_MIXL  0   // Initial distribution (left side) of tracer: Horizontal mixing (if TRACER == 1)
#define HORIZONTAL_MIXR  0   // Initial distribution (right side) of tracer: Horizontal mixing (if TRACER == 1)
#define VERTICAL_MIXUP   1   // Initial distribution (top side) of tracer: Vertical mixing (if TRACER == 1)
#define VERTICAL_MIXDOWN 0   // Initial distribution (bottom side) of tracer: Vertical mixing (if TRACER == 1)

// Other simulation options
#define ACCELERATION     1   // Enable acceleration (rocking motion)
#define AMR              0   // Enable adaptive mesh refinement
#define REMOVE_DROP      0   // Enable automatic droplet removal
#define CFL_COND         0   // Use custom CFL number
#define DUMP             0   // Save dump output
#define NORMCAL          1   // Calculate statistics (norms)
#define FIGURES          0   // Figures: enable only for diagnostics; output dirs created by simulate.py
#define VIDEOS           0   // Videos: enable only for diagnostics; not needed in optimization loop

// Output options
#define OUT_FILES         0   // Full-field dumps: enable only for diagnostics; output dir created by simulate.py
#define OUT_SPECIFIC_TIME 0   // Output data at specific time ranges
#define OUT_INTERFACE     1   // Save interface geometry



// ================================================================== //
//                       SIMULATION SETUP                             //
// ================================================================== //
int NN;  // Grid resolution: set from params.fidelity as 1<<fidelity in main()
#define N_RAMP_CYCLES 3            // Ramp duration in rocking cycles; t_change_st = N_RAMP_CYCLES * T_per_st
const double th_cont = 90;        // Contact angle for wetting conditions (degrees)
double t_mix,t_dump;              // Time at which tracer is released, and dump file is saved (computed later)
const double nMix_cycle = 80;     // Number of cycles for tracer release (used to compute t_mix)
double t_end;                      // Final simulation time (simulation time unit); set from params.t_end in main()

// Output time intervals (derived from experimental timing)
const double dt_file = 0.1519*7;  // Interval for saving data to file
const double dt_video= 0.6074/10; // Interval for video frames
const double dt_Fig  = 0.1519*7;  // Interval for figure output
double t_spec_init, t_spec_end;   // Specific output window for focused data extraction
const double dt_spec = 0.000530525;  // Very high frequency sampling for specific data

const int    i_fig   = 5000;      // Output interval for figures
const double t_out   = 0.1;       // Output interval for statistics [non-dim time]
const double CFL_num = 0.01	;     // CFL number for time-step stability (used only if CFL_COND is enabled)
const double N_output= 128;       // Resolution for output file if needed

// Parameters for drop/bubble removal (if enabled)
const double remove_minsize   = 20;      // the Minimum number of grids for removal
const double remove_threshold = 1.0e-4;  // Threshold for identifying disconnected fluid elements


// ================================================================== //
//                       REACTOR GEOMETRY                             //
// ================================================================== //
double LL = 1.0;        // Width of the reactor domain (non-dimensionalized by L_bio, so always 1)
double Ly = 0.286;      // Height of the reactor domain (non-dimensionalized by L_bio)
double y_init = 0.0;    // Initial liquid height (defines volume fraction f = 1 below this height)
double L_piv  = 0.143;  // Distance from center to pivot point (rocking axis), affects acceleration computation
double L_bio  = 0.25;   // Characteristic length (bag half-length, m); global so events can non-dimensionalize with it


// ================================================================== //
//                     MATERIAL PROPERTIES                            //
// ================================================================== //
const double rho_w = 1.0e3;   // Water density (kg/m^3)
const double rho_a = 1.225;   // Air density at 20°C (kg/m^3)
const double mu_w  = 1.0e-3;  // Water viscosity at 20°C (Pa·s)
const double mu_a  = 1.81e-5; // Air viscosity at 20°C (Pa·s)
const double grav  = 9.8;     // Gravitational acceleration (m/s^2)
const double sigma = 0.0728;  // Surface tension at the liquid-gas interface (N/m)

// Note: Phase 1 = liquid, Phase 0 = air
const double D_tracer_1 = 0.44e-9; // Diffusivity of tracers in water (m²/s)
const double D_tracer_2 = 1.0e-30; // Diffusivity of tracers in air (m²/s); near-zero in air
const double D_oxy_1 = 1.90e-9;    // Diffusivity of oxygen in water (m²/s)
const double D_oxy_2 = 1.98e-5;    // Diffusivity of oxygen in air (m²/s)
const double c_tracer_alpha = 1.0e30;    // Tracer solubility: c_water = c_tracer_alpha*c_air
const double c_oxy_alpha    = 1./30;     // Oxygen solubility: oxy_water = c_oxy_alpha*oxy_air (=1/30 tpiycal for O2 in water)

// Change contact angle
#if CONTACT
vector h[];

// Apply static contact angle boundary conditions (in radians)
h.t[left]  = contact_angle(th_cont*pi/180);  // Convert degrees to radians
h.t[right] = contact_angle(th_cont*pi/180);
#endif

// ================================================================== //
//                 SCALARS AND FIELD POINTERS                         //
// ================================================================== //
// Parameters related to oscillation and tracer release
double Th_max, T_per, R_tr, x_tr, y_tr;

// Scalars to track tracer concentrations and oxygen
scalar c[], oxy[], c1[], c2[], c3[];   // for tracer and oxygen transfer
scalar * stracers = {c,oxy,c1,c2,c3};
double (* gradient) (double, double, double) = minmod2;   // Custom slope limiter function (used for scalar gradients)

// Buffers for file naming and output file pointers for statistics
char buf1[100], buf2[100], buf3[100], buf4[100];
FILE * fp_stats, * fp_norm, * fp_stats2, * fp_stats3;

// Key physical and dimensionless parameters (computed in main)
double U0, Ub, Re_w, Re_a, We_w, Fr, rhor, mur, Pe_tracer_1, Pe_tracer_2, Pe_oxy_1, Pe_oxy_2, Th, Th_d, Th_2d, U_bio, w_bio, w_bio_st, T_per_st, T_bio, Th_max2, D_in_non, U_in_non, t_change_st, t_mix_st;
BioreactorParams params;  // global so all events (acceleration, init) can access it
int MINLEVEL, MAXLEVEL;    // Mesh refinement levels


// ================================================================== //
//                       MAIN FUNCTION                                //
// ================================================================== //
int main(int argc, char * argv[]){

  if (argc < 2) { fprintf(stderr, "Usage: BioReactor params.json\n"); return 1; }
  params = params_read(argv[1]);

  // Derive legacy scalars from params for the rest of main() unchanged.
  // omega_b (rad/s) and theta_max[0] (deg) replace the old ANGLE/RPM CLI args.
  // Multi-harmonic forcing and superellipse geometry will consume params directly
  // once the acceleration event and init event are extended (next steps).
  L_bio = params.geometry_a;  // bag half-length (m); sets all non-dim scales
  double ANGLE = params.theta_max[0];       // Fundamental rocking amplitude (degrees)
  double RPM   = params.omega_b * 60. / (2.*M_PI);  // Convert rad/s → RPM
  NN    = 1 << params.fidelity;             // fidelity → grid cells per side (4→16, 7→128, 9→512)
  t_end = params.t_end;                    // non-dimensional sim time; 1 unit = T_bio seconds

  L0 = 1. [0];  // [0] declares space dimensionless: simulation is fully non-dimensionalized (scaled by L_bio, T_bio, U_bio); Basilisk dimensional analysis requires annotations on literals, not variables
  DT = HUGE [0];  // [0] declares time dimensionless; older qcc omitted the dimensional(u.x[]==Delta/DT) constraint from centered.h that now makes this required
  origin(-L0/2., -L0/2.);   // Set coordinate origin to domain center

#if !AMR
  init_grid(NN);            // Initialize uniform grid if AMR is disabled
#endif

#if AMR
  MAXLEVEL = params.fidelity;
  MINLEVEL = params.fidelity - 2;
  double F_MAX = 1e-6;    // Refinement threshold for volume fraction
  double U_MAX = 0;       // Refinement threshold for velocity refinement
#endif

  // Bag geometry from params (must precede anything that uses Ly or H_bio)
  Ly = params.geometry_b / L_bio;  // dimensionless half-height; overwrites hardcoded 0.286

  // Rocking motion parameters
  R_tr  = 0.0084/L_bio;    // Tracer radius (scaled with domain)
  x_tr  = 0.;              // Tracer x-position (center)
  y_tr  = -Ly*0.5*0.5;     // y-position: halfway into the bottom half of domain
  Th_max = ANGLE*pi/180;   // Max angle in radians
  T_per  = 60./RPM;        // Period in seconds
  
  // Characteristic scales
  double H_bio,V_bio;
  H_bio  = L_bio*Ly;
  V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
  U_bio  = V_bio/(H_bio*0.5)/T_per; // Characteristic velocity scale
  T_bio  = L_bio/U_bio;             // Characteristic time scale
  w_bio  = 2*pi/T_per;              // Angular velocity (rad/s)
  w_bio_st = w_bio*T_bio;           // Dimensionless angular velocity
  T_per_st = T_per/T_bio;           // Dimensionless period
  U0     = w_bio_st*Th_max;         // Initial rotational velocity
  t_change_st = N_RAMP_CYCLES * T_per_st;  // ramp over 3 cycles regardless of omega_b
  t_mix      = T_per_st*nMix_cycle; // Time to release tracers (simulation time)
  t_dump = t_mix;                   // Time to dump data (simulation time)
  t_spec_init= t_mix;                       // Initial time to save data in the specific range
  t_spec_end = T_per_st*(nMix_cycle+10);    // End time to save data in the specific range

  // Dimensionless numbers
  Re_w = rho_w*U_bio*L_bio/mu_w;            // Reynolds number of water
  Re_a = rho_a*U_bio*L_bio/mu_a;            // Reynolds number of air
  We_w = rho_w*U_bio*U_bio*L_bio/sigma;     // Weber number of water
  Fr   = U_bio/sqrt(grav*L_bio);            // Froude number
  rhor = rho_a/rho_w;                       // Density ratio
  mur  = mu_a/mu_w;                         // Viscosity ratio
  Pe_tracer_1 = U_bio*L_bio/D_tracer_1;     // Peclet number of tracer at phase 1 (water)
  Pe_tracer_2 = U_bio*L_bio/D_tracer_2;     // Peclet number of tracer at phase 2 (air)
  Pe_oxy_1    = U_bio*L_bio/D_oxy_1;        // Peclet number of oxygen at phase 1 (water)
  Pe_oxy_2    = U_bio*L_bio/D_oxy_2;        // Peclet number of oxygen at phase 1 (air) 

  rho1 = 1.0;             // Reference density of phase 1 (water) (scaled to 1)
  rho2 = rho1*rhor;       // Scaled density of phase 2 (air)
  mu1  = 1.0/Re_w;        // Dimensionless viscosity of phase1 (water)
  mu2  = mur*mu1;         // Scaled viscosity of phase 2 (air)
  f.sigma = 1.0 / We_w;   // Dimensionless surface tension
  
  // Contact angle setup (if enabled)
  #if CONTACT
    f.height = h;
  #endif

  // Tracer field configuration
  #if TRACER  
  // Define diffusivities in each phase (1 = liquid, 0 = air)
    c.D1 = 1./Pe_tracer_1;
    c.D2 = 1./Pe_tracer_2;
  
  // Tracer concentration ratio across interface: c_liquid = alpha*c_air; c1=alpha*c2
    c.alpha = c_tracer_alpha;

  // Apply slope limiter to gradients for numerical stability
    c.gradient = minmod2;

  // --- Optional extra tracers for comparisons --- //
    c1.D1 = 1./Pe_tracer_1;
    c1.D2 = 1./Pe_tracer_2;
    c1.alpha = c_tracer_alpha;
    c1.gradient = minmod2;
    c2.D1 = 1./Pe_tracer_1;
    c2.D2 = 1./Pe_tracer_2;
    c2.alpha = c_tracer_alpha;
    c2.gradient = minmod2;
    c3.D1 = 1./Pe_tracer_1;
    c3.D2 = 1./Pe_tracer_2;
    c3.alpha = c_tracer_alpha;
    c3.gradient = minmod2;
#endif

// Oxygen field configuration
#if OXYGEN
  oxy.D1 = 1./Pe_oxy_1;     // Oxygen diffusivities (water) 
  oxy.D2 = 1./Pe_oxy_2;     // Oxygen diffusivities (air) 
  oxy.alpha = c_oxy_alpha;  // Henry's law partition coefficient

  // Apply slope limiter for smooth gradients
  oxy.gradient = minmod2;
#endif

// Boundary conditions
  u.n[left]  = dirichlet(0.);  // Set no-slip velocity boundary conditions (zero velocity)
  u.t[left]  = dirichlet(0.);
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

// Solid geometry Boundary (if using embedded solids)
#if EMBED
  u.n[embed] = dirichlet(0.); // Set no-slip velocity boundary conditions (zero velocity)
  u.t[embed] = dirichlet(0.);
#endif

// CFL control (if enabled)
#if CFL_COND
  CFL = CFL_num;
#endif
  

// Output Files
  char name[200],name2[200],name3[200],name4[200];

  sprintf(name, "logstats.dat");        // Performance & time log
  sprintf(name2,"normf.dat");           // Norms and statistics of fields
  sprintf(name3,"vol_frac_interf.dat"); // Volume fraction & interface stats
  sprintf(name4,"tr_oxy.dat");          // Tracer and oxygen stats
  fp_stats = fopen(name, "w");   
  fp_norm  = fopen(name2,"w");
  fp_stats2= fopen(name3,"w");
  fp_stats3= fopen(name4,"w");

  fprintf(fp_norm, "i t Omega_liq_avg Omega_liq_rms Omega_liq_vol Omega_liq_max ux_liq_avg ux_liq_rms ux_liq_vol ux_liq_max uy_liq_avg uy_liq_rms uy_liq_vol uy_liq_max \n");
  fprintf(fp_stats2, "i t f_liq_sum f_liq_interf posY_max posY_min \n");
  fprintf(fp_stats3, "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 c1_liq_sum c1_liq_sum2 c2_liq_sum c2_liq_sum2 c3_liq_sum c3_liq_sum2 \n");

  NITERMAX = 1000;     // Max iterations per timestep
  TOLERANCE = 5.0e-4;  // // Solver tolerance (convergence criterion)
  
  // Run the simulation — Basilisk manages events from here
  run();
  
  // Close all output files
  fclose(fp_stats); fclose(fp_norm); fclose(fp_stats2); fclose(fp_stats3);
}


// ================================================================== //
//                      INITIAL CONDITIONS                            //
// ================================================================== //
event init (t = 0)
{
  // Parametric bag geometry (dimensionless semi-axes)
  double a_nd = params.geometry_a / L_bio;
  double b_nd = params.geometry_b / L_bio;  // == Ly

  // Fill level: liquid occupies fill_level fraction of bag height, measured from bottom
  double y_fill = b_nd * (2.*params.fill_level - 1.);
  fraction(f, y_fill - y);

  // Superellipse solid: |x/a|^n + |y/b|^n = 1
  // n >= 8 → perfect rectangle (avoids pow() singularities at sharp corners)
  #if EMBED
  if (params.geometry_n >= 8.)
    solid(cs, fs, intersection(a_nd - fabs(x), b_nd - fabs(y)));
  else
    solid(cs, fs, 1. - pow(fabs(x/a_nd), params.geometry_n)
                     - pow(fabs(y/b_nd), params.geometry_n));
  #endif
}


// ================================================================== //
//                           TRACER SETUP                             //
// ================================================================== //
#if TRACER
event tracer(t = t_mix){

  double h_tr;
  h_tr = (M_PI*R_tr*R_tr);  // Area of circular tracer patch

  // circular shape tracer at the center of the liquid
  // fraction(c, -(sq(x-x_tr) + sq(y-y_tr) - sq(R_tr)) ); 

  // tracer released as a line (same area)
  // fraction(c, intersection( -(y-y_tr - 0.5*h_tr), -(-(y-y_tr + 0.5*h_tr)) ));

  // Horizontal mixing-left side
  #if HORIZONTAL_MIXL
  fraction (c,0-x);  // Fill left half of domain
  foreach(){
    if ((cs[] == 0) || (f[] < 1))
      c[] = 0;        // Don't place tracer in gas or outside solid
  }
  #endif

  // Horizontal mixing-right side
  #if HORIZONTAL_MIXR
  fraction (c1,0+x);  // Fill right half of domain
  foreach(){
    if ((cs[] == 0) || (f[] < 1))
      c1[] = 0;       // Don't place tracer in gas or outside solid
  }
  #endif

  // Vertical mixing-top side
  #if VERTICAL_MIXUP
  foreach(){
    if ((f[] == 1) && (cs[]==1) && (y >= -Ly*0.5*0.5))
      c2[] = 1.0;    // Upper half of liquid
  }
  #endif

  // Vertical mixing-down side
  #if VERTICAL_MIXDOWN
  foreach(){
    if ((f[] == 1) && (cs[]==1) && (y <= -Ly*0.5*0.5))
      c3[] = 1.0;    // Lower half of liquid
  }
  #endif
}
#endif


// ================================================================== //
//                           OXYGEN SETUP                             //
// ================================================================== //
#if OXYGEN
event oxygen (t=t_mix; i++){

#if OXYGEN_CIRCLE
  fraction(oxy, -(sq(x-0) + sq(y-Ly*0.5*0.5) - sq(0.084*Ly)) );
#endif

#if OXYGEN_AIR
  foreach(){
    if ((f[] == 0) && (cs[]==1))
      oxy[] = 1.;     // Oxygen in gas regions only
  }
#endif
}
#endif


// ================================================================== //
//                           ACCELERATION                             //
// ================================================================== //
#if ACCELERATION
event acceleration(i++)
{
  // Linear amplitude ramp over t_change_st to avoid impulsive start
  double ramp = (t < t_change_st) ? t / t_change_st : 1.0;

  // Multi-harmonic angular forcing: θ(t) = Σ_k θ_k · sin(k·ω_b·t + φ_k)
  // k is harmonic number (dimensionless); k·w_bio_st is the k-th frequency in sim time
  Th = 0;  Th_d = 0;  Th_2d = 0;
  for (int k = 1; k <= params.n_harmonics; k++) {
    double wk  = k * w_bio_st;
    double Ak  = params.theta_max[k-1] * pi / 180.;  // deg → rad
    double phk = params.phi_angular[k-1];
    Th    +=  Ak * sin(wk*t + phk);
    Th_d  +=  Ak * wk * cos(wk*t + phk);
    Th_2d += -Ak * wk*wk * sin(wk*t + phk);
  }
  Th *= ramp;  Th_d *= ramp;  Th_2d *= ramp;

  // Multi-harmonic horizontal forcing: x(t) = Σ_k A_k · sin(k·ω_h·t + ψ_k)
  // x_acc is the non-inertial horizontal body force (dimensionless: A/L_bio · (k·w_h_st)²)
  double x_acc = 0.;
  if (params.omega_h > 0.) {
    double w_h_st = params.omega_h * T_bio;  // ω_h non-dimensionalized
    for (int k = 1; k <= params.n_harmonics; k++) {
      double wk_h = k * w_h_st;
      x_acc += (params.amplitude_h[k-1] / L_bio) * sq(wk_h)
               * sin(wk_h*t + params.phi_horizontal[k-1]);
    }
    x_acc *= ramp;
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


// ================================================================== //
//                           OTHER OPTIONS                            //
// ================================================================== //

// Droplet and bubble removal
#if REMOVE_DROP
event remove_drop(i++){
  // Remove disconnected small fluid regions (droplets/bubbles)
  remove_droplets(f,remove_minsize,remove_threshold,false);  // remove droplets
  remove_droplets(f,remove_minsize,remove_threshold,true);   // remove bubbles
}
#endif

// Cumstom CFL Condition
#if CFL_COND
event CFL_cond(i++){
  CFL = CFL_num;   // Force custom CFL number (defined earlier)
}
#endif

// Adaptive Mesh Refinement (AMR)
#if AMR
event adapt( t=0 ){  
  // Refine inside the bioreactor domain only
  refine(level<MAXLEVEL && (  (y > -0.7*Ly) && (y < 0.7*Ly) ));
}
#endif

// Dump raw wata to files
#if DUMP
event dump(t=t_dump){
  dump(file="dump");    // Save entire simulation state

  snprintf(buf1, sizeof(buf1), "Dump_%d_%g_%d.txt",N,t,pid());
  FILE * out_all = fopen(buf1,"w");
  foreach(){
    fprintf(out_all,"%g %g %g %g %g %g \n",x,y,u.x[],u.y[],f[],c[]);
  }
  fclose(out_all);
}
#endif

//  Log performance and runtime
event logstats (t+=0.1; t <= t_end) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    if (pid() == 0){
      fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
      fflush(fp_stats);
    }
}

// Compute simulation statistics
#if NORMCAL
event normcal (t+=t_out; t<=t_end){
    //timing s = timer_timing (perf.gt, i, perf.tnc, NULL);

    scalar ux_liq[],uy_liq[],ux_liq_abs[],omega[],omega_liq[],oxy_liq[],f_liq[],posY[],c_liq[],c1_liq[],c2_liq[],c3_liq[];
    double omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max;
    double ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max;
    double f_liq_sum,f_liq_interf,posY_max,posY_min,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2;
    
    vorticity (u, omega); // vorticity

    // only liquid velocity
    foreach(){
      ux_liq[]  = u.x[]*f[];
      uy_liq[]  = u.y[]*f[];
      oxy_liq[] = oxy[]*f[];
      omega_liq[] = omega[]*f[];
      f_liq[]   =  (1-cs[])*f[];
      c_liq[]   = c[]*f[];
      ///*
      c1_liq[]  = c1[]*f[];
      c2_liq[]  = c2[]*f[];
      c3_liq[]  = c3[]*f[];
      //*/
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
    posY_max      = statsf(posY).max;
    posY_min      = statsf(posY).min;

    oxy_liq_sum   = statsf2(oxy_liq).sum;
    oxy_liq_sum2  = statsf2(oxy_liq).sum2;
    c_liq_sum     = statsf2(c_liq).sum;
    c_liq_sum2    = statsf2(c_liq).sum2;
    ///*
    c1_liq_sum    = statsf(c1_liq).sum;
    c1_liq_sum2   = statsf2(c1_liq).sum2;
    c2_liq_sum    = statsf(c2_liq).sum;
    c2_liq_sum2   = statsf2(c2_liq).sum2;
    c3_liq_sum    = statsf(c3_liq).sum;
    c3_liq_sum2   = statsf2(c3_liq).sum2;
    //*/

   // i, timestep, no of cells, real time elapsed, cpu time
   if (pid() == 0){
      
      fprintf(fp_norm, "%i %g %g %g %g %g %g %g %g %g %g %g %g %g \n",i,t,omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max,ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max);
      fflush(fp_norm);

      fprintf(fp_stats2, "%i %g %g %g %g %g \n",i,t,f_liq_sum,f_liq_interf,posY_max,posY_min);
      fflush(fp_stats2);

      fprintf(fp_stats3, "%i %g %g %g %g %g %g %g %g %g %g %g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2);
      //fprintf(fp_stats3, "%i %g %g %g %g %g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2);
      fflush(fp_stats3);
   }
}
#endif

// Make videos
#if VIDEOS
event movies_output(t = t_mix; t += dt_video; t<=t_end)
{
  scalar omega[], ff[], cc[], oxyy[];
  char timestring[100];

  vorticity (u,omega);

  // vorticity
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);    // outlines interfaces
  squares("omega",map=cool_warm,min=-50.0,max=50.0);  // maps scalar fields to color
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);   // overlays timestamp
  save("vorticity3.mp4");
  save("vorticity.png");

  // volume fraction
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("f",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  //cells();
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("volume_fraction3.mp4");
  save("volume_fraction.png");

#if TRACER
  // tracer
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("tracer.mp4");
  save("tracer.png");

  // tracer1
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c1",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("tracer1.mp4");
  save("tracer1.png");

  // tracer2
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c2",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("tracer2.mp4");
  save("tracer2.png");

  // tracer3
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c3",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("tracer3.mp4");
  save("tracer3.png");
#endif
// oxygen
#if OXYGEN
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("oxy",map=cool_warm,min=0.0,max=0.033);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save("oxygen3.mp4");
  save("oxygen.png");
#endif
}
#endif

#if FIGURES
event Figures(t=t_mix; t<=t_end; t += dt_Fig)
{
  scalar omega[];
  char timestring[100],figN1[100],figN2[100],figN3[100],figN4[100],figN5[100],figN6[100],figN7[100];
  
  vorticity (u,omega);
  
  snprintf(figN1, sizeof(figN1), "Fig_vor/vor_%d_%.12g.png",N,t);
  snprintf(figN2, sizeof(figN2), "Fig_vol/vol_%d_%.12g.png",N,t);
  snprintf(figN3, sizeof(figN3), "Fig_tr/tr_%d_%.12g.png",N,t);
  snprintf(figN4, sizeof(figN4), "Fig_tr/tr1_%d_%.12g.png",N,t);
  snprintf(figN5, sizeof(figN5), "Fig_tr/tr2_%d_%.12g.png",N,t);
  snprintf(figN6, sizeof(figN6), "Fig_tr/tr3_%d_%.12g.png",N,t);
  snprintf(figN7, sizeof(figN7), "Fig_oxy/oxy_%d_%.12g.png",N,t);

  // vorticity
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("omega",map=cool_warm,min=-50.0,max=50.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN1);

  // volume fraction
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("f",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  //cells();
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN2);

  // tracer
#if TRACER
  #if HORIZONTAL_MIXL
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN3);
  #endif

  #if HORIZONTAL_MIXR
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c1",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN4);
  #endif

  #if VERTICAL_MIXUP
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c2",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN5);
  #endif

  #if VERTICAL_MIXDOWN
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("c3",map=cool_warm,min=0.0,max=1.0);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN6);
  #endif
#endif

  // oxygen
  #if OXYGEN
  clear();
  view(width=1200,height=1200,fov=24.0,ty=0.0);
  draw_vof("f",lw=2);
  squares("oxy",map=cool_warm,min=0.0,max=0.033);
  draw_vof("cs","fs");
  sprintf(timestring,"t=%2.03fs",t*T_bio);
  draw_string(timestring,pos=4,lc={0,0,0},lw=2);
  save(figN7);
  #endif
}
#endif

// Export field data
#if OUT_FILES
event out_files(t=t_mix; t<=t_end; t+=dt_file)
{
  scalar omega[];
  vorticity(u,omega);

  snprintf(buf1, sizeof(buf1), "Data_all/Data_all_%d_%.12g_%d.txt",N,t,pid());
  FILE * out_all = fopen(buf1,"wb");  

  fprintf(out_all,"x y ux uy vol_frac tracer solid oxygen vorticity tracer1-3 \n");
  foreach()
    fprintf(out_all,"%g %g %g %g %g %g %g %g %g %g %g %g\n",x,y,u.x[],u.y[],f[],c[],cs[],oxy[],omega[],c1[],c2[],c3[]);
  fclose(out_all);

  #if OUT_INTERFACE
    snprintf(buf4, sizeof(buf4), "Data_all/Interf_%d_%.12g_%d.txt",N,t,pid());
    FILE * out_interf = fopen(buf4,"wb");
    output_facets(f,out_interf);   // Interface extraction
    fclose(out_interf);
  #endif
}

// Export field data until t_mix from t =0
event out_files_initial(t=0; t<=t_mix; t+=dt_file)
{
  scalar omega[];
  vorticity(u,omega);

  snprintf(buf1, sizeof(buf1), "Data_all/Data_all_%d_%.12g_%d.txt",N,t,pid());
  FILE * out_all = fopen(buf1,"wb");  

  fprintf(out_all,"x y ux uy vol_frac tracer solid oxygen vorticity tracer1-3 \n");
  foreach()
    fprintf(out_all,"%g %g %g %g %g %g %g %g %g %g %g %g\n",x,y,u.x[],u.y[],f[],c[],cs[],oxy[],omega[],c1[],c2[],c3[]);
  fclose(out_all);

  #if OUT_INTERFACE
    snprintf(buf4, sizeof(buf4), "Data_all/Interf_%d_%.12g_%d.txt",N,t,pid());
    FILE * out_interf = fopen(buf4,"wb");
    output_facets(f,out_interf);   // Interface extraction
    fclose(out_interf);
  #endif
}
#endif

// Export field data during specific time ranges with higher sampling rates
#if OUT_SPECIFIC_TIME
event out_spec_time(t=t_spec_init; t<=t_spec_end; t+=dt_spec)
{
  scalar omega[];
  vorticity(u,omega);

  snprintf(buf2, sizeof(buf2), "Data_specific/Data_all_%d_%.12g_%d.txt",N,t,pid());
  FILE * out_all_spec = fopen(buf2,"wb");  

  fprintf(out_all_spec,"x y ux uy vol_frac tracer solid oxygen vorticity tracer1-3 \n");
  foreach()
  fprintf(out_all_spec,"%g %g %g %g %g %g %g %g %g %g %g %g \n",x,y,u.x[],u.y[],f[],c[],cs[],oxy[],omega[],c1[],c2[],c3[]);
  fclose(out_all_spec);

  #if OUT_INTERFACE
    snprintf(buf3, sizeof(buf3), "Data_specific/Interf_%d_%.12g_%d.txt",N,t,pid());
    FILE * out_interf_spec = fopen(buf3,"wb");
    output_facets(f,out_interf_spec);   // Interface extraction
    fclose(out_interf_spec);
  #endif
}
#endif
