# =============================================================================
# One run = one copy of this file.
#
#   cp cases/template.sh cases/myrun.sh
#   nano cases/myrun.sh
#   ./submit_case.sh cases/myrun.sh
#
# submit_case.sh builds from these settings, submits the producer, then submits
# the offline solver with a dependency so it starts only if the producer
# succeeded.  A copy of this file is kept beside the build and beside the
# output, so every run records exactly how it was configured.
#
# This file is sourced by bash: plain `VAR=value`, no exports needed.
# =============================================================================

# ---- name ------------------------------------------------------------------
# The output lands in $OUTROOT/$RUNTAG and the build in oscar/builds/$RUNTAG.
# Leave empty to get N<N>_nmix<NMIX>_<timestamp>, or set a name of your own.
RUNTAG="pre_B_N256"

# ---- producer source -------------------------------------------------------
# codes_bio_run = the reference producer (default).  codes_bm = the same plus
# the biomass model; set it together with BR_BIOMASS=1 below.
PRODUCER_SRC=codes_bio_run

# ---- physical case ---------------------------------------------------------
ARGS="0.25 7 20"          # L_bio (m), rocking angle (deg), rpm

# ---- resolution and dump schedule ------------------------------------------
N=256
NDUMP=2000

# NMIX = rocking cycles of hydrodynamics simulated BEFORE dumping starts.
# This is the whole cost of the producer.  0 for a pipeline test; 5 is the
# value BioReactor.c ships with; 80 is what the stored N=512 dataset used.
NMIX=20

# NSNAP_MAX caps the dump window to this many dt_dump INTERVALS, so you get
# NSNAP_MAX+1 snapshots.  This is the only knob that shortens the run --
# t_end does not.
#   NSNAP_MAX=5      6 snapshots: pipeline/timing test only
#   NSNAP_MAX=99     ONE FULL ROCKING CYCLE for replay (NDUMP=100): phases
#                    0 .. 99/100 of the period, so the offline wrap from the
#                    last snapshot back to the first is one ordinary dt_dump
#                    step and the replay period is exactly 0.6074.  Use
#                    NDUMP-1 if you change NDUMP.
#   NSNAP_MAX=       (empty) NDUMP+1 snapshots: the last one repeats the first
#                    phase, so a replay holds that phase twice per cycle.
# For a replay to be periodic the flow must also BE periodic: set NMIX >= 1
# so the stored cycle is not the start-up transient from rest.
NSNAP_MAX=5999

# Final simulation time (simulation units; one rocking period = 0.6074).
# The run stops at the LATER of t_end and the end of the dump window, so:
#   T_END_PROD=        (empty) stop as soon as the last snapshot is written,
#                      i.e. (NMIX + NSNAP_MAX/NDUMP) x 0.6074
#   T_END_PROD=30      keep running the full simulation to t=30 after dumping
#                      (logstats, figures and scalar outputs continue)
# BioReactor.c's own value is 50 (~82 cycles); it is no longer used as-is.
T_END_PROD=14.0

# Interval of the Data_all field dumps (needs BR_OUT_FILES=1 below), in
# simulation units.  Empty = BioReactor.c's 0.1519*7 ~ 1.06 (~1.75 periods),
# too coarse for mixing curves.  0.06074 = 10 per period; each dump at N=512
# is ~25 MB of text per output time, so mind the total.
DT_FILE_PROD=0.06074

# =============================================================================
# BioReactor.c compile-time switches.
#
# Every BR_<NAME> below rewrites "#define <NAME>" in the staged copy of
# BioReactor.c.  Any switch in that file can be driven from here just by adding
# a BR_ line -- the build script looks them up by name and errors out if a name
# does not exist, so a typo is caught rather than silently ignored.
#
# Comment a line out to leave BioReactor.c's own default in place.
# =============================================================================

# ---- physics ---------------------------------------------------------------
BR_EMBED=1                # embedded boundary for the solid walls
BR_CONTACT=0              # contact-angle boundary condition
BR_ACCELERATION=1         # rocking motion

# ---- oxygen ----------------------------------------------------------------
BR_OXYGEN=1               # transport dissolved oxygen
BR_OXYGEN_AIR=1           # seed (and replenish) oxygen in the gas phase
BR_OXYGEN_CIRCLE=0        # alternative: seed oxygen in a circular blob

# ---- passive tracers -------------------------------------------------------
# NOTE: with BR_TRACER=0 the producer transports oxygen only, so the offline
# run's tracers have nothing to be compared against.  Turn these on to make
# the full run and the offline run solve the same scalars -- at the cost of
# five scalar solves per step in the producer instead of one.
BR_TRACER=1
BR_HORIZONTAL_MIXL=1      # tracer c  : left half
BR_HORIZONTAL_MIXR=0      # tracer c1 : right half
BR_VERTICAL_MIXUP=1       # tracer c2 : upper half
BR_VERTICAL_MIXDOWN=0     # tracer c3 : lower half

# ---- biomass ---------------------------------------------------------------
BR_BIOMASS=0              # needs PRODUCER_SRC=codes_bm (and BR_OXYGEN=1)
BR_BIOMASS_REPLE=0        # keep biomass out of the gas phase (codes_bm: must stay 0)
# Kinetics, codes_bm only (physical units; defaults in codes_bm/BioReactor.c).
# Uncomment to override:
#BR_BM_MU_MAX=0.04        # max specific growth rate [1/h]
#BR_BM_K_O2=0.1           # O2 half-saturation [mg/L]
#BR_BM_Q_O2=0.3           # max specific O2 uptake [pmol/cell/h]
#BR_BM_X0=1e6             # initial biomass [cells/mL]
#BR_BM_C_SAT=8.0          # liquid O2 saturation [mg/L]

# ---- mesh and numerics -----------------------------------------------------
# BR_AMR=1 makes the mesh adaptive, which breaks the offline import: the
# stored cells no longer form a uniform N x N grid, so init_grid(N) in the
# offline solver cannot reproduce them and cells go unmatched.  Leave at 0
# unless you are not importing the result.
BR_AMR=0
BR_REMOVE_DROP=0          # automatic droplet removal
BR_CFL_COND=0             # custom CFL number

# ---- producer output -------------------------------------------------------
BR_OUT_DUMP_BINARY=1      # 1 = binary flow fields (what the offline run reads)
BR_OUT_DUMP_DATA=0
BR_OUT_DUMP=0             # Basilisk dump() restart files
BR_OUT_FILES=0
BR_OUT_SPECIFIC_TIME=0    # high-frequency output over a focused window
BR_OUT_INTERFACE=0        # interface facets
BR_NORMCAL=1              # norm / statistics output
BR_FIGURES=0              # PNGs      (off: unrelated to the import, and slow)
BR_VIDEOS=0               # MP4s      (off: same)

# ---- where things go -------------------------------------------------------
OUTROOT="/oscar/scratch/eaguerov/kim_pipeline/prereg_out"
BASILISK="/oscar/data/dharri15/eaguerov/basilisk/src"
TPP="/oscar/scratch/eaguerov/kim_pipeline/TPP_headers"

# ---- producer job (A2) -----------------------------------------------------
# Sizing reference: N=256 on 4 ranks ran at 1.68 it/s with dt~2.8e-4 on a
# laptop, i.e. ~2100 iterations and ~22 min per rocking cycle.  N=512 is 4x the
# cells with about half the timestep, so ~8x that at the same rank count.
P_NODES=1
P_NTASKS=8
P_MEM=16G
P_TIME=24:00:00
P_PARTITION=batch

# ---- offline job (A3) ------------------------------------------------------
# AT LEAST 2.  The scalar solver fails to converge on a SINGLE MPI rank
# (see OPTIMISATION.md section F) -- never set this to 1.
O_NODES=1
O_NTASKS=16
O_MEM=32G
O_TIME=2:00:00
O_PARTITION=batch

# Set to 0 to submit the producer only (e.g. to inspect the dumps first).
RUN_OFFLINE=0
BR_PREREG_REC=0
BR_PREREG_DIAG=1
