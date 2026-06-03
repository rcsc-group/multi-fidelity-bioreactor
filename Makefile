# ────────────────────────────────────────────────────────────────────────────────
#  Rocking Bioreactor 2D  ◇  Unified Makefile
# ────────────────────────────────────────────────────────────────────────────────
#
#  QUICK START ────────────────────────────────────────────────────────────────
#    $ make                    # build the simulation executable
#    $ make run PARAMS=runs/my_run/params.json   # run one configuration
#
#  COMMON TASKS ────────────────────────────────────────────────────────────────
#  1) Build the simulation binary
#       $ make build
#
#     with custom compile-time flags:
#       $ make build LEVEL=7 DEBUG=1
#
#  2) Build and run with a params.json
#       $ make run PARAMS=runs/my_run/params.json
#
#  3) Submit a SLURM job for one params.json
#       $ make submit PARAMS=runs/my_run/params.json
#       $ make submit PARAMS=runs/my_run/params.json DRYRUN=1   # dry-run only
#
#  HOUSE-KEEPING ───────────────────────────────────────────────────────────────
#       $ make clean        # remove build artifacts
#       $ make deepclean    # nuke build/, patched headers, logs
#
#  OVERRIDABLE VARIABLES ───────────────────────────────────────────────────────
#       LEVEL       mesh refinement level (default 8; use 4–5 for quick tests)
#       DEBUG=1     add -g -O0
#       DRYRUN=1    print sbatch command without submitting
#       PARAMS      path to params.json (required for 'run' and 'submit')
#
#  ENVIRONMENT ─────────────────────────────────────────────────────────────────
#       Uses qcc from $BASILISK/.. or ~/scratch/basilisk/src/qcc as fallback.
#       Set $BASILISK before running, or let the Makefile find the scratch build.
#
#  Maintainer: Elvis Aguero    │   Last update: <2026-05-04>
# ────────────────────────────────────────────────────────────────────────────────


# ==========================================================
#  Locate qcc
#  Prefer the scratch-built qcc over any spack module because
#  the spack-installed qcc on OSCAR has a hardcoded dead build-
#  time path (/tmp/yliu385/spack-stage/…) and cannot find its
#  own headers.  The scratch build uses the real $BASILISK path.
# ==========================================================
SCRATCH_QCC := $(HOME)/scratch/basilisk/src/qcc
ifneq ($(wildcard $(SCRATCH_QCC)),)
  QCC := $(SCRATCH_QCC)
else
  QCC := $(shell command -v qcc 2>/dev/null)
  ifeq ($(QCC),)
    $(error "qcc not found. Build from ~/scratch/basilisk/src or add to PATH.")
  endif
endif

# Derive $BASILISK from qcc location if not already set
ifeq ($(BASILISK),)
  BASILISK := $(dir $(QCC))
  export BASILISK
endif


# ==========================================================
#  Compiler flags
# ==========================================================
DEBUG ?= 0
ifeq ($(DEBUG),1)
  CFLAGS = -g -O0
else
  CFLAGS = -O2
endif
CFLAGS += -w -fopenmp -Wall -autolink -lm

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  OPENGLIBS = -lfb_tiny -framework OpenGL
else
  OPENGLIBS = -lfb_tiny -lGL
endif
LDFLAGS = -L$(BASILISK)/gl -lglutils $(OPENGLIBS)


# ==========================================================
#  Paths
# ==========================================================
SRC_DIR    = src
BUILD_DIR  = build

SIM_SRC    := $(SRC_DIR)/BioReactor.c
EXECUTABLE := $(BUILD_DIR)/BioReactor

SRC_HEADERS := $(wildcard $(SRC_DIR)/*.h)


# ==========================================================
#  Build targets
# ==========================================================
.PHONY: all build
all: build

build: $(EXECUTABLE)

$(EXECUTABLE): $(SIM_SRC) $(SRC_HEADERS)
	@mkdir -p $(BUILD_DIR)
	$(QCC) $(CFLAGS) $< -o $@ $(LDFLAGS)

# Health-check binary: identical to production but compiled with DIAGNOSTICS=1.
# Writes pressure_diag.dat per run; not for production use.
.PHONY: build-health
build-health: $(BUILD_DIR)/BioReactor-health

$(BUILD_DIR)/BioReactor-health: $(SIM_SRC) $(SRC_HEADERS)
	@mkdir -p $(BUILD_DIR)
	$(QCC) $(CFLAGS) -DDIAGNOSTICS=1 $< -o $@ $(LDFLAGS)

# Video binary: compiled with VIDEOS=1. Produces mp4s in the run directory.
# Requires a display/OpenGL context; use the SLURM video template (offscreen via Mesa).
.PHONY: build-video
build-video: $(BUILD_DIR)/BioReactor-video

$(BUILD_DIR)/BioReactor-video: $(SIM_SRC) $(SRC_HEADERS)
	@mkdir -p $(BUILD_DIR)
	$(QCC) $(CFLAGS) -DVIDEOS=1 $< -o $@ -L$(BASILISK)/gl -lglutils -lfb_tiny


# MPI binary: compiled with _MPI=1 using mpicc as the CC99 backend.
# Requires OpenMPI to be loaded: module load openmpi
# Submit via slurm_mpi_template.sh (uses srun --mpi=pmix).
# Do NOT use OMP_NUM_THREADS with this binary; parallelism is MPI-only.
.PHONY: build-mpi
build-mpi: $(BUILD_DIR)/BioReactor-mpi

$(BUILD_DIR)/BioReactor-mpi: $(SIM_SRC) $(SRC_HEADERS)
	@mkdir -p $(BUILD_DIR)
	module load openmpi && \
	CC99='mpicc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1' \
	$(QCC) $(CFLAGS) -D_MPI=1 $< -o $@ -L$(BASILISK)/gl -lglutils -lfb_tiny

# Gold-standard production binary: MPI parallelism + inline video generation.
# This is the default binary for all sweeps and SLURM submissions.
.PHONY: build-mpi-video
build-mpi-video: $(BUILD_DIR)/BioReactor-mpi-video

$(BUILD_DIR)/BioReactor-mpi-video: $(SIM_SRC) $(SRC_HEADERS)
	@mkdir -p $(BUILD_DIR)
	module load openmpi && \
	CC99='mpicc -std=c99 -D_XOPEN_SOURCE=700 -D_GNU_SOURCE=1' \
	$(QCC) $(CFLAGS) -D_MPI=1 -DVIDEOS=1 $< -o $@ -L$(BASILISK)/gl -lglutils -lfb_tiny


# ==========================================================
#  Run / submit targets
# ==========================================================
PARAMS ?=

.PHONY: run
run: $(EXECUTABLE)
ifndef PARAMS
	$(error PARAMS is not set. Usage: make run PARAMS=runs/my_run/params.json)
endif
	cd $(dir $(PARAMS)) && $(abspath $(EXECUTABLE)) $(notdir $(PARAMS))

.PHONY: submit
submit:
ifndef PARAMS
	$(error PARAMS is not set. Usage: make submit PARAMS=runs/my_run/params.json)
endif
	@RUN_DIR=$$(dirname $(PARAMS)); \
	mkdir -p $$RUN_DIR logs; \
	if [ "$(DRYRUN)" = "1" ]; then \
	    echo "[DRYRUN] sbatch --job-name=BioReactor --export=PARAMS=$(PARAMS) config/slurm_template.sh"; \
	else \
	    sbatch --job-name=BioReactor \
	           --export=PARAMS=$(PARAMS) \
	           config/slurm_template.sh; \
	fi


# ==========================================================
#  Clean targets
# ==========================================================
.PHONY: clean deepclean
clean:
	rm -rf $(BUILD_DIR)

deepclean: clean
	rm -rf logs runs/*/Data_all runs/*/*.dat
