# GPU port of the bioreactor solver — design

Status: spec only. Work starts after APS DFD 2026 (user decision, 2026-11-14).
Queues: `gpu-debug` and `gpu` only (user decision, 2026-11-14). No other GPU condo.

## 1. Intent

Make the bag solver run on one GPU with no loss of physics we cannot measure. The prize is not a faster L9 in 2D:
512² is too small to fill a modern GPU. The prize is L10+ in 2D and 3D at a cost we can pay, so that the
multi-fidelity ladder can grow upward and into 3D.

The user's question was "is a GPU-friendly extension of Basilisk titanic, or doable by an AI?" Answer from the
survey below: **doable on a uniform grid; titanic on an adaptive tree.** We stay on the uniform grid.

### Goals
- G1. `BioReactor.c` builds with `-grid=multigrid` and runs on one H100 (CUDA driver, FP64).
- G2. Every GPU loop is checked against the same loop on the CPU (CPU Basilisk is the oracle).
- G3. QoIs (kLa_25, tau_95_qss, ediss_mean_qss, dtmix_0.50) match the CPU solver inside a pre-registered
  tolerance (the σ-test protocol: paired runs, 90 % CI inside ±δ, δ = 5 %).
- G4. Embedded boundaries stay. Non-rectangular bags are coming (user, 2026-11-14).
- G5. A measured cost per period (GPU-h and CPU-core-h equivalent) at L9, L10, L11 in 2D, and one 3D level.

### Non-goals
- Adaptive mesh refinement on the GPU (quadtree/octree). Upstream defers it; so do we.
- Multi-GPU. One GPU per run first (FluTAS 2204.08834: best efficiency at the fewest GPUs).
- Changing the CPU solver used for the talk data.

### Anti-goals (failure even if it "works")
- A silent physics change. VOF → CLSVOF is a physics change and needs its own powered test (WP3).
- Silent CPU fallback in the time loop. Every fallback loop is listed; any fallback in the step loop is a bug.
- Patching Basilisk headers in place. New Basilisk goes in a new directory; the installed one stays.

## 2. Quadtree vs multigrid — what the distinction is

Both are grid *data structures* in Basilisk. Neither one means "adaptive" by itself.

| | quadtree (`grid/quadtree.h`) | multigrid (`grid/multigrid.h`) |
|---|---|---|
| What it stores | a tree: each cell can have 4 children, levels can differ in space | one full uniform array per level (L, L−1, …, 0) |
| Adaptive refinement | yes (`adapt_wavelet`, `refine`) | no; every cell is at the finest level |
| GPU support upstream | **none** ("more difficult") | **yes** (2D and 3D, CUDA/HIP/OpenCL/GLSL) |
| What we use now | **yes** — qcc's default 2D grid; `-grid=octree` for 3D | no |

Our runs are uniform: no `adapt_wavelet` call anywhere in `src/`. The talk L9 grid is a uniform 512² *stored as a
quadtree*. So moving to multigrid loses no resolution. What it does lose is every feature that uses the tree:

- `refine (level < L)` in `src/init.h:85` (cross-level warm start) and `src/init.h:127` (replay at a finer level).
  Multigrid cannot refine. Cross-level warm start needs a new path: restore at level L−1, interpolate onto a new
  level-L multigrid (WP5).
- `#if TREE` branches in `embed.h:28,922`, `fractions.h:38,479`, `vof.h:82,343` (refinement/prolongation of
  fractions). On multigrid, `embed.h:752` uses its own `bilinear_embed` prolongation for the Poisson solver.
  So `embed.h` already has a multigrid path on the CPU.
- Dumps: a quadtree dump is not a multigrid dump. Old checkpoints do not restart a GPU run (WP5).
- MPI layout: tree MPI (`tree-mpi.h`) vs multigrid MPI (`multigrid-mpi.h`). Irrelevant for one GPU.

## 3. Survey facts (read 2026-11-14)

Upstream `src/grid/gpu/grid.h` (fetched from basilisk.fr, diffed against our 2026-07-03 install):
- New since our install: 3D multigrid on GPU; CUDA/HIP/OpenCL drivers no longer "experimental"; Metal driver;
  "On Nvidia cards the CUDA driver generally gives better performance".
- CUDA kernels are general C++: "these restrictions [pointers, variable arrays] can be lifted on CUDA, HIP etc.
  drivers". GLSL keeps them.
- Still listed limits: no quadtree/octree; BCs only for 3×3 stencils; "surface tension will not work yet because
  the estimation of curvature relies on code which is not portable to GPUs (parabola.h)"; only simple external
  types; FP64 via `-DDOUBLE_PRECISION`, single precision "can have an important impact on the convergence and
  accuracy of multigrid solvers".
- Function pointers: only field attributes (custom BCs, coarsening/refinement). `embed.h` uses exactly these.
- I/O only on the CPU; frequent CPU/GPU syncs "kill performance". Profiling: `-DTRACE=3` per foreach loop.

Upstream GPU test list (`src/test/README`): multigrid tests include **`rising.c` (rising bubble, CLSVOF)** and
`lid.c`, `axi.c`. `rising.c` uses `two-phase-clsvof.h` + `integral.h` + `tension.h` on `grid/multigrid.h`. So
two-phase flow with surface tension already runs on the GPU upstream, through CLSVOF, not VOF + height functions.
Our install has `rising-clsvof.gpu.ref`.

No upstream GPU test uses `embed.h`. No paper found on Basilisk + GPU + two-phase + embed (arXiv search
2026-11-14). FluTAS (2204.08834, read in full): uniform Cartesian, algebraic VoF (MTHINC), FFT Poisson; 8 A100 ≈ 359
CPU cores; transposes > 50 % of GPU time. Kumar & Vanka 2104.01636 not read yet (WP0).

Our code (`src/`): `#include "embed.h"`, `vof.h`, `curvature.h` (for `position()`), `two-phase.h`, `tension.h`,
`navier-stokes/conserving.h`, `henry_oxy2.h`, `tag.h` (`remove_droplets`). No contact-angle condition on the embed
wall. foreach count / printf count per file: henry_oxy2.h 19/0, replay.h 17/16, diagnostics.h 12/10, init.h 11/9,
forcing.h 6/1, scalars.h 4/1, others ≤ 1.

## 4. Embedded boundaries — are they avoidable?

User claim: embed is unavoidable for non-rectangular bags. Result: **not falsified.** Alternatives checked:
- Volume penalization (Brinkman): trivial on GPU, but its error sits at the wall, O(h) or O(√η) in the wall shear
  stress. τ is our constraint QoI. Rejected for production; kept as a fallback if WP2 fails, with a τ test.
- Curvilinear/mapped grid: no VOF support in Basilisk. Rejected.
- Rectangular multigrid domain only (`dimensions(nx=7, ny=2)`): removes embed and ~71 % solid cells, but only
  for rectangle bags. Not general. Kept as a cheap special case, not as the plan.

## 5. Hardware

FP64 is required (density ratio 1000:1 with multigrid; upstream warns about single precision). Consumer and
graphics cards on Oscar (RTX 3090, A5000, A5500, L40S, A40, A6000) have weak FP64. FP64-capable cards on the
allowed queues: **H100** on `gpu-debug` and `gpu-he`, and what `gpu` offers (check at WP1). B200, A100 and GH200
exist on other queues; not allowed now.

## 6. Work packages (each with a stop rule)

The oracle for every package: the same source built for the CPU (`-grid=multigrid`, no GPU) on a compute node.
"Match" = max abs difference of each field after N steps ≤ 1e-10 relative in FP64, unless stated.

- **WP0 — reading (no compute).** Read Kumar & Vanka 2104.01636 in full; read upstream `cuda.c`,
  `two-phase-clsvof.h`, `integral.h`, `ast/kernels.c` (what turns a foreach into a kernel and when it falls back).
  Output: a list of constructs in our `src/` that each driver rejects.
- **WP1 — upstream build.** Fetch current Basilisk into `/oscar/data/dharri15/eaguerov/basilisk-gpu/` (new dir;
  installed tree untouched). Build `rising.c` CLSVOF with CUDA, FP64, on one gpu-debug H100. Check it against
  `rising-clsvof.gpu.ref`. Time it vs the CPU build on one Genoa node.
  *Stop if* it does not build or does not match the reference.
- **WP2 — embed on the GPU.** Smallest embed case: a CPU-tested upstream embed Stokes/Poiseuille test, built with
  `-grid=multigrid` and CUDA. Run with `-DTRACE=3`, list every loop that falls back to the CPU.
  *Stop and size* if any loop inside the time step falls back: each one is a port item (fix in a project header
  that wraps the upstream function, never an in-place patch). This is the main risk of the project.
- **WP3 — VOF → CLSVOF on the CPU first.** Build our solver on the CPU with `two-phase-clsvof.h` + `integral.h`
  in place of `two-phase.h` + `tension.h`. Check `conserving.h` and `henry_oxy2.h` work with it (`henry_oxy2.h`
  needs a sharp VOF fraction `f`). Pre-registered paired QoI test vs the current solver, σ-test protocol, δ = 5 %.
  *Stop* if FAIL: then the alternative is porting height-function curvature (`parabola.h` → fixed-size arrays, which
  CUDA allows), and that becomes WP3b.
- **WP4 — our headers on the GPU.** Move I/O out of foreach loops (replay.h, diagnostics.h, init.h, restart.h);
  reductions stay (supported). `tag.h` (`remove_droplets`) stays on the CPU at a low frequency or goes off.
  Oracle: CPU multigrid build, field-by-field.
- **WP5 — restarts.** GPU dumps restart GPU runs. Cross-level warm start without `refine()`: restore at L−1 into a
  multigrid, interpolate onto level L. Oracle: the current quadtree cross-level warm start (same QoI test).
- **WP6 — cost.** Cost per period at L9, L10, L11 (2D) and one 3D level, GPU vs the scaling matrix (CPU).
  Output: where the GPU wins, as a figure.

## 7. Effort estimate (a guess; WP1–WP2 replace it with a measurement)

WP0 1 day, WP1 1–2 days, WP2 2 days to a size, WP3 ~1 week incl. queue time, WP4 ~1 week, WP5 3–5 days,
WP6 2 days. Total ~4–5 weeks if WP2 shows no deep embed problem. If embed needs real kernel work, add a size from WP2.

## 8. Open questions

1. Does `two-phase-clsvof.h` work with `navier-stokes/conserving.h` and with embed? (WP3)
2. Does `henry_oxy2.h` need `f` as a VOF tracer, or can it use the CLSVOF fraction? (WP3)
3. Which FP64 GPUs does the `gpu` queue give us, and the queue wait? (WP1)
4. Interface outputs use `position()` from `curvature.h`: CPU-side, at output time only? (WP4)
