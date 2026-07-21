# Glossary

Plain-language definitions of every technical term used in this documentation.
If you are new to CFD or HPC, read this section first.

**kLa (mass-transfer coefficient)**
The rate at which oxygen dissolves from the air into the liquid in the bag,
measured per unit time. Higher kLa = better mixing and faster oxygenation.
`kLa_25` is the value when the liquid has reached 25 % of full oxygen saturation
and is the standard industrial metric. Units: h⁻¹ (converted from the solver's
internal non-dimensional rate via `3600 / T_bio`, so it is directly comparable
to literature values such as Kim et al. 2024).

**Rocking period (T)**
The time it takes the bag to complete one full back-and-forth rock.
Related to rocking frequency: T = 2π / ω_b. Shorter period = faster rocking.

**Angular frequency (ω_b, omega\_b)**
How fast the bag rocks, measured in radians per second (rad/s).
1 Hz = 2π ≈ 6.28 rad/s. A typical bioreactor runs at 0.3–2 Hz.

**Non-dimensional time**
Simulation time is scaled by the characteristic sloshing timescale T_bio = L / U_bio,
where L is the bag length and U_bio is the typical sloshing velocity.
This makes results independent of the exact bag size or fluid speed, so
one non-dim time unit means roughly "one sloshing timescale has passed."
Physical time in seconds = t × T_bio.

**Fidelity / grid level**
How fine the computational mesh (grid) is. Grid size = 2^fidelity × 2^fidelity cells.
Higher fidelity = more cells = more accurate results, but slower and more memory-intensive.
Fidelity 3 (8×8) takes seconds; fidelity 7 (128×128) takes hours.
See the [Fidelity guide](reference/fidelity-guide.md) for a full table.

**Checkpoint / checkpoint restart**
A snapshot of the full simulation state (fluid velocity, volume fraction, dissolved oxygen,
etc.) saved to a binary file called `checkpoint.dump`.
A later simulation can *restore* this snapshot and continue from that exact point,
skipping the warm-up phase. This is used in chained sweeps to save 70–90 % of compute time.

**Segment**
One SLURM job in a chained sweep. Each segment is a complete, self-contained simulation
that starts either from scratch (segment 0) or from the checkpoint left by the previous segment.

**Chain**
A sequence of segments linked by checkpoint restart. Each segment automatically submits
its successor when it completes successfully (self-submitting chain). No external scheduler
bookkeeping is required.

**n_mix_cycles**
The number of rocking cycles to run *before* injecting oxygen. Used to let the flow
reach a steady state before taking measurements. Typical value: 80 cycles for a fresh start,
10 for a checkpoint restart (flow is already developed).

**t_buffer**
The duration (in non-dim time) of the kLa measurement window after oxygen injection.
Larger t_buffer gives a longer average and more stable kLa estimate, at the cost of more
compute time. Rule of thumb: t_buffer > ln(2) / kLa_nd_expected, where kLa_nd is the
solver's internal *non-dimensional* rate (not the h⁻¹ value reported in `results.json` —
t_buffer is itself a non-dim duration, so it must be sized against the non-dim rate).
At fidelity 7, kLa_nd ≈ 0.03–0.1, so t_buffer = 30–50 is sufficient. At fidelity 5,
kLa_nd ≈ 0.1–0.5, so t_buffer = 10–30 suffices.

**Superellipse (geometry.n)**
The mathematical shape of the bag cross-section. n=2 is a standard ellipse;
n=8 or more looks like a rounded rectangle. The shape is controlled by
`geometry.a` (half-width), `geometry.b` (half-height), and `geometry.n` (roundness).

**Surrogate model**
A fast mathematical approximation of kLa built by fitting to existing simulation data.
Once trained, it can predict kLa for any parameter combination in milliseconds
(vs. hours for a real simulation). Used in Bayesian optimisation to guide where to sample next.

**Bayesian optimisation (BO)**
An algorithm that iteratively chooses which parameter combination to simulate next,
based on a balance between exploring unknown regions (*exploration*) and refining
around the best results found so far (*exploitation*). Uses the surrogate to evaluate
candidates cheaply. The expected improvement (*EI*) acquisition function quantifies
how promising each candidate is.

**DoE (Design of Experiments)**
An initial set of simulation runs that covers the parameter space broadly before
optimisation begins. Used to build the first surrogate model.
Common strategy: Latin Hypercube Sampling (LHS), which spreads points evenly.

**SLURM**
The job scheduler on OSCAR (Brown University HPC). You submit a job with `sbatch`;
SLURM queues it and runs it on a compute node when resources are available.
Jobs communicate via environment variables and output files — never interact
with SLURM interactively from a login node.

**HPC / OSCAR**
High-Performance Computing cluster at Brown University. Always verify you are on a
*compute node* (not a login node) before running expensive simulations.
Login nodes are shared; compute nodes are allocated exclusively per job.
