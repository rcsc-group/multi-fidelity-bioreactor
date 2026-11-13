# Solver source guide

The guide to the solver code lives next to the code, so that it changes with it:
[`src/README.md`](https://github.com/rcsc-group/multi-fidelity-bioreactor/blob/main/src/README.md).
It has:

- a map of the files in `src/` (one header per concern);
- the order of the events in one time step, and where each event lives;
- every compile switch, its default, and who uses it;
- how to check that a change does not change results (the golden-run gate below).

Before you change restart code, face fields or events, read [Solver pitfalls](solver-pitfalls.md).

## Golden-run gate

`tools/golden/golden.py` builds the solver from a source directory (default, `DIAGNOSTICS`,
`EXTRA_TRACERS=0`, `RECORD_FLUX`, `REPLAY_FLUX`), runs six short level-6 cases (fresh run, restart,
diagnostics, lean, record, replay of that record) and stores a SHA-256 of every output file. Two
source trees pass when every output is byte-identical.

```bash
# on a compute node, through SLURM (about 5 minutes each)
sbatch tools/golden/golden.sbatch /path/to/base_src  /scratch/golden/base1
sbatch tools/golden/golden.sbatch /path/to/base_src  /scratch/golden/base2   # to find unstable files
sbatch tools/golden/golden.sbatch /path/to/new_src   /scratch/golden/new
sbatch tools/golden/golden.sbatch /path/to/new_src   /scratch/golden/new_c compile   # research switches build

python tools/golden/golden.py compare /scratch/golden/base1/manifest.json /scratch/golden/new/manifest.json \
    --drop-unstable /scratch/golden/base2/manifest.json
```

`logstats.dat` holds wall-clock times, so it differs between two runs of the same source; the
`--drop-unstable` run finds such files and leaves them out of the comparison.

The cases run serially with `OMP_NUM_THREADS=1`. MPI builds and `VIDEOS` runs are not covered.
