# Archived files

Files removed in the 2026-11 clean-up because no code, test or config used them. They are kept in the
git tag `archive/pre-lean-2026-11`. To get one back:

```bash
git show archive/pre-lean-2026-11:config/chain_mpi_f6.yaml > config/chain_mpi_f6.yaml
```

## config/

| File | What it was |
|---|---|
| `chain_f4_bug2_diag.yaml` | Two-segment chain that reproduced a blow-up at the mixing-time restart. |
| `chain_mpi_f6.yaml` | Fidelity-6 MPI chain, two segments with checkpoint restart. |
| `chain_mpi_fresh_diag.yaml` | Fresh MPI run past the old crash time t = 13.3. |
| `chain_mpi_smoke_diag.yaml` | The smoke-test parameters, run fresh past t = 13.3. |
| `chain_mpi_smoke_f3.yaml`, `chain_mpi_smoke_f3_4ranks.yaml` | MPI chain smoke test at fidelity 3 (16 and 4 ranks). |
| `chain_mpi_smoke_f4.yaml`, `chain_mpi_smoke_f4_n5.yaml` | MPI chain smoke test at fidelity 4 (and with 5 transition cycles). |
| `chain_restart_recovery_baseline_L8.yaml`, `chain_restart_recovery_test_L8.yaml` | Level-8 version of the restart-recovery pair (the level-7 pair is still in `config/`). |
| `sweep_fb_fill_l7.json`, `sweep_fb_fill_l7_v3_f03.json` | Early fill-fraction sweeps at level 7 (superseded by `_v2` and the `_mpi_ckpt` sweeps). |
| `tau_sweep_sentinel.sh` | Watcher script for one old shear-stress sweep. |
