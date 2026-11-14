// checkpoint.h -- checkpoint.dump at the end of the run, and an early checkpoint before the wall clock ends.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

// checkpoint.phase: the LABEL of this checkpoint (forcing_state.h, label_write): its time, forcing phase, the
// units of the dumped fields and the forcing at that instant. A restart reads it and ramps from that forcing
// to its own params.json; without it the restart stops. Written for every dump (planned and early).
// A side file, not a field in the dump, so the dump format stays the same.
static void write_phase_sidecar (double t_now)
{
  if (pid() != 0)
    return;
  FILE *fp = fopen ("checkpoint.phase", "w");
  if (!fp) {
    fprintf (ferr, "WARNING: cannot write checkpoint.phase; a restart from this dump will be refused\n");
    return;
  }
  label_write (fp, t_now);
  fclose (fp);
}

#if !REPLAY_FLUX
// Early checkpoint: when less than MAXRUNTIME_RESERVE_S of wall clock is left, dump and stop. Checked
// every 10 steps; all MPI ranks take the max of their clocks, so they all dump on the same step (dump()
// is collective). Returning 1 stops run(); main() still closes the output files.
event maxruntime_stop (i += 10) {
  if (_maxruntime_s == HUGE || _emergency_dump)
    return 0;
  double wall = perf.t;
  mpi_all_reduce (wall, MPI_DOUBLE, MPI_MAX);
  if (wall < _maxruntime_s - MAXRUNTIME_RESERVE_S)
    return 0;
  _emergency_dump = 1;
  if (pid() == 0)
    fprintf (ferr, "maxruntime: %.0f s of %.0f s used -- checkpointing at "
             "t=%.6f (phase %.4f of a period) and stopping\n",
             wall, _maxruntime_s, t,
             fmod (_forcing_last.psi/(2.*pi), 1.0));
  write_uv_field (t);
  p.nodump = pf.nodump = false;
  dump (file = "checkpoint.dump");
  p.nodump = pf.nodump = true;
  write_phase_sidecar (t);
  return 1;
}

// Final checkpoint, at t_dump_checkpoint (= t_end, a period boundary: theta = 0).
// p and pf are `nodump` by default; we force them into the dump so a restart starts from the right
// pressure. Returning 1 ends the run here: the forcing and oxygen events have no end time, so without
// this run() would keep stepping after the last output.
event dump_checkpoint (t = t_dump_checkpoint) {
  write_uv_field (t);
  if (pid() == 0)
    fprintf (ferr, "checkpoint: writing checkpoint.dump at t=%.4g\n", t);
  p.nodump = pf.nodump = false;
  write_restart_diagnostic ("pre_dump");
  dump (file = "checkpoint.dump");
  p.nodump = pf.nodump = true;
  write_phase_sidecar (t);
  return 1;
}
#else
// A replay ends here, at the same place in the event order as dump_checkpoint ends a coupled run. The
// condition (i < 2e9) is needed: Basilisk stops the loop when no event with a condition is left, and every
// other replay event has none, so without it the run would end after one step.
event end_replay (i++; i < 2000000000) {
  if (replay_stop || t >= t_end - REC_EPS*(1. + fabs (t)))
    return 1;
  return 0;
}
#endif // !REPLAY_FLUX
