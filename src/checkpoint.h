// checkpoint.h -- Checkpoint: phase sidecar, wall-clock stop, event dump_checkpoint.
// Part of BioReactor.c: included at one fixed place; the place sets the event order.

// [PROJECT ADDED, 2026-09-21] Record the phase this checkpoint was written
// at, next to the checkpoint itself.
//
// The restart path INFERS the writer's phase by snapping to the nearest
// period boundary (t_phase_offset = T_per_st*round(t_ck/T_per_st) - t_ck).
// That is correct only because every checkpoint so far is written AT a
// zero-crossing. An emergency dump lands wherever the clock ran out, and
// snapping it moves the forcing by up to half a period -- the bag resumes at
// a different angle than it stopped at, which is the "escape" failure
// documented above at the t_phase_offset derivation.
//
// So the writer writes its phase down. A sidecar rather than a field inside
// the dump: `dump()` serialises the grid, and adding a scalar to it would
// change the checkpoint format and break every existing checkpoint on disk.
// JSON because params_read.h already parses JSON and the file is read once.
static void write_phase_sidecar (double t_now)
{
  if (pid() != 0)
    return;
  FILE *fp = fopen ("checkpoint.phase", "w");
  if (!fp) {
    fprintf (ferr, "WARNING: cannot write checkpoint.phase; a restart from "
             "this dump will have to infer its phase\n");
    return;
  }
  // t_phase_offset is THIS run's offset; a reader resuming from here needs
  // the total, so it can reproduce the same forcing argument w*(t + offset).
  fprintf (fp, "{\"t\": %.17g, \"t_phase_offset\": %.17g, "
           "\"T_per_st\": %.17g}\n", t_now, t_phase_offset, T_per_st);
  fclose (fp);
}

#if !REPLAY_FLUX
// Stop cleanly when the wall clock is nearly gone.
//
// i += 10 rather than i++ so the reduction is not paid every step; at L10 a
// step is seconds, so ten steps is far finer than the 300 s reserve.
event maxruntime_stop (i += 10) {
  if (_maxruntime_s == HUGE || _emergency_dump)
    return 0;
  double wall = perf.t;
  // Every rank must reach the same verdict on the same step: dump() is
  // collective.
  mpi_all_reduce (wall, MPI_DOUBLE, MPI_MAX);
  if (wall < _maxruntime_s - MAXRUNTIME_RESERVE_S)
    return 0;

  _emergency_dump = 1;
  if (pid() == 0)
    fprintf (ferr, "maxruntime: %.0f s of %.0f s used -- checkpointing at "
             "t=%.6f (phase %.4f of a period) and stopping\n",
             wall, _maxruntime_s, t,
             fmod ((t + t_phase_offset) / T_per_st, 1.0));
  write_uv_field (t);
  p.nodump = pf.nodump = false;
  dump (file = "checkpoint.dump");
  p.nodump = pf.nodump = true;
  write_phase_sidecar (t);
  // Nonzero stops run()'s event loop; main()'s fclose() calls still run, so
  // every .dat file is flushed and the run is postprocessable as far as it
  // got. Never a bare `return;` -- see CLAUDE.md.
  return 1;
}

// Write a Basilisk checkpoint at the first complete period boundary after t_end.
// The checkpoint is always at θ=0 (zero-crossing) — clean phase alignment for
// the next segment's soft-start ramp.  Controlled by t_dump_checkpoint global.
event dump_checkpoint (t = t_dump_checkpoint) {
  write_uv_field (t);
  if (pid() == 0)
    fprintf (ferr, "checkpoint: writing checkpoint.dump at t=%.4g\n", t);
  // p and pf have nodump=true by default in centered.h (they're reconstructed
  // from scratch after a restore, causing a first-step velocity crash).
  // Force them into the checkpoint so the restart Poisson solve starts from
  // the correct pressure and applies only a tiny correction.
  p.nodump = pf.nodump = false;
  write_restart_diagnostic ("pre_dump");
  dump (file = "checkpoint.dump");
  p.nodump = pf.nodump = true;
  // [PROJECT ADDED, 2026-09-21] A planned dump is at a zero-crossing, so the
  // reader's nearest-boundary inference would get it right -- but only if it
  // knows which kind of dump this is, and it cannot. Writing the sidecar
  // unconditionally makes the recorded phase the single path for every
  // checkpoint, so the inference is a fallback for OLD dumps alone rather
  // than a branch that has to be chosen correctly at read time.
  write_phase_sidecar (t);
  // acceleration(i++) and oxygen(t=t_mix;i++) are unconditional (no t<=t_end
  // bound), so Basilisk's events() never marks them done and run() keeps
  // stepping forever past t_dump_checkpoint once every t<=t_end-bounded event
  // (logstats/normcal/movies_output) has expired -- confirmed by direct
  // instrumentation: iter/t kept climbing indefinitely with zero further
  // file output after the checkpoint write, for both fidelity 6 and 7, with
  // both oversubscribed and dedicated MPI ranks.  This is why every SHORT
  // fresh run (short relative to typical L7 t_end) was walltime-killed with
  // no results.json, while long L7 runs happened to fit under the walltime.
  // Returning 1 from an event action is the documented Basilisk idiom for
  // stopping run()'s event loop immediately and cleanly (main()'s fclose()
  // calls after run() still execute) -- see grid/events.h: event_do() treats
  // a nonzero action return as event_stop, and events() then returns 0.
  return 1;
}

#else
// Ends the run, at the same place in the event order as event dump_checkpoint ends the coupled run (before
// the outputs of the step that would start at the end time). The condition (i < INT_MAX) is not decoration: Basilisk's events() stops the loop when no
// event with a condition (or with a time and no increment) is left alive, and every other event of a
// replay is unconditional, so without one the run would end after its first step.
event end_replay (i++; i < 2000000000) {
  if (replay_stop || t >= t_end - REC_EPS*(1. + fabs (t)))
    return 1;
  return 0;
}
#endif // !REPLAY_FLUX
