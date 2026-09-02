## Figures

Figures are minimalistic.

- **Axes**: label with the physical quantity (and units). Good, keep doing
  this.
- **Title**: the quantity/condition only (e.g. `theta=7deg`), or omit it
  entirely if the axes already say it. Never put status, provenance, or
  caveat text in a title -- "L6/L10 pending rerun", "partial", "as of
  2026-09-02" etc. belong in the diary, the commit message, or a doc
  caption next to the image, not baked into the image itself. The
  filename and its directory already carry versioning context
  (`figure_replicas/replicated_Fig13.png` already says "replica of Fig
  13" -- don't restate "replica" inside the plot too).
- **Legend**: label each series once, tersely. Don't annotate a series
  with a qualifier that is always true of it -- "(published)" after
  "Kim" is redundant, Kim's numbers are always published; "(ramp-matched,
  current driver)" on our own series is redundant, that's the assumed
  baseline for every current run, not a special condition. Only add a
  qualifier when it's the actual point of the comparison (e.g. plotting
  ramp-matched against not-ramp-matched in the same figure).
