# docs/ — internal notes (not published)

This folder is for the project team: specs and plans, literature reviews, meeting notes, and the history
of the solver (`history/`). It is not part of the documentation site.

The user documentation is in `docs_site/` (built with `mkdocs build`; see `mkdocs.yml`). Pages in
`docs_site/` must not link into this folder.

- `history/solver-comment-archive.md` — the long comments of the solver as they were before the
  2026-11 clean-up (dates, measurements, investigations), verbatim. The code keeps a short "why"; the
  public, neutral version of each lesson is `docs_site/developer/solver-pitfalls.md`.
- `specs/` — design specs, for example the lean-repo plan.
