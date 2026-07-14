# Codex Asset Governance

## Scope and separation

Shennong is an R package repository, not a user analysis project. Its Codex
assets deliberately occupy three separate layers:

1. Package-maintainer guidance stays in repository-level `AGENTS.md` and
   `docs/codex/`.
2. Initialized analysis-project governance ships under
   `inst/codex/project-template/`.
3. Package-usage skills ship under `inst/codex/package-skills/`.

The package root must not be reorganized into analysis-project directories such
as `runs/`, `results/`, or `memory/`. Those belong only in initialized projects.

## Initialized-project template

`inst/codex/project-template/` is the source of truth used by Shennong project
initialization. It remains minimal and explicitly separates governance, source
data, derived data, runs, exported results, and durable memory. Its required
areas are:

- `AGENTS.md`
- `memory/`
- `docs/standards/`
- `skills/`
- `config/`
- `data/`
- `scripts/`
- `notebooks/`
- `runs/`
- `results/`

When project governance changes, edit the template files directly rather than
reconstructing them in R code.

## Skill layers

Package-usage skills under `inst/codex/package-skills/` teach initialization,
the Shennong API, workflows, stored-result management, and interpretation.
Project-governance skills under `inst/codex/project-template/skills/` teach how
to organize runs, scripts, metadata, exports, and memory inside an initialized
analysis project.

Do not duplicate project-governance rules in package-usage skills, or package
API tutorials in project-governance skills. Decide which layer owns a change
before editing a skill.

## Maintainer checklist

- Keep repository-only context out of installed assets.
- Update shipped templates and skills whenever a relevant user-facing contract
  changes.
- Validate installed asset paths through the package APIs and tests that consume
  them.
