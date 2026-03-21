# Project Template Design

## Purpose

`inst/codex/project-template/` is the shipped source of truth for initialized
analysis projects created from Shennong.

See also:

- `AGENTS.md`
- `docs/codex/package-governance.md`
- `docs/codex/skills-design.md`
- `inst/codex/project-template/AGENTS.md`
- `inst/codex/project-template/docs/standards/BioinformaticsAnalysisConventions.md`

## Design Principles

- the template must be minimal
- the template must be explicit
- the template must separate governance, source data, derived data, runs, and
  exported results
- the template must persist durable knowledge for later sessions

## Required Areas

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

## Maintainer Rule

When the initialized-project governance changes, maintainers must update the
template files directly under `inst/codex/project-template/` rather than
reconstructing them in R code.
