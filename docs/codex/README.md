# Codex Maintainer Notes

This directory contains durable repository-maintainer context. It is excluded
from package builds and must not be confused with the Codex assets shipped to
Shennong users under `inst/codex/`.

## Active documents

- `Status.md`: current repository state and latest validation.
- `Decisions.md`: chronological architecture and compatibility decisions.
- `Roadmap.md`: remaining structural work, ordered by risk.
- `Governance.md`: boundaries between package maintenance, initialized-project
  governance, and package-usage skills.

## Historical snapshots

Point-in-time plans, prompts, audits, and superseded status logs live under
`archive/`. They remain useful as historical evidence but are not current
instructions.

## Source of truth

- Repository rules: `AGENTS.md`
- Package source: `R/`, `tests/testthat/`, `vignettes/`, and generated `man/`
- Shipped project template: `inst/codex/project-template/`
- Shipped package skills: `inst/codex/package-skills/`

Use CodeGraph before text search for R, Python, or workflow-code questions.
Use exact text search for Markdown, R Markdown, Rd, data, and other formats that
the current CodeGraph index does not cover.
