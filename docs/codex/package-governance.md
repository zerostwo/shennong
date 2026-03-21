# Package Governance

## Scope

The Shennong repository is an R package repository. It is not a user analysis
project.

See also:

- `AGENTS.md`
- `docs/codex/project-template-design.md`
- `docs/codex/skills-design.md`

## Separation Rule

- package-maintainer guidance stays in repository-level `AGENTS.md` and
  `docs/codex/`
- initialized analysis-project governance ships under
  `inst/codex/project-template/`
- package-usage skills ship under `inst/codex/package-skills/`

## Why This Separation Exists

- package development needs package-maintainer governance
- user projects need project-governance assets that do not belong in the package root
- Codex should distinguish package usage from project organization

## Maintainer Rule

Maintainers must not refactor the package root into `runs/`, `results/`,
`memory/`, or other analysis-project directories. Those belong only in the
initialized project template.

Maintainers should evolve initialized-project governance by editing
`inst/codex/project-template/` directly.
