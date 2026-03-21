# Skills Design

## Two Skill Layers

Shennong ships two distinct Codex-facing skill layers.

See also:

- `AGENTS.md`
- `docs/codex/package-governance.md`
- `docs/codex/project-template-design.md`
- `inst/codex/package-skills/`
- `inst/codex/project-template/skills/`

### Package-Usage Skills

Location: `inst/codex/package-skills/`

Purpose:
- teach the agent how to use the Shennong package API
- teach initialization, workflows, stored-result management, and interpretation

### Project Governance Skills

Location: `inst/codex/project-template/skills/`

Purpose:
- teach the agent how to behave inside an initialized analysis project
- standardize runs, scripts, metadata, result export, and memory updates

## Non-Duplication Rule

Package-usage skills must not duplicate project-governance rules.
Project-governance skills must not duplicate package API tutorials.

## Maintainer Rule

When editing a skill, maintainers must decide first whether it belongs to the
package layer or the initialized-project layer.
