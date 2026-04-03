---
name: use-shennong-project-init
description: Use when initializing a governed bioinformatics analysis project from Shennong package assets, including project templates, governance files, and Codex-readable project structure.
---

# use-shennong-project-init

## Purpose

Teach the agent how to initialize a governed analysis project from the Shennong
package assets.

## When To Use

- starting a new analysis repository with Shennong
- installing or locating project-template governance assets

## Required Inputs

- target project path
- project name
- project objective

## Required Outputs

- initialized analysis project
- populated governance files
- project-ready directory structure

## Rules

- use `sn_initialize_project()` rather than recreating the project structure
  manually when Shennong is available
- distinguish package assets from initialized project files
- do not confuse package skills with project governance skills
- use `sn_get_codex_skill_path()` or `sn_install_codex_skill()` when the task
  is to locate or install shipped package assets rather than initialize a
  project
- once the project is initialized, follow the created `AGENTS.md`,
  `memory/`, and `docs/standards/BioinformaticsAnalysisConventions.md`

## References

- `../_shared/references/package_api_map.md`
- `../_shared/references/workflow_recipes.md`

## Procedure

1. Confirm that `Shennong` is installed.
2. Use `sn_get_codex_skill_path(component = "project_template")` if the
   packaged project template path is needed.
3. Initialize the project with `sn_initialize_project(with_agent = TRUE)`.
4. Use `sn_install_codex_skill()` when the user needs the shipped package
   skills copied into a local agent skill directory.
5. After initialization, follow the created `AGENTS.md` and `memory/`.

## Common Mistakes

- treating the package repository itself as the analysis project
- rebuilding template files by hand
- forgetting to populate `config/default.yaml`

## Examples

```r
library(Shennong)

sn_initialize_project(
  path = "analysis-project",
  project_name = "PBMC pilot",
  objective = "Build a reproducible PBMC analysis workflow.",
  with_agent = TRUE,
  overwrite = TRUE
)
```
