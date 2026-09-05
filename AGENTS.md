# Repository Guidelines

## Project Structure & Module Organization

- `R/` is organized by durable workflow domains (77 files as of 2026-09-05).
  Analysis domains live in one `analysis_<domain>.R` each: preprocessing,
  clustering (`analysis_clustering.R` plus extracted `analysis_label_transfer.R`,
  `analysis_celltypist.R`, `analysis_integration_backends.R`,
  `analysis_rare_cells.R`, and `analysis_integration_benchmark.R`),
  metrics (`analysis_metrics.R` orchestrators plus extracted
  `analysis_integration_metrics.R`, `analysis_composition.R`,
  `analysis_rogue.R`), annotation
  (`annotation.R`, `annotation_ontology.R`, `feature_annotation.R`),
  differential expression (`analysis_de.R`), enrichment, bulk, priority,
  trajectory and velocity/fate, communication/regulatory,
  spatial, CNV, metabolism/grn/program discovery/scoring, adapters, abundance,
  registry, result contract (`analysis_result.R`, `result_bundle.R`), and
  simulation. Presentation lives in per-domain `plot_<domain>.R` files such as
  `plot_composition.R`, `plot_association.R` (shared numeric geometries), and
  `plot_registry.R` (canonical result/view dispatch), plus `visualization.R` (orchestrators;
  theme/palette/density aspects extracted) and
  `figure_spec.R`/`figure_export.R`. Cross-cutting modules: interpretation
  (`interpretation.R` orchestration plus `interpretation_evidence.R` and
  `interpretation_backend.R`), runtime plumbing (`package_tools.R` plus
  `pixi_runtime.R`, `python_bridge.R`, `codex_project.R`), usage tracking,
  acceleration, IO (`data_io.R`, `data.R`), signatures, MCP,
  `utils.R`, and `zzz.R`. Prefer the domain-specific `analysis_*` / `plot_*`
  file already owning a workflow over adding another generic module; consult
  the current `R/` inventory because the module set intentionally evolves.
- Public naming is strictly `sn_verb_noun`. The verb families carry exact
  contracts: `sn_list_*` enumerates a compact tibble while `sn_get_*`
  materializes a stored object/result; `sn_check_*` diagnoses environment or
  readiness without asserting, `sn_validate_*` asserts a schema contract and
  fails fast, `sn_assess_*` returns a QC analysis, and `sn_calculate_*`
  returns a metric value. `sn_run_*` executes a workflow that stores a unified
  result envelope; raw CLI/runtime adapters use `sn_call_*` and never take a
  Seurat object. `sn_with_*` is reserved for scoped execution helpers.
  Renames must ship behind `.Deprecated()` forwarding shims registered as
  `deprecated_alias` exclusions in the runtime-coverage inventory.
- `man/` contains roxygen2-generated `.Rd` files. Treat it as generated output and keep it synchronized with the roxygen comments in `R/`.
- `tests/testthat/` currently has a small unit-test surface; add focused tests near the behavior you change.
- `vignettes/` contains longer workflows. Keep chunks check-safe and avoid unconditional network access or heavyweight setup in examples.
- `data/` stores package datasets. `docs/codex/` stores package-maintainer Codex docs and modernization memory; it is already excluded from package builds via `.Rbuildignore`. Do not create new files under `docs/codex/` by default: the active-document allowlist is enforced by the architecture-gate test, and historical material belongs in `docs/codex/archive/`. Do not append implementation history to `docs/codex/Status.md`; it describes current state only.
- `inst/architecture/` holds committed architecture-gate baselines (public API, dependencies, source-file sizes) enforced by `tests/testthat/test-architecture-gates.R`. Intentional growth of exports, dependencies, or oversized files must update the baseline in the same change set with rationale in `docs/codex/Decisions.md`.
- `inst/codex/project-template/` stores the shipped initialized-project governance template. `inst/codex/package-skills/` stores the shipped package-usage Codex skills. Keep repository-only planning and modernization memory out of those installed user assets.
- `_pkgdown.yml` and `.github/workflows/` define the package website and CI entry points.

## Build, Test, and Development Commands

- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` runs the local test suite.
- `Rscript -e 'testthat::test_local(filter = "composition", stop_on_failure = TRUE)'` runs the current focused tests while iterating on composition-related changes.
- `Rscript -e 'testthat::test_local(filter = "backend-conformance", stop_on_failure = TRUE)'` runs the static method-admission gate and current micro-differential backend contracts.
- `Rscript -e 'if (requireNamespace("devtools", quietly = TRUE)) devtools::document() else stop("devtools not installed")'` regenerates `NAMESPACE` and `man/` after roxygen changes.
- `R CMD build .` builds the package tarball.
- `R CMD check --no-manual Shennong_*.tar.gz` is the full package check after a successful build. In local environments without optional Suggests, use `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual Shennong_*.tar.gz`.
- `Rscript scripts/check-prepush.R --filter="deconvolution" --quick` runs a fast edit-loop check: targeted tests, source build, and structural `R CMD check` without re-running examples/vignettes/full tests.
- `Rscript scripts/check-prepush.R --filter="deconvolution"` runs the standard local pre-push path. It avoids re-running tests inside `R CMD check` after `test_local()` unless `--check-tests` is supplied.

## Coding Style & Naming Conventions

- Preserve public behavior unless a change is explicitly justified, documented, and validated.
- Prefer small, reviewable refactors over repo-wide rewrites.
- Keep exported user-facing functions in the strict `sn_verb_noun` naming family. Do not introduce camelCase, dot.case, or mixed prefixes. Reuse existing `sn_*` naming families whenever possible.
- Use `sn_list_*` for compact enumeration and `sn_get_*` for materialized
  retrieval; `sn_check_*` diagnoses environment/readiness, `sn_validate_*`
  asserts a contract, `sn_assess_*` performs a QC analysis, and
  `sn_calculate_*` returns a metric value. Reserve `sn_with_*` for scoped
  execution helpers.
- Keep internal helpers unexported and clearly named.
- Prefer explicit namespace usage or roxygen `@importFrom` entries over hidden dependencies.
- Keep examples deterministic and safe for package checks; wrap network- or dependency-heavy examples in `\dontrun{}` or make them conditional.
- Avoid introducing new package dependencies unless they remove a concrete maintenance problem that cannot be solved more simply.

## Testing Guidelines

- Add or update tests before changing behavior in risky areas.
- Any new `implemented: true` analysis method or external backend must include an admitted machine-readable contract under `inst/conformance/contracts/`; do not extend `tests/conformance/legacy-methods.txt`. Follow `docs/codex/BackendConformance.md` and compare the Shennong call with a direct upstream reference using the same inputs, parameters, seed, threads, and dependency version.
- Every public analysis parameter must appear in its backend contract as pass-through, transformed, wrapper-only, unsupported/fallback, or waived with a reason. A new or changed parameter without an executable case or explicit waiver fails admission.
- Prefer lightweight tests that do not require external downloads or optional heavyweight packages unless the function contract truly depends on them.
- Run the narrowest relevant tests first, then rerun the full local suite before closing a milestone.
- If roxygen, exports, or package metadata change, regenerate documentation and rerun the relevant validation commands.
- If you add or change any user-facing function, parameter, stored-result schema, or workflow, you must also update the relevant pkgdown article(s) and the shipped Codex assets under `inst/codex/project-template/` and `inst/codex/package-skills/` in the same change set when they are affected.
- The same user-facing change set must also update `NEWS.md` so the release notes reflect the shipped behavior.
- New stored-result retrieval or interpretation features must be documented at the user level with at least one concrete example showing how to discover the stored result and how to retrieve it from a Seurat object.
- After updating any user-facing feature, rebuild pkgdown locally so the rendered site matches the current package sources. Updating only the vignette/reference source files is not sufficient.

## Commit & Pull Request Guidelines

- Keep each modernization step self-contained: one small change set, one validation pass, and one update to `docs/codex/Status.md` and `docs/codex/Decisions.md`.
- Do not overwrite unrelated working tree changes; this repository may be dirty.
- Record compatibility notes, validation commands, and rationale for non-obvious changes in `docs/codex/`.
- Breaking changes require explicit documentation in `docs/codex/Decisions.md` and `NEWS.md`.
- Use Conventional Commits for all commit messages.
- Preferred commit format: `<type>(<scope>): <summary>`.
- Typical types in this repo: `feat`, `fix`, `refactor`, `docs`, `test`, `build`, `ci`, `chore`.
- Keep the summary imperative and specific, for example: `feat(clustering): consolidate quick and integration workflows`.
- If a change is breaking, add a `!` after the type or scope and explain the break in the commit body and `NEWS.md`.

## Architecture Notes (Quick Map)

- Likely public API surface: IO (`sn_read`, `sn_write`, dataset loading), Seurat initialization/normalization/QC, clustering/integration/annotation, plotting, signatures/enrichment, and composition/metrics.
- Current high-risk zones are data import/export, heavyweight optional-package integrations, and any schema stored under `object@misc`.
- Current modernization order is: inventory API, stabilize tests, normalize metadata and namespace handling, reorganize source files without changing behavior, then tighten internals and documentation.
- Repository-internal developer memory stays in `docs/codex/` and `AGENTS.md`. Installed-package user assets belong under `inst/codex/project-template/` and `inst/codex/package-skills/`.

## Ecosystem Coordination

- Use `docs/codex/Ecosystem.md` as the starting point for work that crosses
  Shennong and any of the four sibling repositories.
- This repository is the coordination center, not a monorepo. The sibling
  roots are `/home/duansq/dev/packages/shennong-data`,
  `/home/duansq/dev/services/shennong.one/shennong-os`,
  `/home/duansq/dev/services/shennong.one/shennong-runtime`, and
  `/home/duansq/dev/services/shennong.one/shennong-db`.
- `/home/duansq/dev/services/shennong.one` is an aggregate architecture
  directory, not a Git repository. Inspect and close each repository
  independently.
- Before cross-repository edits, record each repository's branch, HEAD,
  upstream divergence, dirty files, and published/deployed version. Preserve
  unrelated work in every tree.
- A cross-repository behavior change requires a versioned producer/consumer
  contract, focused tests in each affected repository, and a real integration
  fixture. Update `docs/codex/Ecosystem.md` when interface or modality status
  changes.
- Distinguish source, committed, published, deployed, and end-to-end status.
  Do not infer platform modality support from a file format, registry entry,
  discoverable MCP method, or package-only analysis.
