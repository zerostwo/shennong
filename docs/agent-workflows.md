# Package development reference

Use the map when locating a workflow and commands when selecting checks.
These are references, not a prerequisite reading list.

## Project Structure & Module Organization

- `R/` is organized by durable workflow domains (83 files as of 2026-09-07).
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
