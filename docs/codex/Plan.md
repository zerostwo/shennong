# Shennong Modernization Plan

Last updated: 2026-03-18

## Milestone 0: Audit and memory scaffolding

Status: Completed

Scope:
- Audit repository structure, public API surface, tests, docs, metadata, and risk areas.
- Create and maintain the codex memory files.
- Replace the placeholder `AGENTS.md` with repository-specific guidance.

Files affected:
- `AGENTS.md`
- `docs/codex/Prompt.md`
- `docs/codex/Plan.md`
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`

Acceptance criteria:
- Audit findings are captured in the codex memory files.
- Repository-specific working guidance exists in `AGENTS.md`.
- No package behavior changes are introduced in this milestone.

Validation commands:
- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`

## Milestone 1: Stabilize critical behavior with focused tests

Status: Completed

Scope:
- Expand tests around the most exposed and lightweight behavior first.
- Prioritize `sn_calculate_composition()`, small utilities, and other functions that can be validated without heavyweight external dependencies.

Files affected:
- `tests/testthat/`
- Selected `R/` files only when the tests reveal necessary behavior fixes
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`

Acceptance criteria:
- Critical behavior is covered by tests for success paths and input validation.
- Any bug fix introduced here is behavior-preserving or explicitly documented.
- The local test suite passes.

Validation commands:
- `Rscript -e 'testthat::test_local(filter = "composition", stop_on_failure = TRUE)'`
- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`

## Milestone 1A: Data loading and clustering API consolidation

Status: Completed

Scope:
- Introduce `sn_load_data()` as the generalized example-data loader and re-home current PBMC loading behind it.
- Add a unified multi-method ambient RNA interface to `sn_remove_ambient_contamination()`.
- Remove the separate `sn_quick_cluster()` implementation and merge its single-dataset clustering workflow into `sn_run_cluster()`.
- Add regression tests and smoke validation for `pbmc1k` and `pbmc3k`.

Files affected:
- `R/data.R`
- `R/quality_control.R`
- `R/integrate.R`
- `tests/testthat/`
- `man/`
- `NAMESPACE`
- `DESCRIPTION` if dependency declarations must be corrected for the touched paths
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`

Acceptance criteria:
- `sn_load_data()` can load current PBMC example datasets without regressing existing behavior.
- `sn_remove_ambient_contamination()` supports SoupX and decontX through one documented API.
- `sn_run_cluster()` supports both single-dataset clustering and Harmony-based integration without relying on `sn_quick_cluster()`.
- Local tests pass, and explicit smoke runs with `pbmc1k`/`pbmc3k` complete successfully in this environment.

Validation commands:
- `Rscript -e 'testthat::test_local(filter = "data|cluster|ambient", stop_on_failure = TRUE)'`
- `Rscript -e 'if (requireNamespace("devtools", quietly = TRUE)) devtools::document() else stop("devtools not installed")'`
- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`
- `Rscript inst/scripts/smoke_cluster_pbmc.R`

## Milestone 2: Normalize package metadata and namespace handling

Status: Completed

Scope:
- Reconcile `DESCRIPTION`, roxygen imports, and `NAMESPACE` with actual package usage.
- Address missing or misclassified dependencies without changing runtime behavior.
- Resolve documentation/export mismatches where safe.

Files affected:
- `DESCRIPTION`
- `NAMESPACE`
- `R/*.R` roxygen blocks as needed
- `man/`
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`

Acceptance criteria:
- Declared dependencies match package usage more closely.
- Roxygen-generated outputs are consistent with source annotations.
- Exported functions either have generated documentation or are deliberately reclassified and documented in decisions.

Validation commands:
- `Rscript -e 'if (requireNamespace("devtools", quietly = TRUE)) devtools::document() else stop("devtools not installed")'`
- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`

## Milestone 3: Reorganize source files without behavior changes

Status: Pending

Scope:
- Improve file/module organization while keeping function behavior stable.
- Separate public APIs from internal helpers where it reduces maintenance risk.

Files affected:
- Selected files under `R/`
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`

Acceptance criteria:
- File organization is clearer without changing exported function contracts.
- Tests continue to pass after moves or internal refactors.

Validation commands:
- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`
- `R CMD build .`

## Milestone 4: Harden internal utilities and validation paths

Status: Completed

Scope:
- Improve error handling, argument validation, and hidden side effects in internal helpers.
- Reduce fragile internal coupling in high-risk helpers such as IO and package-install checks.

Files affected:
- `R/utils.R`
- `R/io.R`
- Other selected internal helper files
- `tests/testthat/`
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`

Acceptance criteria:
- Internal helpers have clearer contracts and fewer hidden side effects.
- Behavior changes, if any, are explicitly documented and validated.
- Targeted tests for the touched helpers pass.

Validation commands:
- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`

## Milestone 5: Modernize public functions incrementally

Status: Pending

Scope:
- Update exported functions in small slices, starting with the safest and most testable surfaces.
- Avoid simultaneous broad rewrites across unrelated domains.

Files affected:
- Selected exported functions under `R/`
- Matching tests and docs
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`

Acceptance criteria:
- Each touched public function has preserved or explicitly documented behavior.
- New or updated tests validate each incremental change.

Validation commands:
- `Rscript -e 'testthat::test_local(filter = "<target>", stop_on_failure = TRUE)'`
- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`

## Milestone 6: Clean documentation, examples, and vignettes

Status: Completed

Scope:
- Make examples and vignette chunks check-safe.
- Improve public documentation consistency and pkgdown readiness.

Files affected:
- `R/` roxygen comments
- `man/`
- `README.Rmd`
- `vignettes/`
- `_pkgdown.yml`
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`

Acceptance criteria:
- Documentation matches the implemented API.
- Examples and vignette code can run safely in non-interactive package workflows or are properly guarded.

Validation commands:
- `Rscript -e 'if (requireNamespace("devtools", quietly = TRUE)) devtools::document() else stop("devtools not installed")'`
- `R CMD build .`
- `R CMD check --no-manual Shennong_*.tar.gz`

## Milestone 7: Strengthen CI and package checks

Status: Completed

Scope:
- Add or repair repository-local CI and package-check workflows if they are meant to exist in this repo state.
- Align badges and automation with the repository as checked in.

Files affected:
- `.github/` if introduced or restored
- README and pkgdown config as needed
- `docs/codex/Status.md`
- `docs/codex/Decisions.md`

Acceptance criteria:
- CI configuration matches the repository state and validation expectations.
- README badges do not point to missing or misleading automation.

Validation commands:
- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`
- `R CMD build .`
- `R CMD check --no-manual Shennong_*.tar.gz`
