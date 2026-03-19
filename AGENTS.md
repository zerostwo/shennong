# Repository Guidelines

## Project Structure & Module Organization

- `R/` is organized by analysis domain rather than by strict object type. Current modules include IO (`io.R`), preprocessing (`preprocessing.R`), QC (`quality_control.R`), integration (`integrate.R`), metrics/composition (`metrics.R`, `composition.R`), visualization (`visualization.R`), enrichment/GRN (`enrich.R`, `grn.R`), datasets (`data.R`), and small utilities (`utils.R`, `misc.R`).
- `man/` contains roxygen2-generated `.Rd` files. Treat it as generated output and keep it synchronized with the roxygen comments in `R/`.
- `tests/testthat/` currently has a small unit-test surface; add focused tests near the behavior you change.
- `vignettes/` contains longer workflows. Keep chunks check-safe and avoid unconditional network access or heavyweight setup in examples.
- `data/` stores package datasets. `docs/codex/` stores modernization memory files and is already excluded from package builds via `.Rbuildignore`.
- `inst/codex/skills/shennong/` stores the distributable end-user Codex skill that ships with the installed package. Keep it user-facing and analysis-oriented; do not copy repository-only planning or modernization memory into that tree.
- `_pkgdown.yml` and `.github/workflows/` define the package website and CI entry points.

## Build, Test, and Development Commands

- `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` runs the local test suite.
- `Rscript -e 'testthat::test_local(filter = "composition", stop_on_failure = TRUE)'` runs the current focused tests while iterating on composition-related changes.
- `Rscript -e 'if (requireNamespace("devtools", quietly = TRUE)) devtools::document() else stop("devtools not installed")'` regenerates `NAMESPACE` and `man/` after roxygen changes.
- `R CMD build .` builds the package tarball.
- `R CMD check --no-manual Shennong_*.tar.gz` is the preferred package check after a successful build.

## Coding Style & Naming Conventions

- Preserve public behavior unless a change is explicitly justified, documented, and validated.
- Prefer small, reviewable refactors over repo-wide rewrites.
- Keep exported user-facing functions in the strict `sn_verb_noun` naming family. Do not introduce camelCase, dot.case, or mixed prefixes. Reuse existing `sn_*` naming families whenever possible.
- Keep internal helpers unexported and clearly named.
- Prefer explicit namespace usage or roxygen `@importFrom` entries over hidden dependencies.
- Keep examples deterministic and safe for package checks; wrap network- or dependency-heavy examples in `\dontrun{}` or make them conditional.
- Avoid introducing new package dependencies unless they remove a concrete maintenance problem that cannot be solved more simply.

## Testing Guidelines

- Add or update tests before changing behavior in risky areas.
- Prefer lightweight tests that do not require external downloads or optional heavyweight packages unless the function contract truly depends on them.
- Run the narrowest relevant tests first, then rerun the full local suite before closing a milestone.
- If roxygen, exports, or package metadata change, regenerate documentation and rerun the relevant validation commands.

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
- Current high-risk zones are `R/io.R` (heavy use of `rio` internals and dynamic global assignment), package dependency declarations in `DESCRIPTION`, and sparse test coverage.
- Current modernization order is: inventory API, stabilize tests, normalize metadata and namespace handling, reorganize source files without changing behavior, then tighten internals and documentation.
- Repository-internal developer memory stays in `docs/codex/` and `AGENTS.md`. Installed-package user skills belong under `inst/codex/skills/shennong/`.
