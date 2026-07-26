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
- `ResultContractAudit.md`: current analytical-result coverage, corrected
  inconsistencies, and intentional runtime-artifact exclusions.
- `Ecosystem.md`: current five-repository development state, interface ledger,
  modality compatibility matrix, shared-contract targets, and cross-repository
  release gates.
- `ecosystem-lock.json`: machine-readable compatibility lock for the five
  implementation revisions, package/source and image digests, shared contract
  schemas, CI evidence, and deployed status. All implementation/image fields
  are populated for the current milestone; `deployed` remains false until live
  revision checks pass.

## Real public-data validation

The tracked validation surface lives under `scripts/real-data/`:

- `sources.yml` records provenance, licenses, deterministic subset rules, and
  the four logical data bundles.
- `prepare-real-data.R` materializes those bundles below the local data root.
- `validate-real-data.R` verifies manifests, checksums, biological dimensions,
  finite matrices, coverage declarations, and source-control boundaries.
- `coverage.csv` maps analysis and visualization functions to a bundle, core or
  extended execution level, backend, and pkgdown article.
- `coverage-exclusions.csv` accounts for every remaining public `sn_*` export
  as IO, runtime, result governance, interpretation, signature management,
  acceleration control, or a low-level backend adapter. Validation fails when
  any new export is neither runtime-covered nor explicitly classified.
- `run-runtime-coverage.R` evaluates mapped articles with tracing and writes a
  local CSV/JSON audit without allowing downloads.
- `benchmark-autozyme.R` measures baseline versus scoped AutoZyme execution and
  verifies numerical equivalence and patch restoration.

The default root is `data-local/pkgdown-real`, overridable with
`SHENNONG_REAL_DATA_DIR`. Everything below that root is ignored by Git and
excluded from package builds; it is validation evidence, never package data.
Use `Rscript scripts/build-pkgdown.R --full --real` for the clean real-output
pkgdown gate after the matrix has been prepared and validated.

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
