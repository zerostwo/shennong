# Shennong Maintainer Status

Last updated: 2026-08-30

This file describes what is true now. It is not a change log: Git history and
`docs/codex/archive/` hold point-in-time evidence, `Decisions.md` holds durable
rationale, and `NEWS.md` records user-visible changes.

## Current validation

- Full local testthat suite: `FAIL 0 | WARN 8 | SKIP 1 | PASS 4436` in the
  current local dependency profile (most optional backends installed; skip
  counts therefore differ from profiles with fewer packages). Warnings are
  environmental diagnostics from optional-backend paths, not failures.
- Architecture gates (`test-architecture-gates.R`): passing against the
  committed `inst/architecture/` baselines.
- Backend-conformance tests pass as part of the full local suite.
- Source build and `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual
  --no-tests` complete with no errors, 1 WARNING and 1 NOTE. The warning is
  the pre-existing `sn_simulate_scdesign3` Rd link; the note names `cell` and
  `setNames` in unchanged annotation helpers. Tests run separately above.
  The pre-push script completes successfully; this is not a clean check.
- Last full real-data pkgdown audit (2026-08-26) under `site/dev`: 24 article pages, 65 audited
  figure assets; runtime tracing observes 94/94 declared core functions across
  15/15 mapped articles with zero download attempts. This full runtime audit
  was not rerun for the QC helper extraction.
- Standalone count QC: 32 focused assertions pass (157 across QC,
  preprocessing, BPCells, and architecture tests); the 2,000-cell / 32,738-gene
  local PBMC fixture matches direct column-sum percentages exactly before and
  after a controlled count perturbation, without changing the source object.
  This verifies recalculation, not an ambient-correction algorithm.
- Incremental pkgdown rebuild and rendered QC reference/article checks pass;
  the package is installed locally and a fresh-session exported-helper call
  passes.

## Current architecture state

- `sn_add_qc_metrics()` writes or refreshes three count-based QC percentages
  independently of initialization. Explicit assay/layer selection and optional
  suffixes support corrected-count comparisons without changing expression
  layers, default assay, or `nCount_*`/`nFeature_*` metadata.

- `R/` is organized as the 74-file domain module map documented in
  `AGENTS.md`. The 2026-08-26 decomposition pass split `analysis_metrics.R`
  (into `analysis_integration_metrics.R`, `analysis_composition.R`,
  `analysis_rogue.R`), `interpretation.R` (into
  `interpretation_evidence.R`, `interpretation_backend.R`),
  `package_tools.R` (into `pixi_runtime.R`, `python_bridge.R`,
  `codex_project.R`), and `visualization.R` (theme/palette/density
  extractions) via pure moves with byte-identical function blocks.
  Earlier on 2026-08-25, `analysis_clustering.R` was decomposed into
  `analysis_integration_backends.R` and `analysis_rare_cells.R`, four shared
  pixi/control helpers moved to `utils.R`, the three redundant internal
  helper pairs from the redundancy audit were consolidated, and four
  duplicated function definitions (scArches/scPoli shims, CellPhoneDB,
  infercnvpy) were removed. Public naming follows the strict `sn_verb_noun`
  verb-family contracts; off-family names are deprecated shims. The
  environment-specific `sn_call_*()` pixi aliases are deprecated forwarding
  shims over `sn_call_pixi_environment()`, the single supported runtime
  primitive.
- The method registry (`inst/methods/*.yml`) is load-bearing: a two-directional
  parity test pins registry entries to workflow choice sets with explicit
  waivers. New methods cannot claim `implemented: true` without an admitted
  executable contract under `inst/conformance/contracts/`.
- Analytical results follow the canonical result contract (schema 1.0.0) with
  `sn_audit_results()` / `sn_upgrade_results()` governance; the durable
  coverage inventory lives in `ResultContractAudit.md`.
- Usage tracking is opt-in SQLite-outbox-first infrastructure covering 257 of
  267 exports after explicit enablement; architecture in `UsageTracking.md`.
- ShennongOpt is the acceleration engine behind the same scientific API, with
  an explicit nine-patch automatic subset behind operation/input guards.
- Data distribution belongs to the sibling `shennong-data` repository; real
  public-data validation fixtures live under the ignored `SHENNONG_REAL_DATA_DIR`
  boundary governed by `scripts/real-data/`.

## Known open boundaries

- BPCells layers cross into in-memory sparse `Matrix` at O(nnz); the
  zero-materialization writer is separate work.
- Backend conformance covers six executable pilots against the frozen
  64-method historical inventory; growing that coverage follows the priority
  ladder in `BackendConformanceAudit.md`.
- Coralysis capacity-benchmark inputs (~1.86 GiB retained) are not yet
  regenerable from versioned sources; see `Roadmap.md`.
- No ecosystem service is deployed; `deployed` remains false in
  `ecosystem-lock.json` until the end-to-end loop passes.

## Architecture gates

Repository-level regression gates live in `tests/testthat/test-architecture-gates.R`
with committed baselines under `inst/architecture/`: export growth, invalid
export names, maintainer-document allowlist, source-file size guard, and
dependency growth. Intentional changes update the baseline files in the same
commit set with rationale in `Decisions.md`.
