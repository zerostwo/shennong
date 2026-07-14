# CodeGraph audit — 2026-06-22 (archived snapshot)

## Scope

This audit continues the requested CodeGraph-based review of the Shennong R package. It uses the repository state in `/home/duansq/dev/packages/shennong` on branch `main`.

## Evidence collected

- CodeGraph index is present and usable under `.codegraph/`.
- CodeGraph database quick counts:
  - files: 63
  - nodes: 1076
  - edges: 3220
  - unresolved_refs: 0
- Current exported API surface from `NAMESPACE`:
  - total exports: 149
  - dot-prefixed exports: 14
- Current source/test size inventory:
  - largest R modules: `R/analysis_clustering.R` (4804 lines), `R/interpretation.R` (4176), `R/analysis_metrics.R` (3938), `R/package_tools.R` (3114), `R/visualization.R` (2796)
  - test files: 18
- Documentation coverage:
  - all exported symbols currently have `.Rd` aliases
  - 24 non-dot exported symbols are missing from `_pkgdown.yml` reference groups
- Validation attempted:
  - `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`
  - result: FAIL 2, WARN 23, SKIP 2, PASS 1399
  - focused repeat: `Rscript -e 'testthat::test_local(filter = "data_io_contracts|signatures", stop_on_failure = FALSE)'`
  - result: the same two failures reproduce

## Current high-level logic

### Data and IO

`sn_load_data()` is the primary example/public-data entry point. It dispatches legacy bundled examples to `.sn_load_example_data()` and public collection samples to `.sn_load_public_data()`. The legacy example path downloads from Zenodo with `sn_download_zenodo()`, reads matrices with `sn_read()`, and initializes Seurat objects with `sn_initialize_seurat_object()`.

`sn_read()` / `sn_write()` wrap `rio` for ordinary formats and dispatch custom formats such as 10x, STARsolo, H5/H5AD, GMT, BPCells, qs, and qs2 through Shennong-specific adapters.

### Preprocessing and clustering

`sn_initialize_seurat_object()` creates Seurat objects with species-aware QC metadata when possible.

`sn_run_cluster()` is the central clustering/integration workflow. It supports RNA and CITE-seq modes, multiple normalization methods, multiple integration backends, stage reuse signatures under `object@misc$sn_run_cluster$stages`, block-gene removal, rare-feature additions, PCA/neighbors/clustering/UMAP, and label-transfer preparation paths.

### Differential expression, enrichment, and stored results

`sn_find_de()` runs marker, contrast, or pseudobulk DE and stores a list under `object@misc$de_results[[store_name]]` containing schema metadata, the result table, analysis/method fields, ranking column metadata, and source parameters.

`sn_enrich()` resolves vector/data-frame/Seurat inputs, runs GO/KEGG/MSigDB enrichment via `clusterProfiler`, and stores results through `sn_store_enrichment()` under `object@misc$enrichment_results[[store_name]]`.

`sn_list_results()` and `sn_get_*_result()` expose a partial retrieval layer for selected `object@misc` collections.

### LLM interpretation

`sn_interpret_annotation()`, `sn_interpret_de()`, and `sn_interpret_enrichment()` prepare evidence, build prompts, call `sn_run_llm()`, parse structured annotation when needed, and store results under `object@misc$interpretation_results[[store_name]]`.

## Findings

### P0 — Current full local test suite is not green

The current full local test command fails:

```text
Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'
[ FAIL 2 | WARN 23 | SKIP 2 | PASS 1399 ]
```

Failures:

1. `tests/testthat/test_data_io_contracts.R:289`
   - Test: `sn_write creates parent directories before dispatching writers`
   - Cause in this environment: `.qs` writer path tries to auto-install `qs` from GitHub (`qsbase/qs`) before reaching the mocked `.export.rio_qs()` writer. GitHub/TLS failure then breaks the unit test.
   - Relevant code: `R/data_io.R:681-774`, especially `.sn_ensure_writer_dependencies()` and `.sn_install_qs_serializer()`.
   - Risk: this test is not hermetic; it can pass only when `qs` is already installed or network/GitHub succeeds.

2. `tests/testthat/test_signatures.R:46`
   - Test: `signature catalog helpers add, update, and delete custom signatures`
   - Cause in this environment: `SignatuR` is not installed, but the test calls `sn_add_signature()` directly. That helper converts the bundled tree back into an upstream SignatuR object through `.sn_signature_tree_to_signatur_db()`, which calls `check_installed_github("SignatuR", "carmonalab/SignatuR")`.
   - Relevant code: `R/signatures.R:105-126`, `R/signatures.R:424-459`.
   - Risk: the test suite requires optional GitHub packages unless the test is skipped or the dependency is mocked.

Recommendation: fix the tests first before deeper refactoring. Either make these tests mock optional dependency installers/checks fully, add explicit skips for optional GitHub dependencies, or refactor the helpers so bundled catalog editing does not require reconstructing the upstream SignatuR API object.

### P1 — Public API surface is very large and still includes internal-looking adapters

`NAMESPACE` exports 149 symbols. Fourteen exported functions are dot-prefixed rio adapters:

```text
.export.rio_bpcells
.export.rio_h5
.export.rio_h5ad
.export.rio_qs
.export.rio_qs2
.import.rio_10x
.import.rio_10x_spatial
.import.rio_bpcells
.import.rio_gmt
.import.rio_h5
.import.rio_h5ad
.import.rio_qs
.import.rio_qs2
.import.rio_starsolo
```

These are documented through `sn_read` / `sn_write` aliases, but they look internal and conflict with the repository guideline that IO adapters should remain internal unless there is a documented reason to keep them exported.

Recommendation: decide whether rio requires these adapter functions to be exported. If yes, document that design explicitly and keep them out of user-facing reference sections. If not, remove `@export` from adapter internals and keep `sn_read()` / `sn_write()` as the public API.

### P1 — `_pkgdown.yml` misses 24 exported non-dot APIs

The following exported symbols are absent from `_pkgdown.yml` reference groups:

```text
sn_call_cell2location
sn_call_cellphonedb
sn_call_infercnvpy
sn_call_mmochi
sn_call_pixi_environment
sn_call_scanvi
sn_call_scarches
sn_call_scpoli
sn_call_scvi
sn_call_spatialdata
sn_call_squidpy
sn_call_stlearn
sn_call_tangram
sn_ensure_pixi
sn_install_pixi
sn_run_cell2location
sn_run_cellphonedb
sn_run_infercnvpy
sn_run_scanvi
sn_run_scpoli
sn_run_spatialdata
sn_run_squidpy
sn_run_stlearn
sn_run_tangram
```

These are mostly Python/pixi backend helpers. They have `.Rd` aliases, but they are not grouped on the reference site.

Recommendation: either add a dedicated `Python and Pixi Backends` reference section or intentionally mark the low-level `sn_call_*` helpers as advanced/internal and narrow the public surface.

### P1 — Large modules concentrate too much behavior

The largest files now carry multiple independent responsibilities:

- `R/analysis_clustering.R` — 4804 lines, 88 functions, 9 exports
- `R/interpretation.R` — 4176 lines, 95 functions, 19 exports
- `R/analysis_metrics.R` — 3938 lines, 73 functions, 20 exports
- `R/package_tools.R` — 3114 lines, 112 functions, 38 exports
- `R/visualization.R` — 2796 lines, 53 functions, 11 exports

Risk: CodeGraph can navigate them, but human review and targeted regression testing are harder because unrelated workflows share files and hidden helpers.

Recommendation: do not do a blind rewrite. Split only around stable seams: clustering stage-cache helpers, Python/pixi backend wrappers, stored-result helpers, and visualization palette/scaffold helpers.

### P1 — Stored-result schema is only partially centralized

There is a useful central storage pair:

- `.sn_store_misc_result()` in `R/interpretation.R:1-7`
- `.sn_get_misc_result()` in `R/interpretation.R:9-16`

It is used for DE, enrichment, interpretation, deconvolution, Milo, communication, regulatory activity, feature annotation, and QC assessment paths.

However, other important workflows write directly into `object@misc`, for example:

- `object@misc$sn_run_cluster$stages`
- `object@misc$integration`
- `object@misc$mmochi`
- `object@misc$bpcells_layers`
- `object@misc$infercnvpy`
- method-named Python backend manifests in `R/package_tools.R`

Risk: `sn_list_results()` exposes only selected collections and can miss integration/backend artifacts. There is no single registry of `object@misc` collection names, schema versions, retrieval helpers, or migration rules.

Recommendation: create a small stored-result registry documenting collection name, schema, writer, reader, and whether it appears in `sn_list_results()`.

### P2 — Optional dependencies and auto-install behavior make tests and workflows non-hermetic

`DESCRIPTION` has a broad `Suggests` surface spanning CRAN, Bioconductor, and GitHub packages. Some workflows also auto-install missing dependencies from inside analysis functions, for example `sn_write()` and Leiden clustering.

Risk: unit tests can trigger network installs or GitHub/Bioconductor resolution unless each optional path is skipped or mocked carefully. The current `.qs` and `SignatuR` failures demonstrate this.

Recommendation: set a rule that package tests must not require live package installation. Tests for installers should mock installers. Runtime workflows can auto-install, but test paths should not.

### P2 — CIBERSORTx local command execution is shell-string based and returns credentials in the command string

`.sn_create_cibersortx_command()` constructs a shell command string with `shQuote()`, and `.sn_run_command()` executes it through `system(command, ...)`.

Relevant code:

- `R/analysis_deconvolution.R:417-423`
- `R/analysis_deconvolution.R:444-513`

Risk is reduced by `shQuote()`, but shell-string execution is still more fragile than `system2(command, args = ...)`. The returned `command` string also includes `--username` and `--token`, which can leak credentials if stored, printed, or embedded in dry-run output.

Recommendation: migrate this path to `system2()` with argument vectors and mask token values in returned metadata/logging.

### P2 — LLM provider execution has limited direct coverage

CodeGraph reports `sn_run_llm()` as called by interpretation functions but without direct test coverage for provider argument filtering and response normalization. Its behavior is central to all LLM interpretation helpers.

Relevant code: `R/interpretation.R:3127-3160`.

Recommendation: add direct unit tests for providers with and without `...`, character responses, structured-only list responses, malformed responses, and tool/structured_type argument filtering.

### P2 — Some small helpers have no targeted tests despite controlling public behavior

Examples from CodeGraph blast-radius output:

- `.sn_enrich_normalize_database_labels()` in `R/analysis_enrichment.R:83-85`
- `.sn_write_zenodo_manifest()` in `R/data_zenodo.R:365-386`
- `.sn_bioconductor_packages()` and `.sn_find_missing_packages()` in `R/package_tools.R`

Recommendation: add narrow tests around these helpers when changing enrichment database semantics, Zenodo manifest schema, or dependency classification.

### P2 — CodeGraph artifacts are currently untracked

Current working tree includes an untracked `.codegraph/` directory. This is useful local state but should usually not be committed into the R package source.

Recommendation: add `.codegraph/` to `.gitignore` and `.Rbuildignore` unless the project intentionally vendors a CodeGraph database.

### P3 — `docs/codex/Status.md` contains outdated snapshot text

The bottom snapshot still mentions older file names and state, such as `integrate.R`, `quality_control.R`, and a single composition test file. The current repository has a broader domain layout and 18 test files.

Recommendation: refresh the snapshot after the next concrete fix so project memory reflects the current architecture.

## Suggested next sequence

Follow-up status: the first three items below were addressed in the same work
session after this audit was written. The remaining items are still useful
medium-term cleanup targets.

1. Fix the two current test failures without broad refactoring:
   - make the `.qs` writer parent-directory test fully mock installer/missing-package checks or skip when `qs` is absent
   - skip/mock the SignatuR-dependent edit test or remove the SignatuR dependency from custom catalog edits
2. Add `.codegraph/` to ignore/buildignore rules if not intended for version control.
3. Add `_pkgdown.yml` coverage for exported Python/pixi APIs or reduce their public exposure.
4. Decide and document whether dot-prefixed rio adapters must be exported.
5. Create a stored-result registry and align direct `object@misc` writers with it.
6. Add direct tests for `sn_run_llm()` provider normalization and selected small helpers with high workflow impact.

## Validation log

Commands run during this audit:

```bash
Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'
Rscript -e 'testthat::test_local(filter = "data_io_contracts|signatures", stop_on_failure = FALSE)'
```

Observed package dependency state:

```text
SignatuR installed: FALSE
qs installed: FALSE
qs2 installed: TRUE
rio installed: TRUE
```
