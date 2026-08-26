# Comprehensive Package Audit — 2026-08-24

Six parallel audits covering module organization, unified analysis interface,
API naming, stored-result system, backend conformance, and test/documentation
health. This document is written as a work reference: every finding carries a
file:line anchor and a priority, and the action list at the end uses stable
IDs (`AUDIT-xx`) that later change sets should cite.

Scope at audit time: `R/` has 59 files, ~46k lines, 1,122 function
definitions, 267 exports (253 `sn_*` + 14 rio S3 hooks). Source version
0.3.0.9000.

## Overall verdict

The architecture is sound and discipline is visible in the newest generation
of code (`backend_control`, unified result envelope, pixi call layer,
fail-closed conformance gate). The dominant problem is a **generation gap**:
functions written before the conventions converged keep idiosyncratic
signatures, and several mechanisms exist but are decorative rather than
load-bearing. Exactly one finding is bug-class today; everything else is
consistency debt addressable without breaking changes.

---

## 1. Module organization & code health

1. **AGENTS.md documents 12 modules; `R/` actually contains 59 files.**
   Documented drift: usage_tracking.R, acceleration.R, analysis_bulk.R,
   analysis_priority.R, analysis_result.R, result_bundle.R, annotation.R,
   analysis_communication.R, analysis_deconvolution.R, mcp_server.R,
   spatial/cnv/grn/metabolism/velocity families, plot_*.R convention, zzz.R.
2. **Duplicate definition, divergent behavior (bug-class).**
   `.sn_score_programs_aucell` is defined twice:
   - R/analysis_metabolism.R:314 — hardcoded defaults, transposed output,
     filtered signatures;
   - R/program_scoring.R:139 — control-list driven, different `aucMaxRank`
     default, raw `getAUC` orientation.
   Same name + different behavior. Keep the program_scoring.R version as
   canonical only after adding an AUCell differential case.
3. **Misplaced functions** (suggested targets):
   - Python-backed runners in R/package_tools.R:857–1159:
     `sn_run_scarches`/`sn_run_scpoli` → analysis_clustering;
     `sn_run_infercnvpy` → cnv; `sn_run_cellphonedb` → communication;
     `sn_run_cell2location`, `sn_run_tangram`, `sn_run_squidpy`,
     `sn_run_spatialdata`, `sn_run_stlearn` → spatial.
   - `sn_store_enrichment` R/interpretation.R:2683 → analysis_enrichment.R.
   - Retrieval family `sn_list_results` / `sn_get_de_result` /
     `sn_get_enrichment_result` / `sn_get_interpretation_result`
     (R/interpretation.R:411–542) → analysis_result.R.
   - `.sn_get_misc_result` R/interpretation.R:228 and shared validator
     `.sn_validate_seurat_object` R/analysis_communication.R:1 → utils.
   - `sn_simulate` / `sn_simulate_scdesign3` R/analysis_clustering.R:6479,
     6536 → data_examples or a simulation module.
   - `sn_detect_accelerator` R/package_tools.R:1915 → acceleration.R.
4. **Seurat-guard boilerplate ×31**: inline
   `inherits(object, "Seurat")` stops in 14 files while the shared validator
   exists once (R/analysis_communication.R:1) and is used only ~6 times.
5. **Oversized code**: `.sn_run_cluster_impl` is 1,469 lines
   (R/analysis_clustering.R:4071) mixing multi-backend arms
   (coralysis :2200, scanvi/scarches :2216, switch points :1650/:1831/:2413);
   24 functions exceed 150 lines (`sn_enrich` 428, `sn_assess_integration`
   425, `.sn_priority_scissor_impl` 308, ...).
6. **Fragile attach**: `base::attachNamespace("WGCNA")` at
   R/analysis_bulk.R:405 inside `sn_run_wgcna`; prefer explicit
   `getExportedValue` or local wrappers.
7. Duplication hotspots: store/get-result triplets repeated per domain
   (communication.R:797 ≈ deconvolution.R:781 ≈ regulatory.R:203 ≈
   metrics.R:2928) are candidates for one generic parameterized getter.
8. Clean bill: no `library()`/`require()` in R/, no TODO/FIXME/browser()
   residue, no undeclared `pkg::` usage vs DESCRIPTION, rio S3 hook exports
   are intentional, `.onAttach` correctly isolated in R/zzz.R.

## 2. Unified analysis-method interface

Conventions that genuinely hold (preserve): `*_by` grouping suffix family;
`assay`+`layer` pair within a domain; `backend_control=list()` escape hatch;
unified result envelope (`tables$primary`, provenance, `sn_validate_result`);
`return_object=` switch; pixi call layer
(`sn_call_pixi_environment`, R/package_tools.R:645) with hard stop on
non-zero exit (:2420–2422).

Four fault lines:

1. **First argument split**: metrics/plots/IO use `x`; workflow `run_*` use
   `object`. Affected `x`-first: all `sn_calculate_*` (metrics.R), `sn_enrich`
   (enrichment.R:434), `sn_deconvolve_bulk` (deconvolution.R:90),
   `sn_run_milo`, `sn_run_celltypist` (clustering.R:6781), filters
   (preprocessing.R:822,945), all `sn_plot_*`.
2. **`layer` semantic ambiguity**: means `"counts"` in preprocessing.R:1285,
   deconvolution.R:90, clustering.R:3385; means `"data"` in grn.R:165,
   metabolism.R:225, annotation.R:563, program_scoring.R:231,
   communication.R:555; counts-layer arg also spelled `counts_layer`
   (trajectory.R:405) and `spliced_layer/unspliced_layer`
   (velocity_fate.R:246).
3. **Seed/thread/verbosity drift**: seed via `cluster_random_seed=717`
   (clustering.R:3398), `seed=717` (metrics.R:1139),
   `backend_control$seed` fallback NA (trajectory.R:556) vs 717L
   (velocity_fate.R:160, spatial.R:169); threads `ncores` (preprocessing.R:1285)
   / `n_cores` (deconvolution.R:90) / none; verbosity `verbose` /
   `quiet` (clustering.R:6781) / unconditional logging.
4. **Return-shape outliers**: bulk trio returns bare validated lists and does
   not persist (`sn_run_wgcna` R/analysis_bulk.R:398–465,
   `sn_run_survival` bulk.R tail, `sn_assess_bulk_qc`); plot stragglers not on
   the object-or-result pattern: `sn_plot_de`/`sn_plot_enrichment`
   (plot_results.R:49,137), `sn_plot_milo(milo_name=)` (visualization.R:1292),
   `sn_plot_annotation_confidence(store_name=)` (plot_annotation.R:22).
5. **Storage-name argument split**: `store_name=` everywhere except `name=`
   in grn.R:165, program_scoring.R:231, program_discovery.R:188.
6. **Control-bag split**: newer functions use `backend_control`; legacy use
   `transfer_control` (clustering.R:6231), `gibbs_control/opt_control`
   (deconvolution.R:90), `umap_control/wnn_control/cluster_control/
   integration_control` (clustering.R:3385), `method_control`
   (package_tools.R:857+) — same concept, five names.
7. **Method registry is decorative**: all dispatch truth is hardcoded
   `match.arg(method, c(...))` (32 sites); `inst/methods/*.yml` is consumed
   only by `sn_list_methods`/`sn_method_status` and the MCP server
   (mcp_server.R:166–168). No runtime or test asserts parity between registry
   keys and accepted `method` values.
8. Legacy-style stragglers: `sn_run_celltypist` keeps `outdir/prefix/quiet/
   plot_results` with no `return_object`/`result_name`;
   `sn_run_scvi`/`sn_run_scanvi` are thin `sn_run_cluster` delegations whose
   options live in `integration_control`, unlike sibling adapters.

## 3. API naming (`sn_verb_noun`)

239/253 conform (~94%). Violations and friction points:

1. **run/call collision (worst structural overlap)**: 13 `sn_call_*` are raw
   CLI wrappers `(command, args)` sharing one rdname
   (R/package_tools.R:703–779) while identically suffixed `sn_run_*` are
   high-level workflows `(object, ...)`. Candidate: unexport or rename to
   internal `.sn_*`.
2. **Live aliases without lifecycle status**, hardwiring backends against the
   family's own `method=` pattern: `sn_run_spatial_deconvolution` →
   `sn_run_cell2location` (R/analysis_spatial.R:429, dispatched :501);
   `sn_run_spatial_mapping` → `sn_run_tangram` (:434, dispatched :502). Also
   `sn_export_figure` documented alias of `sn_save_figure`
   (figure_export.R:127–130).
3. **Ten noun-first names that roxygen proves are getters/actions**:
   `sn_metabolic_signatures` ("Retrieve..."), `sn_method_status`
   ("Report..."), `sn_annotation_confidence` ("Calibrate..."),
   `sn_annotation_consensus` ("Build..."), `sn_figure_spec` ("Inspect..."),
   `sn_integration_control_template`, `sn_pixi_paths`,
   `sn_pixi_config_path`, `sn_mcp_server(_config)`; plus bare verbs
   `sn_enrich` ("Run gene set enrichment") and `sn_simulate` (no noun).
4. **Bulk subfamily internally inconsistent**: `bulk` is noun-modifier
   everywhere except `sn_deconvolve_bulk`; spatial sibling is
   `sn_run_spatial_deconvolution`.
5. Latent-but-coherent undocumented contracts worth writing down:
   list=enumerate-tibble vs get=materialized-data; check=environment
   diagnostics vs validate=schema assertion vs assess=QC analysis vs
   calculate=metric value; `sn_with_*` = sanctioned withr idiom.
6. Pair asymmetries: `store_X` ⇄ `get_X_result` suffix mismatch across six
   domains; usage store has create/check/flush but no delete;
   `export_result_bundle` deliberately has no import counterpart
   (document it).

## 4. Stored-result system

Architecture: Track A legacy typed collections under registry routing
(`.sn_misc_result_registry()` interpretation.R:1–98; writers
`.sn_store_misc_result` :206; readers `.sn_get_misc_result` :228–240) +
Track B generic `@misc$analysis_results[[type]][[name]]`
(R/analysis_result.R:484–542) under one versioned envelope enforced by
`sn_validate_result` (analysis_result.R:408) with read-time additive upgrade
(`.sn_upgrade_analysis_result`, analysis_result.R:150). Plus Track C:
~21 registered artifact collections (interpretation.R:12–34) outside the
contract. The dual track is deliberate, tested, and symmetric — not drift.

Gaps:

1. `integration_comparison` written unregistered
   (analysis_clustering.R:3961,3977,4043,4059,4063); coralysis stores a full
   SCE under a reduction-name key (analysis_clustering.R:1124).
   `sn_audit_results(include_artifacts=TRUE)` flags them "unregistered".
2. `sn_list_results` (interpretation.R:411–429) enumerates only the 8
   listable collections + generic track; artifacts invisible. Two divergent
   discovery surfaces (`sn_list_results` vs `sn_audit_results`) is itself a
   trap.
3. Half of Track-A writers leave `random_seed = NA` (milo, deconvolution,
   regulatory_activity, qc_assessment, de/enrichment); deconvolution/milo/
   qc lack input/parameters provenance lists entirely.
4. Upgraded Track-A results double-store identical data frames under both
   `table` and `tables$primary` (interpretation.R:163–167,185–187) — memory
   cost for large DE/milo tables.
5. `sn_delete_result` leaves empty generic type containers behind
   (analysis_result.R:842); no artifact deletion path.
6. Generic namespace can be shadowed by artifact-looking types
   (`*_artifact` lands in Track B silently, analysis_result.R:421–429).
7. Bulk results are never persisted (see §2 item 4) — no collision today,
   risk is future-facing.

## 5. Backend conformance system

Flow: hardcoded match.arg → hand-written `inst/methods/*.yml` (64 methods,
12 files, all implemented:true; r=40/pixi=13/adapter=10/cli=1) read by
`.sn_method_registry()` (analysis_registry.R:52) → contracts under
`inst/conformance/contracts/` (**only 6 files**) → gate
tests/testthat/test-backend-conformance-registry.R (C0) + micro-differential
pilots (C1). Baseline files sha256-pinned inside test source (verified
matching). Pilot evidence is real upstream-oracle execution (e.g.
test-backend-conformance-clustering.R builds a genuine sparse Seurat object
and compares a 7-stage Seurat oracle pipeline at 1e-12), not mocks.

Findings:

1. **Coverage hole by construction**: admission equation iterates the YAML
   registry only. Hardcoded surfaces invisible to the gate: 11 integration
   choices at R/analysis_clustering.R:3388 (harmony, unintegrated, coralysis,
   seurat_cca/rpca, scvi, scanvi, scpoli, bbknn, totalvi, mmochi);
   normalization seurat/scran/sctransform (:3387); doublets
   auto/decontx/decontpro/soupx (preprocessing.R:1697); pseudobulk
   DESeq2/edgeR/limma (analysis_de.R:152); plus exported surfaces with no
   registry task (`sn_run_llm`, `sn_run_multimodal`, `sn_run_spatialdata`,
   `sn_find_doublets`, most metrics beyond silhouette).
2. **0/70 subjects admitted; 63/64 registry methods unconformed** (all parked
   in tests/conformance/legacy-pending-methods.txt). Admitted-path machinery
   (fresh-process, fixture hashes, C2 tier) has never executed against a real
   case.
3. Parameter-role taxonomy looser than AGENTS.md claims: gate allows
   `parameter_map` catch-all (clustering contract sweeps ~15 effective params
   into one `"..."` entry); no `unsupported/fallback/waived` roles exist in
   schema despite being documented vocabulary.
4. Five contracts carry `registry_key: null` with domains absent from
   inst/methods/ (permitted by spec, but relationship to inventory is
   implicit).
5. Local softness: without `SHENNONG_CONFORMANCE_STRICT=true`, version pins
   skip and missing backends skip rather than fail
   (helper-backend-conformance.R:67–76) — local green ≠ CI green.
6. No promotion pressure from pilot → admitted.

## 6. Tests & documentation health

Suite: 48 files, 20,084 LOC, 572 `test_that()` blocks. Strengths: stored-result
round-trips with schema checks, fully mocked pixi/python execution (no live
network), heavy deps properly skip-guarded. man/Rd ↔ roxygen sampled 15/15 in
sync; NEWS.md, README.md, pkgdown articles all current as of 2026-08-23.

Gaps:

1. **31 truly untested exports** after excluding rio hooks exercised via
   round-trips. Tier 1: `sn_interpret_de` (interpretation.R:3876),
   `sn_interpret_enrichment` (:3965) — sibling `sn_interpret_annotation` has
   12 refs; `sn_with_usage_tracking` (usage_tracking.R:1916),
   `sn_summarize_usage` (:2165) — headline NEWS features, zero calls in
   tests; `sn_write_figure_legend` (interpretation.R:4171),
   `sn_write_presentation_summary` (:4272); `sn_export_figure`
   (figure_export.R:130).
2. Tier 2 wrapper-only coverage: `sn_run_scvi`/`sn_run_scanvi` public
   signatures never directly exercised; spatial aliases tested only through
   targets; 13-function `sn_call_*` family shares one page, children
   untested.
3. `_pkgdown.yml` missing `sn_call_trajectory` reference entry
   (function defined R/package_tools.R:739).
4. No visual regression testing (no vdiffr anywhere); error-path density
   uneven in spatial/trajectory/deconvolution test files vs metrics/IO.

---

## Prioritized action list (cite these IDs)

### P0 — correctness first
- **AUDIT-01**: Deduplicate `.sn_score_programs_aucell`; make
  program_scoring.R:139 canonical, port needed behavior explicitly, add an
  AUCell differential/conformance case before removing the twin.
- **AUDIT-02**: Give `sn_run_spatial_deconvolution`/`sn_run_spatial_mapping`
  either lifecycle-deprecation badges or real parameterized bodies routed
  through `method=` like every other `run_*` (R/analysis_spatial.R:429–436).
- **AUDIT-03**: Replace `attachNamespace("WGCNA")` (R/analysis_bulk.R:405)
  with explicit `getExportedValue`/local wrappers.

### P1 — non-breaking interface convergence (alias-first)
- **AUDIT-04**: Add `object=` formal alias on every `x`-first analysis
  function; add `store_name=` alias where `name=` is used for storage
  (grn.R:165, program_scoring.R:231, program_discovery.R:188).
- **AUDIT-05**: Introduce top-level `seed=NULL` and `verbose=` on `run_*`
  functions that bury them in control bags; precedence
  `seed > backend_control$seed > task default`; stamp into
  `provenance$random_seed` uniformly (already supported,
  analysis_result.R:307–330).
- **AUDIT-06**: Make the method registry load-bearing: add a two-directional
  test asserting registry `task::name` keys equal each owning function's
  `match.arg` set; then extend the admission gate to harvest hardcoded
  choice sets (integration/normalization/doublet/pseudobulk-DE) so they are
  fail-closed too.
- **AUDIT-07**: Fix plot stragglers: object-or-result pattern for
  `sn_plot_de`, `sn_plot_enrichment`, `sn_plot_milo`, 
  `sn_plot_annotation_confidence`.

### P2 — hygiene and debt
- **AUDIT-08**: Relocate misplaced functions per §1 item 3 (cut-paste moves;
  NAMESPACE stable after `devtools::document()`); route the 31 inline Seurat
  guards through the shared validator.
- **AUDIT-09**: Split analysis_clustering.R: label-transfer block
  (6023–7047) out first; then decompose `.sn_run_cluster_impl` one backend
  arm per PR.
- **AUDIT-10**: Update AGENTS.md module map from 12 files to the real layout;
  document the implicit verb contracts (list/get, check/validate/assess/
  calculate) in AGENTS.md and user docs.
- **AUDIT-11**: Add tests for AUDIT §6 tier-1 untested exports (usage
  tracking wrappers, interpret_de/enrichment, publication writers,
  export_figure); add pkgdown entry for `sn_call_trajectory`.
- **AUDIT-12**: Stored-results hardening: register or redirect stray writers
  (`integration_comparison`, coralysis SCE); capture seed/provenance at weak
  writers; give `sn_list_results` an include-artifacts mode; prune empty
  containers in `sn_delete_result`; reject reserved `*_artifact` types in
  `sn_store_result`.
- **AUDIT-13**: Naming renames behind deprecation shims (major-version
  material): noun-first getters → `sn_get_*`; `sn_enrich` →
  `sn_run_enrichment`; `sn_deconvolve_bulk` → `sn_run_bulk_deconvolution`;
  demote `sn_simulate_scdesign3`; decide run/call family boundary
  (unexport CLI wrappers or rename to internal).
- **AUDIT-14**: Conformance maturation: add `unsupported/fallback/waived`
  roles to contract schema, forbid `parameter_map` catch-alls for named
  formals, promote one subject to admitted to exercise strict machinery,
  consider strict-mode-by-default locally.

## Validation expectations

Any change set working through this list should follow repo governance:
narrowest tests first, then
`Rscript scripts/check-prepush.R --filter="<domain>" --quick`, full suite
before closing a milestone; renames need shims + NEWS.md +
docs/codex/Decisions.md + pkgdown rebuild + inst/codex asset sync per AGENTS.md.
