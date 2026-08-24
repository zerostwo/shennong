# Structural Roadmap

Last updated: 2026-08-22

## Completed since 2026-07-14

- Backend conformance v1 shipped: machine-readable contracts under
  `inst/conformance/contracts/`, a frozen inventory of the 64 historical
  `implemented: true` methods, six executable C0/C1 pilots (Seurat
  log-normalization, silhouette widths, edgeR quasi-likelihood bulk DE, the
  no-batch Seurat clustering pipeline, MSigDB ORA and GSEA), and clustering/
  enrichment drift repairs found by those pilots.
- Usage observability shipped: SQLite schema v2 with 257 of 267 exports
  instrumented after explicit `sn_enable_usage_tracking()`, consent-gated DBI
  delivery through a local outbox, and an overhead benchmark (~2 ms persisted
  median per tiny call).
- AutoZyme acceleration is now the explicit nine-patch automatic subset behind
  operation/input guards; the fresh-process PBMC30k three-arm benchmark
  recorded measured speedups (UCell 33.26x, LISI 9.11x, NormalizeData 4.49x,
  merge 3.04x, JoinLayers 1.69x, scDblFinder 1.75x) with passing comparators.
- Clustering gained conditional parameter grids, resumable RDS checkpoints,
  manifest schema 2.1.0 with runtime/RSS columns, scib-metrics comparison via
  managed JAX environments, and BPCells-backed scran normalization.
- CellTypist Seurat export preserves sparsity end to end; the repository
  redundancy audit removed 20.71 MiB of reproducible products and replaced the
  stale pkgdown tree with the current real-data site build.

## Immediate hygiene

- Land the current uncommitted working tree on `main`
  as self-contained Conventional Commits before starting new structural work.
- Normalize filesystem ownership of the generated `site/` and `pkgdown/`
  trees so pkgdown rebuilds cannot be blocked by cross-user permissions.

## Next safe refactors

- Split large durable domains only when a focused test boundary exists. The
  current size ranking makes `analysis_clustering.R` (~266 KB) the top
  candidate, followed by `interpretation.R`, `analysis_metrics.R`,
  `package_tools.R`, and `visualization.R`.
- Consolidate the three structurally identical internal helper pairs reported
  by the redundancy audit (`scvi`/`python` run-dir defaults, Seurat/result
  object validators, spatial/dynamics result resolvers).
- Deduplicate repeated `.absolute_path`/`.option`/`.flag` argument parsers in
  standalone scripts while preserving direct execution of each script.
- Add focused source-level tests where CodeGraph cannot infer dynamic
  `object@misc` dispatch or test reachability.
- Add an automated pre-push check for the recurring pkgdown reference-index
  omission failure mode instead of relying on review discipline.

## Correctness milestones

- Grow backend-conformance coverage from six pilots toward the frozen 64-method
  implemented inventory; every newly admitted method still requires its own
  executable contract. Priority ladder per `BackendConformanceAudit.md`:
  grouped GSEA, GO/KEGG resource and ID-mapping parity, real msigdbr fixtures,
  batch/Harmony/CCA/RPCA stage comparators, then fresh-process C2 evidence
  with fixture hashes and real/OOD inputs.
- Close the remaining open audit findings (Milo, spatial registry,
  Python-runner adapters, importers, generic provenance).
- Ship the zero-materialization BPCells writer; BPCells currently crosses into
  an in-memory sparse Matrix at O(nnz).

## Ecosystem milestone

- Execute `NextEcosystemMilestone.md`: deploy the pinned immutable OS/Runtime/
  DB images and complete the first real PBMC3K five-repository end-to-end loop
  (analysis, result bundle upload, promotion, readback, lineage verification).
  Modality claims stay conservative and `deployed` remains false in
  `ecosystem-lock.json` until that loop passes.

## Data and artifact governance

- Govern `dev/outputs/` (2.50 GiB): commit scripts, manifests, and summaries
  as tracked source before archiving any binary research output.
- Make the Coralysis capacity benchmark inputs regenerable by versioning their
  caller-supplied base object, or archive them content-addressed; until then
  1.86 GiB remains retained without exact regeneration guarantees.
- Migrate repeated `data-local/` content (~68.7 MiB exact duplicates) to
  canonical content-addressed artifacts with an executable consumer check.
- Keep benchmark inputs and run logs outside git; retain only scripts, compact
  summaries, and reproducible metadata under `benchmarks/`.

## Pending strategy decisions

- Publishing target: GitHub-only versus CRAN/Bioconductor admission, including
  whether optional backends should split into module packages. Record the
  decision in `Decisions.md`; it determines the future of the ~100-package
  Suggests surface.
- Pinning policy for load-bearing unpinned GitHub Remotes (monocle3, CellChat,
  copykat, and similar) once their conformance contracts exist.

## Guardrails

- Preserve exported API behavior during file splits.
- Do not delete generated `man/` files independently of roxygen sources.
- Do not remove `data/*.rda`, `inst/pixi/`, `inst/codex/`, vignette sources, or
  CI workflows merely because CodeGraph does not index them.
- Re-sync CodeGraph after structural changes and verify both new and removed
  symbols explicitly.
