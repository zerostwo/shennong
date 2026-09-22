# Shennong Maintainer Status

Last updated: 2026-09-22

This file describes what is true now. Git history and `docs/codex/archive/`
hold point-in-time evidence, `Decisions.md` holds durable rationale, and
`NEWS.md` records user-visible changes.

## Current validation

- Two- and three-dimensional embedding styles now use a shared local WebGL shader for
  browser and PDF output. Focused tests cover actual headless captures, exact
  output pixel dimensions, locally shipped dependencies, camera round-trip,
  assay/layer values and disconnected density islands. Visual checks include
  the same 2,000-cell PBMC before/after view and a 6,288-cell Blood2 fixture;
  no cells are duplicated to increase apparent density. This revision passes
  807 focused assertions (including 64 embedding-style assertions; zero failures/warnings/skips), source build, and
  structural package check (one existing Sankey `head` NOTE). The local
  pkgdown site was rebuilt. Browser-downloaded angles were reapplied to a
  PDF whose raster resolution and vector text were independently read back.
  Two-dimensional PBMC dim/feature exports also pass 600 dpi/vector-text
  readback, with all cells and contours constrained to zero depth.
  The complete-suite
  figures below are the earlier baseline, not a full rerun for this revision.
- The complete local test suite passes 5,691 assertions with zero failures,
  11 warnings, and one skip. The warnings are explicit compatibility or
  environment boundaries: one BBKNN/Seurat command-log warning, seven
  Scrublet warnings for ignored scDblFinder-only controls, two scmap warnings
  from a deliberately underspecified small fixture, and one local
  TMB/glmmTMB version warning. The sole skip is an unavailable local
  public-data figure fixture.
- Real backend conformance passes 542 assertions with zero failures, two
  explicit Scrublet parameter-boundary warnings, and no skips. This includes
  direct PopV, Scrublet, clustering, enrichment, pilot, registry, and runtime
  comparisons; it does not imply that every optional method has a live
  upstream oracle.
- Python runtime/resource contracts pass 288 assertions, the public parameter
  matrix passes 720 assertions, architecture gates pass six assertions, all
  83 R source files parse, and all packaged Python scripts compile.
- A clean staged `R CMD build` succeeds, including vignettes, and excludes
  live `.codegraph` sockets and generated `.pixi` environments from the source
  archive. With unavailable Suggests permitted, `R CMD check --no-manual`
  completes with zero errors, zero warnings, and zero notes. Its
  installed-package test phase passes 5,111 assertions with zero failures, 11
  warnings, and 21 expected skips for repository-only CI, conformance,
  real-data, platform-shell, and maintainer resources excluded from source
  packages.
- A complete `pkgdown::build_site()` succeeds, including reference pages,
  affected articles, news, sitemap, redirects, and search index.
- These are local source and installed-package results. Remote CI, publication,
  deployment, and ecosystem end-to-end status must be reported separately by
  exact commit SHA.

## Current architecture state

- `sn_plot_dim()` and `sn_plot_feature()` retain classic rendering and add
  nebula/glass scenes. Default two-dimensional inputs use planar KDE contours
  and a front-facing camera; explicit `dims = 1:3` uses real 3D coordinates. Local WebGL supplies smooth normals, per-pixel
  rim lighting, soft points and bloom. ggsave captures that same renderer via
  chromote/Chrome at draw-time physical size times 600 dpi, with vector labels
  and legends. Local-neighbor bandwidth replaces global-SD KDE smoothing to
  avoid inflated empty shells. Export needs Chrome/Chromium; missing WebGL is
  an explicit error. Interactive mode remains single-panel and returns camera
  changes through downloaded JSON/copied R, not a live R-session binding.
- `R/` contains 83 domain-oriented source files. Large mixed modules have been
  split along stable responsibilities, including bulk design, BBKNN,
  result semantics, gene-symbol preparation, Python artifact validation, and
  infercnvpy integration. Public functions remain in the strict
  `sn_verb_noun` families; compatibility names are deprecated shims.
- All 33 built-in analysis-result types use the schema 2.0.0 envelope and an
  explicit semantic contract. Validation covers required alternatives, types,
  finiteness, ranges, key uniqueness, cross-column invariants, fate-probability
  sums, survival interval ordering, and spatial embedding dimensionality.
  Result identity, discovery, retrieval, plotting, interpretation, upgrade,
  audit, and deletion share the same `analysis_type`/`result_id` model.
- Statistical workflows expose the biological unit and expression scale they
  actually use. Composition, differential expression, bulk interaction
  contrasts, enrichment universes, CNV, spatial, trajectory, communication,
  deconvolution, and program analyses retain the inputs and diagnostics needed
  to understand what was calculated. Bulk interaction estimands are aligned
  across edgeR, DESeq2, limma, and dream.
- Python-backed workflows use 15 committed Pixi lockfiles, exact direct
  dependency pins, bounded artifact imports, explicit cell/feature identity
  checks, owned run directories, and conservative cleanup. CellPhoneDB and
  infercnvpy require declared log-normalized input; count-based workflows
  validate integer-like raw counts. Velocity keeps its owned H5AD by default
  when downstream fate analysis needs it.
- Delimited interchange is record-aware, including quoted fields containing
  embedded newlines. PopV and Scrublet use the same owned-run lifecycle as the
  other Python backends. CIBERSORTx and BayesPrism state and validate their
  expression-scale expectations.
- Package-owned runtime and artifact paths use a cross-platform normalized
  representation. Pixi environment variables, generated Python string
  literals, MCP `Rscript` discovery, publication PDF export, and recursive 10x
  discovery have explicit Windows/macOS/Linux-safe contracts.
- The method registry is load-bearing: implemented methods must have an
  admitted machine-readable conformance contract or remain explicitly pending.
  Architecture gates pin export, dependency, maintainer-document, and
  source-file-size growth.

## Known open boundaries

- BPCells-backed scran normalization is bounded by `max.cluster.size` when
  cluster assignments are supplied and keeps normalized data lazy. Automatic
  scran clustering still materializes the selected layer at O(nnz), and some
  other BPCells-backed operations retain similar compatibility boundaries.
- Live upstream conformance is intentionally strongest for admitted and pilot
  methods. Registry and static parameter gates cover the broader optional
  backend surface, but they are not substitutes for a real oracle run.
- Optional methods still depend on external packages, runtimes, data, licenses,
  credentials, or services. A discoverable method is not necessarily runnable
  in every installation.
- Coralysis capacity-benchmark inputs are not yet regenerable from versioned
  public sources.
- No sibling Shennong service is declared deployed by this package change;
  package validation is not ecosystem end-to-end proof.

## Architecture gates

Repository-level regression gates live in
`tests/testthat/test-architecture-gates.R`, with committed baselines under
`inst/architecture/`. Intentional export, dependency, maintainer-document, or
source-size growth requires a same-change rationale in `Decisions.md`.
