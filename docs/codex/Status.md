# Shennong Maintainer Status

Last updated: 2026-10-27

This file describes what is true now. Git history and `docs/codex/archive/`
hold point-in-time evidence, `Decisions.md` holds durable rationale, and
`NEWS.md` records user-visible changes.

## Package audit validation

- The DE/program audit fixes are included in the consolidated source delivery. A fixed
  snapshot of starting `main` (`7b5602b`) plus the original uncommitted work
  and these fixes passes 5,988 assertions with zero failures/errors, 17
  warnings, and one missing local-public-data fixture skip (776 cases across
  64 test files). Later concurrent spatial/atlas development is outside this
  full-suite claim.
- Empty DE subsets and table selection work without fabricating statistics;
  subset metadata cannot overwrite genes, comparisons, or backend columns.
  Missing/blank subset labels fail before dispatch. Mean program scoring
  preserves multiple programs for one cell/group, and malformed signature
  tables report missing columns explicitly. Two unused extraction helpers
  are removed; compatibility aliases remain.
- Source build and structural `R CMD check` return `Status: OK`, including
  the repaired Sankey `utils::head()` namespace lookup. The check skips
  examples/tests; the complete `testthat::test_local()` pass runs separately.
  Pkgdown was rebuilt in the fixed snapshot and the changed guide text was
  independently read back. The consolidated source delivery includes these
  fixes; validation uses an isolated library without modifying the user's
  installed package. Website publication is verified separately.
- One environment warning requires separate attention: installed glmmTMB was
  built with TMB 1.9.21, while the host has TMB 1.9.23. The audit preserves
  shared libraries. The remaining warnings concern small fixtures, missing
  scale data, or method-specific controls; passing tests do not establish
  equivalence for every optional backend.

## Atlas figure styles

- Three shared atlas themes coexist with 23 function-specific structural styles:
  matrix triangles/rings/annotations/dendrograms/tracks, paired and half/split/
  stacked distributions, supplied effect forests/bars/lollipops, UpSet, and
  bipartite/tree/ribbon/rewiring networks. The existing Seurat heatmap and
  violin APIs dispatch matching structures. Two new exports cover effect
  intervals and set membership. No dependencies were added.
- Three installed catalogs cover six papers, 38 main figures and structural
  layouts with panel-level sources and input requirements. T-cell Science
  image verification and supplementary-figure coverage remain incomplete;
  ATAC, FACS, histology, Venn and circular phylogenetic trees remain external.
- `scripts/smoke/atlas-figure-layouts.R` exports all 23 structural examples to
  `data-local/atlas-layout-review/`; demonstrations are labeled synthetic.
  A separate Nature Fig. 2c render uses all 408 original values with exact CSV
  readback and independently recomputed clustering. Source tables, pairing
  IDs, supplied intervals and topology remain inspectable in figure bundles.
- Focused atlas/heatmap/figure/architecture checks passed 198 assertions (one
  existing boxplot deprecation warning and one absent local-fixture skip).
  Structural/plot-registry/backend-admission checks passed 568 assertions
  without warnings or failures; these runs overlap and are not summed.
  Source build and structural R CMD check report OK; examples/full tests are
  not part of that check. A complete local pkgdown rebuild passed and its
  rendered atlas article/API pages were independently read back under
  `site/atlas-layout-review/dev/`. All 23 PDFs have 8 pt text and no text beyond
  page bounds. Logs and source checksums accompany the review gallery.
- The older 6,138-assertion full-suite run predates the structural extension
  and is not evidence for its new geometry. The structural extension is now
  included in the consolidated source delivery.

## Spatial map validation

- Spatial coordinate, feature, domain, and deconvolution maps use shared
  palette/theme/raster helpers and retain equal coordinate units. Native
  Visium/FOV centroids and stored section boundaries are supported; histology
  and segmentation outlines are outside this coordinate-map implementation.
- Focused spatial/figure/architecture/API-registry validation passes (684
  assertions), followed by a final 60-assertion spatial-style pass with no
  failures, warnings or skips. Both real Visium fixture sections (400 spots
  each) render to PNG/PDF with exact native-coordinate readback. Reproduce with
  `Rscript scripts/smoke/spatial-plots.R`; review artifacts are local under
  `data-local/spatial-plot-review/`.
- Local pkgdown rebuild succeeds from an isolated source snapshot; shared-tree
  concurrent builds had collided on `vignettes/--find-assets.html`. Verified
  rendered pages are under `site/spatial-review/`; nothing was deployed.
- Source build and structural R CMD check pass, with the existing Sankey
  unqualified `head` NOTE. This is focused validation, not a new full-suite
  or external-backend run. The spatial changes are now included in the
  consolidated source delivery.

## Current validation

- The remaining embedding-style branch is merged with current API/documentation
  changes. Focused visualization, real WebGL capture, architecture and registry
  tests pass. The generated parameter inventory and roxygen files are refreshed
  from the combined source; full remote checks run on the merge.
- API issues #17–#27 have breaking fixes and public-entry regressions for
  result selection, independent direction/significance filtering, consistent
  result envelopes, safe replacement, explicit aggregation, effective seeds,
  shared controls, and discoverable clustering inputs. The preceding workflow
  fixes are tracked in #8–#15; PRs #16 and #28 contain the delivery history.
- CI issue #29 corrects version-provenance assertions and refreshes strict
  backend records to edgeR 4.10.5 and enrichit 0.2.4. The updated remote
  conformance job passes numerical parity and exact-version checks; all four
  R-CMD-check platforms pass on the workflow-fix branch.
- The consolidated source passes `R CMD build` and
  `R CMD check --no-manual --no-tests` with `Status: OK`, including examples
  and vignette rebuilding. `_R_CHECK_FORCE_SUGGESTS_=false` permits unavailable
  optional packages; the source test suite runs separately.
- The complete source suite is not claimed green: the marker-DE regression
  records an empty tested-feature background (expected 90 genes). The same
  two assertions fail on unchanged baseline `7b5602b` in a fresh process under
  the current installed dependencies. This predates the consolidated changes;
  retain this limitation rather than silently substituting all input genes.
  The extended full-suite run was stopped after reproducing that baseline
  failure; it is not a completed full-suite pass. The separate DE suite has
  two failed background assertions and three dependent enrichment errors;
  the DE-boundary suite has one failed background assertion and one NULL
  background error. Atlas layouts/styles, spatial plotting, heatmaps, program
  scoring and architecture gates complete without failures in the focused run.
  The existing ggplot2 boxplot deprecation warning remains.
- The complete local pkgdown site is rebuilt from an isolated installation of
  the current source. Atlas resources are available in that installation;
  the atlas/spatial/visualization guides, feature reference and NEWS pages are
  independently read back. Rendered output is synchronized to `site/dev`.
  The user's installed library is unchanged; source and website publication
  remain separate checks.
- Pkgdown deployment now sets an explicit bot author/committer identity; the
  previous main workflow built successfully but failed before publishing (#30).
  Deployment publishes the checked `site/dev` output directly and routes the
  website root to current development documentation. It then explicitly requests
  and verifies the latest Pages build of the published commit (#31, #37),
  tolerating a superseded build only while awaiting the exact expected commit.
- Documentation links and package metadata use the reachable GitHub Pages URL
  instead of the previous host without address records.
- No dependencies or exports are added by the documentation refresh. The
  homepage, guide navigation, four introductory/reference articles, related
  workflow examples, and installed usage assets are synchronized. The broader
  plotting and analysis updates are included in the consolidated source delivery.
  Publication is verified separately against the main commit and pkgdown workflow.

## Embedding validation

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

## Current architecture state

- The embedding-style branch is integrated with the shared API and refreshed
  guides. Nebula/glass support planar default dimensions and explicit real 3D
  coordinates, a locally shipped WebGL viewer, transferable cameras, and
  static capture through Chrome/Chromium. Existing classic plots are retained.
- Shared API controls are `batch_by`, `backend_control`, `n_workers`, and
  top-level `seed`; DE/enrichment/Milo return unified results or stored objects.
  Table selection has explicit scope, result IDs resolve only unambiguously,
  and program aggregation names the order of operations. Fate reruns require
  explicit overwrite while retaining metadata ownership checks.
- `sn_plot_dim()` and `sn_plot_feature()` retain classic rendering and add
  nebula/glass scenes. Default two-dimensional inputs use planar KDE contours
  and a front-facing camera; explicit `dims = 1:3` uses real 3D coordinates. Local WebGL supplies smooth normals, per-pixel
  rim lighting, soft points and bloom. ggsave captures that same renderer via
  chromote/Chrome at draw-time physical size times 600 dpi, with vector labels
  and legends. Local-neighbor bandwidth replaces global-SD KDE smoothing to
  avoid inflated empty shells. Export needs Chrome/Chromium; missing WebGL is
  an explicit error. Interactive mode remains single-panel and returns camera
  changes through downloaded JSON/copied R, not a live R-session binding.
- `R/` contains 88 domain-oriented source files. Large mixed modules have been
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
