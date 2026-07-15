# Shennong Maintainer Status

Last updated: 2026-07-15

## Current structure

- CodeGraph is initialized and synchronized for the repository's R, Python, and
  workflow source. The current index is used for symbol discovery and call-path
  review; Markdown, Rmd, Rd, data, and TOML still require exact text search.
- Package-maintainer documentation is indexed by `docs/codex/README.md`.
  Historical prompts, plans, audits, and long status logs are archived under
  `docs/codex/archive/`.
- Communication and regulatory activity implementations are separated into
  `R/analysis_communication.R` and `R/analysis_regulatory.R`.
- Repository benchmarks live together under `benchmarks/`. Development smoke
  scripts live under `scripts/smoke/` and are not installed with the package.

## Latest cleanup

- The accumulated `Unreleased` roadmap, workflow, compatibility, documentation,
  and publication-figure changes are frozen as the `0.2.0` release candidate.
  `DESCRIPTION` and `NEWS.md` now share that version boundary; after the
  release tag is published, `main` will advance to `0.2.0.9000` so subsequent
  features cannot silently enter the released source. The release gate passes
  `FAIL 0 | WARN 0 | SKIP 6 | PASS 1906`, builds
  `Shennong_0.2.0.tar.gz`, completes structural
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual` with `Status: OK`,
  and validates the complete pkgdown reference index.
- The first remote `0.2.0` pkgdown run reached the evaluated annotation article
  but exposed a clean-runner dependency omission: UCell was not installed for
  the default program-scoring example. The website workflow now declares
  `bioc::UCell`, and the corrected clean-runner pkgdown build and deployment
  passes in GitHub Actions run `29400422339`.
- The corrected run passed the UCell-dependent article and then exposed an
  unconditional Zenodo download in `data-io-projects.Rmd`; Zenodo returned HTTP
  504. That article now requires the explicit
  `SHENNONG_RUN_NETWORK_VIGNETTES=true` opt-in so website builds are not coupled
  to third-party availability.
- The same clean-runner gate also exposed an R 4.6 development-runner failure
  while matching the 47-formal `sn_run_cluster()` dispatcher, plus two
  optional-backend checks that ran before dependency-independent input
  validation. Clustering now uses a compact public compatibility wrapper and a
  single-argument private implementation; its allowlisted tail resolver
  preserves named and positional calls, defaults, and explicit `NULL` values
  without any internal long-formal matching. Scissor/Symphony
  validate required inputs before checking their optional packages. The full
  clustering module passes `FAIL 0 | WARN 0 | SKIP 0 | PASS 442`; the combined
  clustering, abundance, and annotation regression previously passed
  `FAIL 0 | WARN 0 | SKIP 1 | PASS 505`, and the data-I/O article renders with
  network examples disabled.
- Full source-package checking exposed that excluding only CodeGraph contents
  could still leave an empty top-level `.codegraph` directory in the tarball.
  The root directory now has an explicit build-ignore rule, the active index
  socket was stopped before packaging, tarball audit confirms the directory is
  absent, and the rebuilt source completes full `R CMD check --no-manual` with
  `Status: OK`.
- The final pre-merge `testthat::test_local(stop_on_failure = TRUE)` run passes
  with `FAIL 0 | WARN 0 | SKIP 6 | PASS 1906`. Skips are limited to unavailable
  optional lisi, ROGUE, scmap, and zen4R dependencies.
- Pre-merge validation exposed one stale registry expectation and one benign
  `enrichit` qvalue fallback warning on a tiny deterministic ORA fixture. The
  registry test now matches the implemented Slingshot state, and the exact
  backend fallback is muffled without hiding other enrichment warnings;
  focused regression tests pass 116 assertions without warnings.
- The analysis/publication roadmap was completed on
  `feat/analysis-publication-roadmap`, merged into `main` as `1f85f80`, and the
  merged source passed `FAIL 0 | WARN 0 | SKIP 6 | PASS 1906` plus full
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual` with `Status: OK`.
  The merged topic branch was then pruned locally.
- `README.Rmd` now exposes a module-by-module one-command software matrix and
  regenerates `README.md`. The matrix covers 33 public workflow entry points
  and all 63 registry methods, and explicitly distinguishes direct R,
  Shennong-managed pixi/CLI, and external runner/result adapters from current
  local dependency availability.
- The final roadmap adapter pass adds the explicit multimodal entry point,
  keeps backend-specific annotation/trajectory/CellRank functions internal,
  adds direct optional Monocle 3 inference, standardizes external Palantir and
  scCODA/pertpy results, and leaves no method-registry entry marked
  unimplemented. Focused adapter tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 34`, including a real Monocle 3 principal
  graph and pseudotime run. The combined abundance/trajectory/adapter quick
  suite passes 97 assertions, its structural `R CMD check` reports
  `Status: OK`, and the complete pkgdown site rebuild includes the multimodal,
  abundance, and trajectory documentation.
- Milestone D now has six generic figure profiles, automatic specifications for
  500 through 5,000,000 simulated points, publication preflight QA,
  PDF/SVG/TIFF/PNG export, reproducible bundles, result-aware DE/enrichment/GSEA
  figures, core QC/clustering/integration diagnostics, and spec migration for
  the established core and bulk plot families. Figure-engine tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 34`; the quick structural check passes with
  `Status: OK`, and the full pkgdown site includes the publication-figure
  article and all new reference pages.
- Milestone C now includes standalone bulk QC, fixed/mixed-design differential
  expression through edgeR, DESeq2, limma, and dream, pathway scores, WGCNA,
  Cox survival, clinical associations, six result-aware plot families, and a
  matrix/list/`SummarizedExperiment` input contract. The focused suite passes
  with `FAIL 0 | WARN 0 | SKIP 0 | PASS 39`, including real DESeq2, dream, and
  WGCNA backend runs. The quick pre-push source build and structural
  `R CMD check` pass with `Status: OK`, and the complete pkgdown site includes
  the new article and all thirteen bulk reference pages.
- Milestone A now has a shipped method registry, availability diagnostics, and
  a common analysis-result schema with generic store/get/list/delete/validate
  APIs. Existing registered result collections are upgraded compatibly rather
  than copied.
- Milestone B1 now has the unified annotation API, marker/reference consensus,
  SingleR/CellTypist/Seurat/Symphony/scmap/scANVI adapters, confidence and
  runner-up calibration, hierarchical labels, a versioned Cell Ontology
  mapping, low-confidence review, and three result-aware diagnostic plots.
- Milestone B2 now has UCell/AUCell/GSVA/ssGSEA/mean program scoring, signature
  coverage diagnostics, cell metadata storage, sample-aware program tests, and
  result-aware activity/heatmap plots.
- Milestone B3 now has Slingshot topology/pseudotime/lineage probabilities,
  terminal-state and curve storage, optional tradeSeq dynamic/branch tests,
  fitted trends and convergence diagnostics, plus six result-aware plots.
- Milestone B4 now has unified Propeller/permutation/Milo abundance tests,
  sample-level evidence and contributions, Augur-style sample-held-out state
  separability, real Scissor bulk integration, RareQ topology discovery plus
  phenotype association, and two result-aware plots.
- Milestone B5 now has one comparable communication schema across LIANA,
  CellChat, CellPhoneDB, NicheNet, and MultiNicheNet; cross-method consensus and
  concordance; sample-level LR evidence and condition comparisons; retained
  ligand-target evidence; and seven result-aware plot modes.
- Milestone B6 now has multi-restart local NMF discovery, cNMF/Hotspot
  adapters, a real optional GENIE3 backend, explicit pySCENIC/SCENIC/GRNBoost2
  adapters, regulon activity and group-specificity summaries, and six
  result-aware program/GRN plot modes.
- Milestone B7 now has a managed scVelo/CellRank pixi environment, real
  spliced/unspliced velocity inference, projected vectors and transition
  evidence, GPCCA terminal states/fate probabilities/lineage-driver import,
  and two result-aware dynamics plots.
- Milestones C1/C2 now have one spatial dispatcher plus explicit feature,
  domain, neighborhood, deconvolution, mapping, integration, and communication
  APIs. The local path includes Moran's I with permutation evidence,
  memory-bounded KNN graphs, neighborhood enrichment/co-occurrence, and
  distance-constrained communication; nnSVG/BANKSY and existing Python
  backends remain discoverable optional paths.
- The CNV/malignancy/metabolism milestone now wraps inferCNVpy and CopyKAT in
  one stored result, exports chromosome evidence from the pixi backend, derives
  reference-calibrated malignancy/subclone/sample diagnostics, and adds a
  curated sample-aware metabolism workflow plus optional heavy-backend adapters.

- Removed internal helpers that had no callers and removed the now-unused
  `data.tree` and `later` dependencies.
- Consolidated three overlapping Codex design notes into `Governance.md`.
- Hardened the pre-push build path so temporary tarballs/check directories are
  cleaned and source-package contents are audited for local CodeGraph or plot
  artifacts.
- Preserved all pre-existing 2026-07-02 SCTransform `block_genes` source,
  documentation, skill, and test changes.
- Preserved the merged grouped-dotplot coordinate fix and the
  minimal-dependency vector fallback for feature plots without `ggrastr`.
- Aligned the pending data-server integration with the installed ShennongData
  0.2 client contract.
- Declared `msigdbr` in the pkgdown workflow's explicit website dependency
  list so the feature-annotation vignette can execute on a clean CI runner.
- Marked the CITE-seq WNN template as non-executing because the article does
  not construct or ship the required multimodal `pbmc_cite` input object.

## Validation

- Spatial workflow tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 46`, covering spatial autocorrelation,
  adapters, domains, metadata storage, graph/enrichment/co-occurrence evidence,
  communication distance filtering, integration, the dispatcher, registry,
  and eight rendered plot types.
- Velocity/fate R contract tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 20`. A real CPU pixi smoke run using
  scVelo 0.3.4 and CellRank 2.3.2 completed end to end on 80 synthetic cells,
  returning 80 velocity vectors, 2,594 transition edges, and 80 CellRank fate
  probabilities.
- Program-discovery and GRN focused tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 79`, covering local NMF, external program
  adapters, GRN edge standardization, derived and supplied regulon activity,
  group specificity, method discovery, metadata storage, and all plots.
- CNV/metabolism tests pass with `FAIL 0 | WARN 0 | SKIP 0 | PASS 56`, covering
  the actual inferCNVpy import adapter, CopyKAT-style predictions, reference
  calibration, all CNV plots, curated and external metabolism backends,
  sample-level contrasts, registry discovery, and all metabolism plots.
- `scripts/check-prepush.R --filter=cnv-metabolism --quick` passes source build,
  structural `R CMD check` with `Status: OK`, and reference-index validation;
  the complete pkgdown site also rebuilds with the new article and five new
  reference pages.
- Communication tests pass with `FAIL 0 | WARN 0 | SKIP 0 | PASS 48` against
  real NicheNet and MultiNicheNet backends,
  synthetic CellChat/LIANA standardization, CellPhoneDB file parsing,
  sample-level contrasts, consensus/concordance, unified result storage, and
  all communication plots.
- Abundance and state-priority tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 34`. The suite runs real speckle Propeller,
  miloR, RareQ, and Scissor backends plus sample-label permutation and
  sample-held-out separability on deterministic synthetic data.
- The quick pre-push source build and structural check pass with `Status: OK`,
  and the complete pkgdown site rebuilt with the abundance/state-priority
  article and four new reference pages.
- Trajectory tests pass with `FAIL 0 | WARN 0 | SKIP 0 | PASS 29`, including
  deterministic branching and linear Slingshot paths, real tradeSeq GAM fits,
  result-contract retrieval, convergence/trend tables, and rendered ggplot
  grobs. The quick pre-push path passed source build and structural check; the
  initial Rd markup note was corrected before the final validation pass.
- The complete pkgdown site rebuilt successfully with the trajectory article
  and all seven new trajectory/dynamic-gene reference pages.
- Program-scoring tests currently pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 50`, including all five scoring backends,
  sparse mean scoring, coverage/drop behavior, sample-level inference, and
  rendered ggplot grobs.
- `scripts/check-prepush.R --filter=program-scoring --quick` passes the source
  build, structural R CMD check, and reference-index validation with
  `Status: OK`; the complete pkgdown site also rebuilds with the program
  scoring/comparison documentation and four new reference pages.
- Annotation workflow tests currently pass with
  `FAIL 0 | WARN 0 | SKIP 1 | PASS 39`; scmap is the single skip because the
  optional package is not installed locally. SingleR and Symphony CPU paths run
  against deterministic synthetic query/reference objects.
- `scripts/check-prepush.R --filter=annotation-workflow --quick` passes the
  annotation tests, source build, tarball audit, structural R CMD check, and
  pkgdown reference-index validation with `Status: OK`.
- The complete pkgdown site was rebuilt from an installed copy of the current
  source; the expanded annotation article and all eight new reference pages
  rendered successfully.
- The Milestone A result/registry plus interpretation compatibility tests pass
  with `FAIL 0 | WARN 0 | SKIP 0 | PASS 142`.
- pkgdown was rebuilt from an installed copy of the current source; the new
  analysis-method/result-contract article and all existing references rendered
  successfully.
- `scripts/check-prepush.R --filter=analysis-registry-result --quick` passed
  its targeted tests, source build, tarball audit, structural
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check`, and pkgdown reference-index
  validation with `Status: OK`.
- `devtools::document()` completed and updated regulatory Rd source pointers
  after the module split.
- Focused tests for signatures, communication/regulatory activity,
  interpretation, visualization, and shipped Codex skills passed with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 308`.
- The full local suite passed with
  `FAIL 0 | WARN 1 | SKIP 5 | PASS 1429`. The warning is the known tiny-input
  enrichment `qvalue` fallback; skips are for unavailable optional packages.
- The temporary-source pre-push build completed, and its tarball-content audit
  confirmed that `.codegraph`, `Rplots.pdf`, and the development smoke script
  were absent. The structural `R CMD check` completed with `Status: OK` after
  the CI compatibility fixes.
- The complete pkgdown site rebuilt successfully with the local pkgdown template
  override.
- CodeGraph re-synced successfully at 65 indexed files, 1,091 nodes, and 3,015
  edges; the new regulatory module is discoverable and the removed helper names
  have no exact source occurrences.
- PR #4 passed `R-CMD-check`, `test-coverage`, and both Codecov statuses before
  merging into `main`.
- Post-merge pkgdown deployments exposed a missing website-only `msigdbr`
  installation and an unseeded CITE-seq example after a transient jsDelivr SSL
  failure; both reproducibility fixes are being validated on GitHub Actions.

## Deferred local data

Ignored `dev/outputs/` and benchmark input/run caches contain several gigabytes
of untracked data and scripts. They were not deleted automatically; see
`Roadmap.md`.
