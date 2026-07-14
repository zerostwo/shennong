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

- Work is active on `feat/analysis-publication-roadmap` against the
  `shennong_analysis_and_publication_figure_roadmap.md` specification.
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
