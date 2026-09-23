# Shennong Maintainer Status

Last updated: 2026-09-23

This file describes what is true now. Git history and `docs/codex/archive/`
hold point-in-time evidence, `Decisions.md` holds durable rationale, and
`NEWS.md` records user-visible changes.

## Current validation

- API issues #17–#27 have breaking fixes and public-entry regressions for
  result selection, independent direction/significance filtering, consistent
  result envelopes, safe replacement, explicit aggregation, effective seeds,
  shared controls, and discoverable clustering inputs. The preceding workflow
  fixes are tracked in #8–#15; PRs #16 and #28 contain the delivery history.
- CI issue #29 corrects version-provenance assertions and refreshes strict
  backend records to edgeR 4.10.5 and enrichit 0.2.4. The updated remote
  conformance job passes numerical parity and exact-version checks; all four
  R-CMD-check platforms pass on the workflow-fix branch.
- The full local pkgdown build executes the new bundled-data examples and
  renders the complete site. All 9,301 internal links across 326 HTML pages
  resolve; 32 public calls in the introductory recipes pass argument checks.
  Source build and structural R CMD check pass with zero errors and warnings
  and the existing Sankey `head` NOTE. Full-suite and remote CI evidence is
  recorded with the PRs rather than inferred from a structural check.
- Pkgdown deployment now sets an explicit bot author/committer identity; the
  previous main workflow built successfully but failed before publishing (#30).
  Deployment publishes the checked `site/dev` output directly and routes the
  website root to current development documentation.
- No dependencies or exports are added by the documentation refresh. The
  homepage, guide navigation, four introductory/reference articles, related
  workflow examples, and installed usage assets are synchronized. Unrelated
  embedding/heatmap changes remain outside this delivery. Publication is
  verified separately against the merged main commit and pkgdown workflow.

## Current architecture state

- Shared API controls are `batch_by`, `backend_control`, `n_workers`, and
  top-level `seed`; DE/enrichment/Milo return unified results or stored objects.
  Table selection has explicit scope, result IDs resolve only unambiguously,
  and program aggregation names the order of operations. Fate reruns require
  explicit overwrite while retaining metadata ownership checks.
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
