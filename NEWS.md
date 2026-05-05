All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# Unreleased

### Fixed

- `sn_write()` now creates missing parent directories before dispatching both
  rio and custom writers, so nested `.qs`, `.h5ad`, `.h5`, and BPCells outputs
  no longer fail only because the containing directory does not exist.

### Added

- `sn_calculate_variance_explained()` now ranks metadata variables such as
  platform, study, tissue, and sample by weighted embedding variance explained,
  with single-variable and partial multi-variable modes for batch-effect
  diagnostics.
- `sn_calculate_roe()` now computes observed-over-expected enrichment for
  categorical composition tables from Seurat metadata or data frames, with long
  table output by default and optional matrix output for heatmaps.
- `sn_transfer_labels()` now wraps reference mapping with a query-first API
  for pipe-friendly workflows. The default Seurat anchor workflow is retained,
  and `method = "coralysis"` now projects queries onto Coralysis-trained
  references with `Coralysis::ReferenceMapping()`. `method = "scanvi"` and
  `method = "scarches"` now provide semi-supervised scVI-family label transfer
  through the pixi-managed scverse backend.
- `sn_upload_zenodo()` now uploads reusable data files to Zenodo through
  `zen4R`, with a simple draft-first interface and an automatically uploaded
  Shennong manifest that records dataset version, package version, file sizes,
  md5 checksums, and sha256 checksums for reproducible reuse.
- `sn_download_zenodo()` now downloads reusable files from public Zenodo
  records without requiring a token, with optional token support for
  restricted records. `sn_load_data()` now uses this download layer and accepts
  multiple example datasets such as `dataset = c("pbmc1k", "pbmc3k")`; filtered
  datasets are returned as one merged Seurat object and raw datasets as a named
  list of sparse matrices.
- `sn_simulate()` now provides a method-based simulation entry point.
  `method = "scdesign3"` wraps `scDesign3::scdesign3()` for Seurat or
  SingleCellExperiment inputs and can return simulated counts as a Seurat
  object, SingleCellExperiment, sparse matrix, or raw scDesign3 result.
  `sn_simulate_scdesign3()` remains as the backend-specific wrapper.
- `sn_plot_heatmap()` now draws focused heatmaps for user-selected genes, with
  cell-level and group-averaged modes, optional grouping/splitting, default
  rasterization, hidden cell names/ticks, 8 pt group labels, Paired group-bar
  colors, and automatic scaling of requested features when needed.
- `sn_run_cell_communication()` now wraps real cell-cell communication
  backends: CellChat, NicheNet (`nichenetr`), and LIANA. Results can be stored
  and retrieved with `sn_store_cell_communication()` and
  `sn_get_cell_communication_result()`.
- `sn_run_regulatory_activity()` now runs fast footprint-style activity
  inference with DoRothEA regulons or PROGENy pathway models through
  `decoupleR::run_ulm()`. Results can be stored and retrieved with
  `sn_store_regulatory_activity()` and `sn_get_regulatory_activity_result()`.
- `sn_run_cluster()` now accepts `integration_method` for batch workflows.
  In addition to the historical Harmony path, users can run Coralysis
  multi-level integration or Seurat layer integration with CCA/RPCA through
  the same clustering entry point.
- `sn_run_cluster()` now accepts `integration_method = "scvi"` and
  `"scanvi"`. These backends export selected count data to a pixi-managed
  scverse runtime under `~/.shennong/pixi/`, run the Python model, import the
  latent representation as a Seurat reduction, and then continue Shennong's
  neighbors/clustering/UMAP workflow. The scANVI path requires
  `integration_control = list(label_by = ...)`. Convenience wrappers
  `sn_run_scvi()` and `sn_run_scanvi()` expose the same workflows directly.
- New pixi helpers `sn_check_pixi()`, `sn_install_pixi()`,
  `sn_ensure_pixi()`, `sn_pixi_paths()`, `sn_list_pixi_environments()`,
  `sn_pixi_config_path()`, `sn_prepare_pixi_environment()`,
  `sn_call_pixi_environment()`, `sn_detect_accelerator()`, and
  `sn_configure_pixi_mirror()` expose the Python-runtime setup used by
  scVI/scANVI and future Python backends. Shennong now keeps pixi workspaces
  under `~/.shennong/pixi/`, renders package-bundled configs from
  `inst/pixi/`, can auto-install pixi when missing, selects CPU or CUDA pixi
  environments automatically, and can write China mirror configuration into
  the Shennong `PIXI_HOME`.
- Bundled pixi configs and command helpers are available for concrete Python
  method families including `scvi` (shared by scVI/scANVI), `scarches`
  (shared by scArches/scPoli), `infercnvpy`, `cellphonedb`, `cell2location`,
  `tangram`, `squidpy`, `spatialdata`, and `stlearn`. Use environment calls
  such as `sn_call_cell2location()` or analysis wrappers such as
  `sn_run_tangram()` to run commands inside those managed environments.
- `sn_run_infercnvpy(object = seurat_obj, ...)` now provides an object-level
  infercnvpy workflow: it exports the selected Seurat assay/layer with gene
  positions, runs infercnvpy in the managed pixi environment, and imports CNV
  metadata and optional CNV reductions back into the Seurat object.
- Packaged Python runner scripts now live with their pixi family configs under
  `inst/pixi/<family>/scripts/`. Object-level Seurat workflows are available
  for scArches/scPoli, CellPhoneDB, cell2location, Tangram, Squidpy,
  SpatialData, and stLearn through their corresponding `sn_run_*()` wrappers;
  method-specific Python settings can be supplied with `method_control`.
- `R CMD check` namespace diagnostics are tighter: optional `qs`/`qs2`
  serialization packages are now declared, and previous `ave`, `tail`,
  `target`, and `mor` code-analysis notes have been resolved.
- `sn_run_cluster()` now accepts `hvg_features`, a user-supplied feature list
  that is validated against the object and merged with internally selected
  HVGs and rare-aware features before scaling/PCA. This lets users force rare
  population marker genes into the clustering feature set when global HVG
  selection misses them.
- `sn_sweep_cluster_resolution()` now provides a formal resolution-sweep
  interface for empirically comparing candidate cluster counts across Seurat
  resolutions with metrics such as silhouette width, graph connectivity,
  cluster purity, clustering agreement, and optional ROGUE summaries.
- `sn_list_dependencies()` now reports the package's required and recommended
  R package surface with install status and expected source, and
  `sn_install_dependencies()` can install missing CRAN, Bioconductor, and
  GitHub dependencies in one step.
- `sn_list_10x_paths()` now scans a root directory for 10x Genomics outputs
  and can return `outs/` directories, filtered matrix paths, raw matrix paths,
  H5 files, or `metrics_summary.csv` paths. The default now returns `outs/`
  paths so the result can be passed directly to `sn_initialize_seurat_object()`.
  Returned vectors are now named with inferred sample identifiers.
- `sn_interpret_annotation()` now accepts `label_candidates` so sorted or
  enriched datasets can constrain annotation toward expected cell-type spaces
  such as `ILC1` / `ILC2` / `ILC3` instead of relying only on free-text
  background notes.
- `sn_interpret_annotation()` now supports `annotation_mode = "agentic"` for
  a two-stage workflow: broad lineage/state annotation followed by focused
  refinement on ambiguous or lineage-sensitive clusters. Annotation evidence
  now also exposes a `canonical_marker_snapshot` table so prompts can compare
  lineage-defining markers across clusters without relying only on top-ranked
  DE hits. The `ellmer` path now also uses native structured output for the
  final annotation table and can run a tool-assisted focused-comparison step
  before the refinement pass.
- Annotation heuristics and prompt priors are now more robust for ILC-rich
  datasets: blood ILC workflows now bias `KIT+ ILCP-like` over premature
  mature `ILC3` calls when the type-3 program is incomplete, mixed `T/NK` and
  `NK/ILC3` transitional states are called out more explicitly, and dominant
  hemoglobin programs can now surface as `erythroid contamination`.
- `sn_interpret_annotation()` and `sn_prepare_annotation_evidence()` now
  support `marker_selection = "specific"` and
  `enrichment_selection = "specific"` so annotation evidence can prefer
  cluster-restricted marker genes and pathway terms over generic top-ranked
  features.
- Annotation evidence now adds concise canonical lineage heuristic hints from
  known marker programs so the LLM can use deterministic guardrails such as
  `ILC2-like`, `KIT+ ILC-like`, `T-cell-like`, or `B-cell-like` when those
  programs are clearly supported.
- `sn_assess_qc()` now summarizes overall and per-sample QC status, reports
  current QC risk signals such as failed-QC fractions, doublet rates, and
  decontamination zero-count rates, and can compare a filtered object against a
  pre-filter reference to quantify low-quality-cell removal, doublet removal,
  and clean-cell retention. Reports can be stored under
  `object@misc$qc_assessments`.
- `sn_plot_composition()` now provides a composition-focused bar plot helper
  for grouped proportions, counts, QC pass/fail summaries, and similar
  categorical tables.
- `sn_compare_composition()` now compares sample-level composition between two
  groups and reports mean proportions, differences, log2 fold changes, and
  optional Wilcoxon/FDR statistics per category.
- `sn_run_milo()` now provides a Shennong wrapper around miloR for
  neighborhood-level differential abundance testing between two sample groups
  from a Seurat embedding.
- `sn_list_palettes()` and `sn_get_palette()` now expose the package palette
  registry directly, and `sn_list_palettes()` now includes preview-oriented
  display output plus the `OkabeIto` palette from `ggokabeito`.

### Changed

- Standardized public metadata-selector arguments on the `*_by` family across
  grouping, sample, batch, label, cluster, annotation, condition, reference,
  and cell-type workflows. Deprecated aliases such as `group`, `group_col`,
  `sample_col`, `batch`, `label`, `label_col`, `labels_key`,
  `annotation_col`, `condition_col`, `cluster_col`, `reference_key`,
  `cell_type_key`, `cell_type_col`, `cell_state_col`, `groupby`, and
  `cnv_score_groupby` remain available for compatibility.
- `sn_run_cluster()` now uses a simpler rare-feature interface. The supported
  automatic rare feature methods are `gini` and `local_markers`; less common
  `local_hvg` and `ciara` modes were removed from the clustering wrapper.
  Advanced thresholds are consolidated into `rare_feature_control =
  list(group_max_fraction = ..., group_max_cells = ..., gene_max_fraction =
  ..., min_cells = ...)`. The old scalar threshold arguments remain as
  deprecated compatibility aliases.
- `sn_plot_feature()` now silently replaces Seurat's default expression color
  scale when a Shennong palette is requested, avoiding the noisy duplicate
  colour-scale message.
- `sn_plot_feature()` now exposes additional Seurat 5.5 `FeaturePlot()`
  arguments including `assay`, `dims`, `cells`, `alpha`, `stroke_size`,
  `min_cutoff`, and `raster_dpi`. With `ggrastr` available, `raster = TRUE`
  rasterizes a regular ggplot point layer so `pt_size` behaves like
  `raster = FALSE`, fixing overly large rasterized feature points.
- `sn_plot_dim()` now exposes `label_halo` so users can disable the white label
  halo/background, and `label = TRUE, repel = TRUE` now keeps a repel-aware
  label layer instead of replacing it with fixed-position shadow text.
- `sn_plot_dot()` now uses black colorbar frame/tick styling and suppresses the
  duplicate colour-scale replacement message when applying Shennong palettes.
- `sn_list_dependencies()` and `sn_install_dependencies()` now classify
  `anndataR` and `tidytemplate` as GitHub-hosted optional dependencies and
  `Nebulosa` as a Bioconductor dependency instead of routing them through CRAN.
  Legacy `.qs` support now remains opportunistic when `qs` is already
  installed, while new one-step dependency installation uses `qs2` and avoids
  attempting the archived, R 4.6-incompatible `qs` package.
- `sn_install_dependencies()` now installs required dependencies for
  GitHub-hosted optional packages by default and stops with the package names
  that remain missing after installer warnings, making partial installation
  failures easier to diagnose.
- Reworked the pkgdown article set around a PBMC3k tutorial path, adding
  explicit data/project and visualization articles and rewriting workflow
  articles to explain why each Shennong function is used before showing the
  code. Heavy or credentialed chunks now remain opt-in through
  `SHENNONG_RUN_VIGNETTES=true` so local website builds stay fast.
- `sn_run_cluster(normalization_method = "sctransform", batch_by = ...)` now runs
  SCTransform followed by Harmony integration instead of rejecting
  SCTransform-based integration workflows.

- `sn_run_cluster()` now applies `rare_feature_n` per selected
  `rare_feature_method` before de-duplicating the combined rare-aware feature
  set, matching the documented contract. The stored
  `object@misc$rare_feature_selection` record now keeps both the requested
  `rare_feature_n`, the resolved `rare_feature_control`, and the realized
  `selected_rare_feature_n`.
- `sn_calculate_rogue()` now avoids materializing the full matrix before
  optional subsampling, skips redundant entropy work when grouped ROGUE scores
  are requested, and returns tidy per-cluster or per-sample-per-cluster tables
  when grouping metadata are supplied.
- `sn_list_palettes()` now renders palette names and swatches without overlap
  in plot mode and includes built-in viridis-family palettes in the shared
  palette registry. `sn_plot_dim()`, `sn_plot_feature()`, and the other
  `sn_plot_*()` wrappers now reconcile `aspect_ratio` with `panel_widths` /
  `panel_heights` automatically instead of erroring when fixed panel sizes are
  requested, and Seurat reduction plots now apply axis hiding more reliably
  while `sn_plot_dim(label = TRUE)` adds a white halo behind labels for better
  legibility. `sn_plot_feature()` now also supports
  `mode = "density"` for Nebulosa-style embedding density maps with a
  galaxy-like default theme and shared colorbar collection across multi-feature
  plots. `sn_find_doublets()` now records skipped cells as
  `unresolved` instead of `NA` in the stored class column and orders doublet
  classes as `singlet`, `doublet`, then other levels.
- `sn_initialize_seurat_object()` now accepts the named character vectors
  returned by `sn_list_10x_paths()` and imports all detected 10x samples in one
  call, returning a named list of Seurat objects. Annotation-aware
  `sn_filter_genes()` warnings now report example unmatched feature names, and
  `sn_plot_*()` legends use safer non-negative spacing so legend text does not
  overlap plotted panels.
- `sn_initialize_seurat_object()` now recognizes typical 10x Genomics `outs/`
  directories, reads the filtered matrix automatically, and stores discovered
  source metadata such as `raw_feature_bc_matrix` paths and
  `metrics_summary.csv` contents in `Seurat::Misc(object, "input_source")`.
  `sn_remove_ambient_contamination()` now reuses that stored raw path
  automatically when the selected method can use background droplets and the
  caller leaves `raw = NULL`. When `sample_name = NULL` and the input path is a
  named 10x path, the inferred sample identifier is now written into
  `meta.data$sample` automatically.
- `sn_enrich()` now reuses in-session caches for repeated MSigDB term-table
  loads and repeated SYMBOL-to-ENTREZ conversions, which reduces repeated
  overhead during enrichment-heavy test and analysis sessions.
- `sn_interpret_annotation()` now supports cluster-level functional evidence
  through `enrichment_name`, can incorporate cluster QC summaries into the
  prompt, requests structured annotation JSON by default, stores the parsed
  cluster annotation table, can map normalized cell-type labels plus
  confidence/risk fields back onto Seurat metadata, and now defaults to an
  `ellmer`-backed provider path rather than Shennong-managed local provider
  config. It also now resolves `de_name` automatically when omitted,
  preferring a stored `default` marker result, then a single available DE
  result, and otherwise the most recent marker result. Annotation prompts now
  include the full cluster_by evidence table instead of truncating at eight rows,
  request one record per cluster, use more conservative evidence-grounded
  label selection, can inject candidate-label priors for sorted datasets, and
  can attach cluster-neighborhood geometry from reductions such as UMAP.
  Prompt assembly is now more explicitly markdown-structured, which keeps the
  system/task/evidence sections easier to iterate on as prompt templates.
  `sn_interpret_de()`, `sn_interpret_enrichment()`, `sn_write_results()`,
  `sn_write_figure_legend()`, and `sn_write_presentation_summary()` now share
  the same step-wise progress logging and elapsed-time reporting surface.
- LLM-provider integration is now centered on `ellmer`. The old
  `sn_configure_llm_provider()`, `sn_list_llm_providers()`,
  `sn_get_llm_provider()`, `sn_make_openai_provider()`, and
  `sn_make_sub2api_provider()` compatibility shims have been removed, along
  with the legacy `~/.shennong` provider/history workflow. The supported entry
  point is now `sn_make_ellmer_provider()`, which can explicitly forward
  `reasoning_effort` to compatible GPT-5 chat-completions endpoints. Default
  environment-variable discovery is now limited to `OPENAI_*`; the temporary
  `SUB2API_*` compatibility path has been removed.
- Annotation metadata write-back is now lean by default. High-level
  interpretation writes only the core fields needed for visualization and
  grouping (`label`, `broad_label`, `confidence`, `status`, `risk_flags`)
  into Seurat metadata, while detailed supporting markers/functions/notes stay
  in the stored interpretation result table under `object@misc`.
- `sn_calculate_composition()` now supports multi-column `group_by` values,
  can return proportions, counts, or both through `measure`, preserves factor
  columns from the source metadata, filters returned composition categories by
  `min_cells`, and can sort a single grouping column by a chosen category level
  such as WT proportion.
- `sn_initialize_project()` is now the single project bootstrap entry point. It
  now also writes a repository `.gitignore` plus a generated project `.Rproj`
  file into initialized analysis repositories.
- The built-in `Paired` discrete palette is now overridden by Shennong so its
  brightest yellow swatch uses `#ECD577` instead of Brewer's `#FFFF99`.
- Visualization helpers now share a common discrete-palette resolver. Named
  palettes such as `\"Paired\"` expand automatically when more categories are
  present than the base palette length, and key plotting helpers now expose
  `panel_widths` / `panel_heights` plus consistent axis-label handling.
- `sn_plot_dim()` and `sn_plot_feature()` now default to hiding coordinate
  axes, choose point sizes automatically when `pt_size = NULL`, and keep a
  shared point-size heuristic for small versus large datasets. `sn_plot_dot()`
  now defaults to a warmer-high / cooler-low color direction, keeps the
  Z-score legend ahead of the percent legend, supports `Min`/`Max` legend
  labels, uses hollow black-edged percent legend dots, accepts
  `legend_position`, and uses thicker black colorbar ticks.
- Continuous-color helpers now use the same palette registry through
  `sn_get_palette(..., palette_type = "continuous")`, and expression-oriented
  plotting functions now apply continuous palettes through the shared internal
  resolver instead of separate ad hoc `scale_*_distiller()` logic.
- `sn_plot_barplot()` now supports automatic summary bars for repeated
  observations, optional SD/SE error bars, and optional jittered raw points,
  making it suitable for sample-level effect summaries as well as simple
  identity bars.
- `sn_compare_composition()` now adds a `change` factor with levels
  `Increase` and `Decrease` derived from the sign of `log2_fc`.
- `sn_find_doublets()` now skips zero-count and low-feature cells before
  running `scDblFinder()` on corrected layers, records corrected-layer results
  with `_corrected` suffixes, and works with the zero-count flags produced by
  ambient-RNA correction.
- `sn_filter_cells()` now validates its `method` argument, keeps constant-value
  QC groups when MAD collapses to zero, and checks plotting dependencies
  explicitly before rendering diagnostics. `sn_filter_genes()` now validates
  `min_cells` and keeps its threshold summary stable when the requested
  threshold exceeds the number of cells in the object.
- `sn_run_cluster()` now defaults `hvg_group_by` to the `batch` column when
  batch integration is requested and the user does not explicitly supply a
  separate HVG grouping variable.
- `sn_standardize_gene_symbols()` now also accepts character vectors of gene
  symbols or gene IDs and returns the standardized vector directly, while
  preserving the existing matrix and Seurat-object behavior.
- `sn_run_cluster(normalization_method = "scran", batch_by = ...)` now runs
  scran normalization before the selected batch-integration backend instead of
  rejecting batch workflows.
- The repository now ships `scripts/check-prepush.R` so maintainers can run
  documentation, targeted tests, the full test suite, `R CMD build`, and
  `R CMD check --no-manual` in one local pre-push command.
- Internal DE result storage now reuses the shared misc-result helper instead
  of maintaining a second collection-specific implementation.
- `sn_install_shennong()` now prefers unified `source` / `ref` arguments for
  GitHub and local installs while keeping `github_repo`, `github_ref`, and
  `local_path` as compatibility aliases.

### Fixed

- Fixed `sn_list_10x_paths()` so detected samples are returned in deterministic
  sample-name order instead of depending on platform-specific filesystem or
  `find` traversal order.
- Fixed `sn_standardize_gene_symbols()` so unresolved or ambiguous
  `HGNChelper` suggestions no longer propagate `NA` row names or drop
  otherwise valid original symbols. Truly missing or empty feature names are
  still removed before duplicate symbols are aggregated.
- Fixed IO edge cases where `sn_read(row_names = "column")` failed for column
  names, detected 10x spatial directories were not dispatched to a custom
  reader, `sn_write()` failed for existing `SingleCellExperiment` h5ad exports,
  and `qs2` serialization called `qs2::qs_save()` with the wrong argument name.
- Fixed `sn_run_celltypist()` path inputs so precomputed CellTypist inputs return
  a prediction table instead of trying to write metadata onto a character path.
- Fixed `sn_remove_ambient_contamination(method = "soupx")` so Seurat returns
  now add `nCount_<assay>_corrected` and `nFeature_<assay>_corrected` metadata
  from the SoupX-corrected layer, matching the decontX writeback behavior.
- Fixed the pkgdown reference index by adding the exported
  `sn_sweep_cluster_resolution()` topic.
- Fixed a malformed hidden R chunk in the clustering vignette that prevented
  pkgdown from rendering articles.
- Fixed local and CI test helpers so Seurat fixtures used in `de_enrich` and
  `utils` tests no longer depend on implicit species inference or non-returned
  normalization calls.
- Fixed pkgdown reference indexing so new helpers such as `sn_assess_qc()` and
  `sn_list_10x_paths()` are included in the generated site configuration.

# Version 0.1.2

Released 2026-03-25.

### Added

- Bundled human and mouse GENCODE gene-annotation data was added for gene-level
  filtering workflows. The snapshot now also stores genome-context fields such
  as sequence name, source, coordinates, strand, gene status/source, and
  annotation level in addition to gene identifiers and gene types.
- Small built-in PBMC example assets were added as `pbmc_small` and
  `pbmc_small_raw`, sampled from the packaged `pbmc1k` / `pbmc3k` references
  for check-safe examples and README workflows.
- Integration and cluster-diagnostics metrics were expanded with
  `sn_calculate_silhouette()`, `sn_calculate_graph_connectivity()`,
  `sn_calculate_pcr_batch()`, `sn_calculate_clustering_agreement()`,
  `sn_calculate_isolated_label_score()`, `sn_calculate_cluster_entropy()`,
  `sn_calculate_cluster_purity()`, `sn_identify_challenging_groups()`, and the
  aggregate `sn_assess_integration()` wrapper.
- Rare-cell-aware clustering support was added through
  `sn_detect_rare_cells()` and new `sn_run_cluster()` parameters that can
  append rare-aware features such as Gini-selected genes, local HVGs, local
  markers, or optional CIARA-derived features before PCA and Harmony.
- Bulk RNA-seq deconvolution support was added through
  `sn_deconvolve_bulk()`, `sn_store_deconvolution()`, and
  `sn_get_deconvolution_result()`, covering local BayesPrism runs plus
  local CIBERSORTx container workflows and result import.
- pkgdown documentation is now reorganized around workflow stages rather than
  source files alone. New end-to-end articles cover preprocessing and QC,
  clustering and integration, metrics and diagnostics, annotation and pathways,
  composition analysis, and interpretation/reporting.
- `sn_find_de()` now supports pseudobulk differential expression with
  `limma` in addition to `DESeq2` and `edgeR`.
- `sn_find_de()` now supports marker discovery with `COSGR` when the optional
  GitHub package is installed.
- `sn_initialize_project()` now scaffolds user analysis repositories from the
  shipped `inst/codex/project-template/` assets, creating a governed
  project layout with `AGENTS.md`, `memory/`, `docs/standards/`, `skills/`,
  `config/`, `data/`, `scripts/`, `notebooks/`, `runs/`, and `results/`.
- `sn_initialize_project()` is now a convenience wrapper over the packaged
  project-template initializer.
- `sn_get_codex_skill_path()` now exposes packaged Codex asset paths for the
  Codex root, package-usage skills, project template, and project-template
  skills.
- `sn_install_codex_skill()` now installs package-usage skills, project
  governance skills, or both from the packaged Codex asset layout.
- Signature catalog helpers were added:
  `sn_list_signatures()`, `sn_add_signature()`, `sn_update_signature()`, and
  `sn_delete_signature()`.
- Stored-result discovery and retrieval helpers were added:
  `sn_list_results()`, `sn_get_de_result()`, `sn_get_enrichment_result()`, and
  `sn_get_interpretation_result()`.
- The interpretation layer now supports user-supplied background context and
  dual output styles for either model-facing prompt bundles or human-readable
  summaries.

### Changed

- Internal sparse-matrix workflows were optimized for better runtime and lower
  peak memory use. Pseudobulk DE aggregation now groups columns without
  materializing dense matrices, split Seurat assay layers are combined through
  sparse triplet assembly instead of repeated indexed writes, exact kNN
  fallbacks now use blockwise distance evaluation instead of constructing a
  full cell-by-cell distance matrix, and gene-symbol standardization now reuses
  bundled annotation data plus grouped row aggregation instead of external CSV
  reads and per-column duplicate collapsing.
- Developer-facing informational notifications are now routed through internal
  package logging helpers instead of a mix of ad hoc `logger`, `message()`,
  and `cli` progress calls. Progress wording is now more consistent across
  initialization, clustering, enrichment, and preprocessing workflows, while
  real `warning()` and `stop()` conditions retain their existing semantics.
- `sn_enrich()` now uses a single `x` input with automatic dispatch across
  gene vectors, ranked named vectors, data frames, and Seurat-stored DE
  results. Its `gene_clusters` formula now drives both grouped ORA
  (`gene ~ cluster`) and ranked GSEA (`gene ~ log2fc`) workflows, and
  multi-database requests such as `database = c("H", "GOBP", "KEGG",
  "C2:CP:REACTOME")` are supported in one call.
- `sn_enrich()` now aligns MSigDB arguments with `msigdbr` by supporting
  `collection` / `subcollection`, while still accepting collection strings in
  `database` such as `H` or `C2:CP:REACTOME`.
- `sn_enrich()` now filters returned enrichment tables by raw p-value rather
  than adjusted p-value and writes one result file per requested database with
  stable `prefix` / `outdir` naming.
- Dataset documentation is now being consolidated into a single `R/data.R`
  source, and user-facing examples are shifting from simulated matrices or
  network-backed PBMC downloads to the built-in small PBMC package data.
- `sn_find_de()` now uses a single `method` argument instead of separate
  `test_use` and `pseudobulk_method` parameters.
- `sn_find_de()` now uses `layer` consistently for Seurat v5 workflows and no
  longer exposes the legacy `slot` argument.
- `sn_filter_genes()` now supports annotation-aware filtering through
  `gene_class` and exact `gene_type` values backed by the bundled GENCODE
  snapshot for human and mouse.
- Signature registry maintenance now delegates add/update/delete operations to
  the upstream `SignatuR` package API instead of maintaining a parallel custom
  tree-editing implementation inside Shennong.
- Harmony-backed integration now targets the `immunogenomics/harmony`
  `harmony2` developer branch in package metadata and CI instead of the CRAN
  release line.
- Integration metrics now prefer stored Seurat neighbor graphs when available
  and otherwise fall back to Annoy-based approximate kNN or exact distance
  search, keeping the default assessment path fast enough for routine use.
- `sn_get_signatures()` now reads from a package-owned signature snapshot built
  from the full `SignatuR` tree during development, so runtime signature
  retrieval is stable, tree-structured, and no longer depends on the installed
  `SignatuR` package version.
- Signature build assets now center on `data/shennong_signature_catalog.rda`,
  which is rebuilt directly from the upstream `SignatuR` dataset during
  development, replacing the opaque `R/sysdata.rda` storage used by the earlier
  snapshot prototype.
- `sn_enrich()` now stores enrichment results in
  `object@misc$enrichment_results[[store_name]]` when a Seurat object is
  supplied, aligning enrichment with the existing stored DE workflow.
- pkgdown articles, shipped Codex skill references, and `NEWS.md` are now
  treated as required deliverables for any user-facing workflow change.

### Fixed

- Fixed `sn_install_shennong(channel = "github")` so it no longer tries to
  resolve runtime-optional GitHub `Suggests` by default during installation.
- Removed the optional rare-cell backends `FiRE`, `CellSIUS`, and `EDGE` from
  Shennong's supported dependency surface and `sn_detect_rare_cells()`
  interface.
- Fixed `sn_plot_dot()` theme handling so the optional `catplot` theme no
  longer tries to impose an additional aspect ratio on top of `coord_fixed()`.
- Fixed Rd example line-width failures in `R CMD check` and declared the
  runtime `data.tree` dependency explicitly.

# Version 0.1.1

Released 2026-03-19.

### Added

- Automatic human/mouse species inference through `sn_get_species()` using
  `hom_genes` and mitochondrial naming patterns.
- A layer-aware pkgdown article covering inferred species, non-default count
  layers, and stored differential-expression metadata.
- Additional regression tests for species inference, DE metadata, and
  BPCells-backed Seurat layers.

### Changed

- `sn_initialize_seurat_object()` now attempts species inference before deciding
  whether to compute species-specific QC metrics.
- `sn_get_signatures()` now uses package-local fallback signatures for the core
  categories needed by Shennong workflows when `SignatuR` is unavailable.
- Stored DE results now include schema and provenance metadata such as package
  version, timestamp, assay/layer context, and threshold settings.
- `sn_plot_*()` helpers now treat `catplot` as an optional enhancement rather
  than a mandatory dependency.

### Fixed

- Fixed `sn_remove_ambient_contamination()` for BPCells-backed Seurat layers by
  materializing BPCells matrices before passing them to `decontX`.
- Fixed optional dependency handling so missing `SignatuR` no longer breaks core
  initialization and clustering paths.
- Fixed SoupX validation order so missing `raw` input is reported before package
  availability issues.
- Fixed the `sn_deconvolve_bulk()` example so the `cibersortx` dry-run path no
  longer fails package examples for missing credentials.
- Fixed GitHub Actions dependency installation so check and coverage jobs keep
  a minimal, solvable dependency set instead of trying to install unsupported
  GitHub-only optional backends during lockfile generation.

# Version 0.1.0

Released 2026-03-18.

### Added

- `sn_load_data()` as the primary example-data loader.
- Consolidated clustering into `sn_run_cluster()`.
- A unified ambient contamination interface for SoupX and decontX.
- Expanded test coverage across composition, data loading, clustering,
  utilities, visualization, and ambient contamination.
- Missing help pages for exported functions, CI scaffolding, a richer README,
  and updated package metadata for current R syntax requirements.

### Changed

- Retained `sn_load_pbmc()` as a deprecated compatibility wrapper around
  `sn_load_data()`.
- Refreshed generated documentation and package metadata as part of the
  modernization effort.

### Fixed

- Addressed multiple modernization issues in the package build, test, and
  documentation pipeline.
