---
name: use-shennong-single-cell-workflows
description: Use when working with current Shennong single-cell, CITE-seq, spatial, deconvolution, communication, CNV, metabolism, regulatory, simulation, visualization, IO, runtime, or reporting workflows built around Seurat objects.
---

# use-shennong-single-cell-workflows

## Purpose

Teach the agent how to use the current Shennong API for preprocessing,
clustering/integration, CITE-seq, label transfer, differential expression,
enrichment, composition, Milo, bulk deconvolution, communication, CNV, metabolism, regulatory
activity, simulation, visualization, IO, Python runtime helpers, and reporting.
This skill is the main entry point for package usage.

## When To Use

- building or extending an analysis workflow with the Shennong API
- replacing raw Seurat boilerplate with Shennong helpers

## Required Inputs

- counts or a Seurat object
- metadata columns
- intended workflow stage

## Required Outputs

- a Shennong-driven workflow step
- reusable Seurat outputs or stored results
- a Zenodo upload plan or record when the user asks to release reusable data

## Rules

- prefer Shennong exported APIs when they already expose the needed capability
- keep the main object as a Seurat object unless another return type is required
- use `sn_upload_zenodo()` for reusable data releases instead of hand-written
  `zen4R` upload code; keep `publish = FALSE` until the user explicitly wants
  to publish the Zenodo draft
- use `sn_download_zenodo()` for reusable data downloads from public Zenodo
  records; do not require a token unless the record is restricted/private
- use `sn_list_datasets()` before loading the Shennong public Zenodo
  collection, then call `sn_load_data(dataset = <sample_id>)` or
  `sn_load_data(dataset = <study_id>, sample_id = <sample_id>)`
- use vectorized `sn_load_data(dataset = c(...))` when several Shennong
  example datasets should be loaded together; filtered inputs return a merged
  Seurat object and raw inputs return a named list
- use `sn_load_data(dataset = ..., backend = "api")` for lazy Shennong Data
  Server resources; select assay/layer views through `api_args`, and set
  `lazy = FALSE` only when the user explicitly wants materialization
- respect strict `sn_verb_noun` naming conventions
- use the shared package API and workflow references instead of inventing
  partial wrapper logic
- use `sn_find_de()` as the common single-cell, pseudobulk, and standalone bulk
  DE entry point; keep `sn_find_bulk_de()` only for compatibility with older
  scripts
- use current public argument names only; do not pass retired compatibility
  names such as `group_col`, `sample_col`, `annotation_col`, `condition_col`,
  `cluster_col`, `label_col`, `labels_key`, `groupby`, `cnv_score_groupby`,
  scalar rare-feature threshold aliases, `slot`, `angle`, `query`, `github_repo`,
  `github_ref`, or `local_path`
- use `batch` for `sn_run_cluster()`, `sn_run_scvi()`, and `sn_run_scanvi()`;
  use the documented `*_by` selectors for functions whose current formals use
  them, such as composition, metrics, scArches, scPoli, and label transfer
- use `layer` for Shennong expression-layer selectors; only internal calls to
  backend packages should use backend-specific names such as Seurat's `slot`
- call `sn_call_*()` helpers for direct managed-Python commands; reserve
  `sn_run_*()` Python wrappers for object-level workflows that export/import a
  Seurat object
- if work is happening inside an initialized project, also respect the project
  `AGENTS.md`, `memory/`, and `docs/standards/`

## References

- `../_shared/references/package_api_map.md`
- `../_shared/references/workflow_recipes.md`

## Procedure

1. Start with data access and preprocessing:
   initialize or load the object, infer species if needed, and run QC with
   `sn_filter_cells()` and `sn_filter_genes()`.
   `sn_filter_genes()` can combine `min_cells` with bundled GENCODE-based
   `gene_class` or exact `gene_type` filtering for human and mouse workflows.
2. Run clustering or batch integration with `sn_run_cluster()`. Use
   `hvg_features = c(...)` when the user has known marker genes that should be
   forced into the backend feature set, use
   `rare_feature_method = "gini"` or `"local_markers"` when Shennong should
   automatically add rare-aware genes, use `block_genes = c(...)` to exclude
   bundled signature queries such as `cellCycle.G2M`, `ribo`, or `mito` and/or
   custom gene symbols from internally selected HVGs in log-normalization or
   SCTransform workflows, and use
	   `integration_method = "harmony"`, `"coralysis"`, `"seurat_cca"`,
	   `"seurat_rpca"`, `"scvi"`, or `"scanvi"` when a specific
	   batch-integration backend is requested. For scVI/scANVI, Shennong manages a
   shared pixi scverse project under `~/.shennong/pixi/scvi/`, writes run
   artifacts under `~/.shennong/runs/`, and imports the latent reduction back
   into Seurat; scANVI requires `integration_control = list(label_by = ...)`.
   Coralysis stores the trained reference SingleCellExperiment by default for
   label transfer; use `integration_control = list(store_sce = FALSE)` only for
   clustering-only runs.
   Use `sn_pixi_paths()` when users ask where Python environments live, use
   `sn_list_pixi_environments()` and `sn_pixi_config_path()` to inspect bundled
   configs under `inst/pixi/`, and pass
   `integration_control = list(accelerator = "auto", mirror = "auto")` when
   GPU/CPU selection and China-friendly mirror configuration should be handled
   by Shennong.
   Re-running `sn_run_cluster()` on its own output reuses matching stages by
   default; use `rerun_from = "integration"` or `reuse = FALSE` when a stage
	   must be forced to recompute. Leiden clustering auto-installs `leidenbase` by
	   default unless `auto_install = FALSE`.
	   Use `umap_control = list(n.neighbors = ..., min.dist = ..., spread = ...)`
	   with `rerun_from = "umap"` when only the two-dimensional embedding needs
	   retuning.
	   Use `normalization_method = "sctransform"` with `batch = ...` only when an
	   SCTransform-normalized Harmony integration is requested.
	   For CITE-seq objects with paired RNA and ADT assays, use
	   `modality = "cite_seq"` plus `multimodal_method = "wnn"` for Seurat
	   weighted nearest-neighbor clustering on `weighted.nn` / `wsnn`,
	   `multimodal_method = "totalvi"` for scvi-tools totalVI RNA+ADT latent
	   integration, or `multimodal_method = "coralysis"` to run native
	   Coralysis on the ADT protein assay. Use `multimodal_method = "mmochi"`
	   when ADT alignment should be driven by MMoCHi landmark registration
	   before protein-only clustering; it can run with `batch = NULL` for a
	   single CITE-seq sample.
3. Assess integration quality or cluster_by structure with
   `sn_assess_integration()`, `sn_calculate_lisi()`,
   `sn_calculate_variance_explained()`,
   `sn_calculate_isolated_label_score()`, or
   `sn_identify_challenging_groups()` when sample mixing, rare groups,
   isolated labels, or difficult-to-separate populations matter.
4. Move to downstream biological interpretation:
   marker discovery with `sn_find_de()`, pathway analysis with `sn_enrich()`,
   marker-class prioritization with `sn_annotate_de_features()` for TF,
   surface/plasma-membrane, cytokine, and chemokine hits,
   traceable annotation with `sn_run_annotation()`; use its default consensus
   for marker evidence with optional SingleR/reference support, then inspect
   low-confidence cells/clusters with `sn_review_annotation()` and retrieve the
   stored result with `sn_get_result(object, "annotation", name)`. Use
   `sn_map_cell_ontology()` for explicit ontology mapping. Lower-level
   reference annotation remains available with `sn_transfer_labels()` or
   `sn_transfer_labels(method = "coralysis")`; use
   `sn_transfer_labels(method = "scanvi")` or
   `sn_transfer_labels(method = "scarches")` when semi-supervised scVI-family
   label transfer is requested. Use
   `sn_prepare_label_transfer_reference()` before saving a reusable reference
   for later transfer. Optional external
   annotation with `sn_run_celltypist()`.
   Score named or bundled gene programs with `sn_score_programs()`; use UCell
   for sparse per-cell scoring, GSVA/ssGSEA with `group_by` for aggregated
   profiles, and `sn_test_programs(sample_by = ...)` for replicate-aware
   condition tests.
   Discover latent programs with `sn_discover_programs()` and review restart
   diagnostics before interpreting NMF factors. Infer regulatory networks with
   `sn_run_grn()`; GENIE3 is directly runnable, while pySCENIC, legacy SCENIC,
   GRNBoost2, cNMF, and Hotspot require explicit runner/result adapters so
   external runtime and database provenance remain visible.
   Infer cluster-aware lineage structure with `sn_run_trajectory()`; provide
   explicit `start`/`end` cluster labels, retrieve the complete result with
   `sn_get_result(object, "trajectory", store_name)`, and inspect per-lineage
   pseudotime/probabilities before using tradeSeq dynamic or branch tables.
   Set `dynamic_features` to an explicit auditable feature set for formal
   analyses, and use `test_dynamic = FALSE` only for topology review.
   Run `sn_run_velocity()` only when raw spliced and unspliced layers are
   present. Select `method = "regvelo"` only with an explicit, versioned
   regulator-target GRN in `backend_control$prior_grn`; otherwise use scVelo.
   Review stored transition/confidence evidence before passing the retained
   H5AD artifact to `sn_run_fate()` for CellRank GPCCA terminal states and fate
   probabilities.
   Test cell-type abundance with `sn_test_abundance()` using `sample_by` as the
   biological replicate; use Propeller by default, permutation for transparent
   validation, and Milo for neighborhood effects. Use
   `sn_prioritize_states()` for perturbation separability, RareQ topology, or
   Scissor only when matched bulk expression and phenotype are supplied.
   Prefer `gene_clusters` formulas such as `gene ~ cluster` for grouped ORA
   or `gene ~ log2fc` for ranked GSEA, and use `database = c(...)` when the
   same input should be tested against multiple databases in one call.
5. Use `sn_calculate_composition()`, `sn_calculate_roe()`,
   `sn_compare_composition()`, `sn_run_milo()`, `sn_plot_composition()`, and
   `sn_deconvolve_bulk()` for
   comparative summaries across samples, conditions, annotations, or paired
   bulk RNA-seq mixtures.
6. Use `sn_run_cell_communication()` for LIANA, CellChat, CellPhoneDB,
   NicheNet, MultiNicheNet, or cross-method consensus. Supply `sample_by` and
   `condition_by` for replicate-aware comparisons, then inspect the stored
   concordance, sample evidence, condition effects, and ligand-target tables.
   Use `sn_run_regulatory_activity()` for DoRothEA TF or PROGENy pathway
   activity workflows.
7. Use `sn_run_cnv()` only with explicit normal references; include
   `sample_by` for multi-patient data and review chromosome evidence,
   malignancy scores, subclones, and sample summaries with `sn_plot_cnv()`.
   Use `sn_run_metabolism()` with UCell by default and `sample_by` before any
   condition claim; scFEA/Compass require an explicit runner or parsed result.
   For spatial data, preserve coordinate metadata and use `sn_run_spatial()`
   or the explicit feature/domain/neighborhood functions. Run communication
   inference before adding distance constraints; proximity is supporting
   evidence, not a substitute interaction score.
8. Use `sn_plot_dim()`, `sn_plot_feature()`, `sn_plot_heatmap()`,
   `sn_plot_violin()`, `sn_plot_dot()`, `sn_plot_boxplot()`,
   `sn_plot_barplot()`, `sn_plot_composition()`, and `sn_plot_milo()` for
   package-style plots; resolve reusable colors with `sn_list_palettes()` and
   `sn_get_palette()`.
9. Build prompts or stored-result summaries with the interpretation helpers
   when a narrative or report-ready output is needed.
10. Simulate from real objects with `sn_simulate(method = "scdesign3")`, or
   `sn_simulate_scdesign3()` when the task needs direct scDesign3 controls.
11. Inspect and reuse bundled signatures with `sn_list_signatures()` and
   `sn_get_signatures()` when workflows need curated blocklists or marker
   programs.
12. Use `sn_check_pixi()`, `sn_ensure_pixi()`, `sn_pixi_paths()`,
   `sn_list_pixi_environments()`, `sn_pixi_config_path()`,
   `sn_prepare_pixi_environment()`, `sn_call_pixi_environment()`, and the
   family-specific `sn_call_*()` helpers when a Python backend must be checked,
   prepared, or invoked directly.
13. Use `sn_check_version()`, `sn_install_shennong()`,
   `sn_list_dependencies()`, and `sn_install_dependencies()` for package
   maintenance tasks.
14. When the correct entry point is unclear, read
   `../_shared/references/package_api_map.md` and choose the exported `sn_*`
   function that matches the task instead of falling back to raw Seurat calls.

## Common Mistakes

- bypassing stored-result helpers
- mixing raw Seurat conventions with Shennong naming in ways that obscure intent
- failing to record assay or layer assumptions
- hardcoding gene blocklists when a shipped signature path already exists

## Examples

- `sn_initialize_seurat_object()`
- `sn_run_cluster()`
- `sn_transfer_labels()`
- `sn_run_annotation()`
- `sn_review_annotation()`
- `sn_plot_annotation_confidence()`
- `sn_plot_annotation_markers()`
- `sn_simulate(method = "scdesign3")`
- `sn_assess_integration()`
- `sn_calculate_variance_explained()`
- `sn_calculate_isolated_label_score()`
- `sn_identify_challenging_groups()`
- `sn_find_de(..., return_object = TRUE)`
- `sn_annotate_de_features(object, de_name = "cluster_markers")`
- `sn_enrich(x = object, source_de_name = "cluster_markers")`
- `sn_score_programs(object, signatures, method = "ucell")`
- `sn_test_programs(object, score_name, condition_by, sample_by)`
- `sn_plot_program_activity()` / `sn_plot_program_heatmap()`
- `sn_calculate_composition()`
- `sn_calculate_roe()`
- `sn_run_milo()`
- `sn_run_cell_communication(method = "cellchat")`
- `sn_run_regulatory_activity(method = "dorothea")`
- `sn_store_cell_communication()` / `sn_get_cell_communication_result()`
- `sn_store_regulatory_activity()` / `sn_get_regulatory_activity_result()`
- `sn_plot_dim()`
- `sn_plot_feature()`
- `sn_plot_heatmap()`
- `sn_deconvolve_bulk(..., method = "cibersortx", cibersortx_dry_run = TRUE)`
- `sn_store_deconvolution()` / `sn_get_deconvolution_result()`
- `sn_check_pixi()` / `sn_call_scvi()`
- `sn_read()` / `sn_write()`
- `sn_list_signatures(species = "human")`
- `sn_get_signatures(species = "human", category = "Compartments/Mito")`
