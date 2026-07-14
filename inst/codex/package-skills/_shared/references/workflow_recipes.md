# Shennong Workflow Recipes

This reference gives compact task-oriented recipes so an agent can move from a
user request to the right Shennong function family quickly.

## Recipe: Start from raw counts or 10x outputs

1. Discover 10x inputs with `sn_list_10x_paths()` when needed.
2. Initialize the Seurat object with `sn_initialize_seurat_object()`.
3. When `sn_list_10x_paths()` returns multiple named paths, pass the whole
   vector to `sn_initialize_seurat_object(x = tenx_paths)` to import all
   samples at once as a named list of Seurat objects.
4. Infer or verify species with `sn_get_species()`.
5. Run QC and filtering with `sn_filter_cells()` and `sn_filter_genes()`.

## Recipe: Open a Shennong Data Server resource

1. Use `sn_load_data(dataset = ..., backend = "api")` for a lazy
   ShennongData 0.2 resource handle.
2. Select an assay or layer with
   `api_args = list(assay = ..., layer = ...)` when the default resource view
   is not sufficient.
3. Keep `lazy = TRUE` for query planning. Set `lazy = FALSE` and supply
   `api_args = list(collect_args = list(...))` only when materialization is
   explicitly required.
4. Use `backend = "auto"` for server-first discovery with local/Zenodo
   fallback, or `backend = "local"` to prohibit server access.

## Recipe: Normalize and cluster a single dataset

1. `sn_normalize_data()` when explicit normalization control is needed.
2. `sn_run_cluster()` for PCA, neighbors, clustering, and embeddings. Add
   `hvg_features = c(...)` when known rare-population markers should be forced
   into the PCA feature set, and use `block_genes = c(...)` when cell-cycle,
   ribosomal, mitochondrial, or other signature/custom genes should be excluded
   from internally selected HVGs in log-normalization or SCTransform workflows.
   Re-run the returned object with a new
   `resolution` to reuse normalization, feature selection, PCA, neighbors, and
   UMAP while recomputing cluster labels.
3. `sn_plot_dim()` for clusters and metadata visualization.

## Recipe: Integrate multiple samples

1. Ensure batch metadata are present.
2. Run `sn_run_cluster(batch = ..., integration_method = "harmony")` for the
   default fast workflow. Use `integration_method = "coralysis"` for
   Coralysis multi-level integration on imbalanced datasets, or
   `"seurat_cca"` / `"seurat_rpca"` to compare Seurat layer-integration
   backends. Use `"scvi"` or `"scanvi"` when the workflow should run through a
   pixi-managed scverse environment under `~/.shennong/pixi/`; scANVI requires
   `integration_control = list(label_by = ...)`. Use `sn_pixi_paths()` to
   inspect where Shennong will create the pixi workspace and
   `sn_list_pixi_environments()` / `sn_pixi_config_path()` to inspect bundled
   configs under `inst/pixi/`. Use
	   `integration_control = list(accelerator = "auto", mirror = "auto")` when
	   CPU/CUDA selection and Shennong-level mirror configuration should be handled
	   automatically. Set
	   `normalization_method = "sctransform"` only with Harmony when that
	   normalization is desired.
	   For CITE-seq, use `modality = "cite_seq"` and select
	   `multimodal_method = "wnn"`, `"totalvi"`, `"coralysis"`, or `"mmochi"`
	   depending on whether the protein signal should enter through Seurat WNN,
	   scvi-tools totalVI, Coralysis on the ADT assay, or MMoCHi landmark
	   registration on the ADT assay. MMoCHi can run as a single-sample CITE-seq
	   workflow with `batch = NULL`.
3. Evaluate with `sn_assess_integration()` and metric helpers.
4. Use `sn_calculate_variance_explained(variables = c(...))` to rank platform,
   study, tissue, sample, or other metadata drivers of residual embedding
   variation.
5. If clusters are stable but the two-dimensional display is hard to read, use
   `umap_control = list(n.neighbors = ..., min.dist = ..., spread = ...)` and
   `rerun_from = "umap"` to retune the UMAP without changing integration or
   clustering.
6. Use `sn_identify_challenging_groups()` or rare-cell helpers when needed.
7. Use `sn_transfer_labels()` when labels should be projected from a
   reference to a query object. Use the default Seurat backend for anchor
   transfer, or `method = "coralysis"` when the reference was trained and
   stored with Coralysis. Use `method = "scanvi"` or `method = "scarches"`
   for semi-supervised scVI-family label transfer.
8. Use `sn_prepare_label_transfer_reference()` before saving a reusable
   reference. This avoids shipping a full analysis Seurat object when label
   transfer only needs a compact backend-specific reference.

## Recipe: Simulate single-cell data

1. Start from a clustered Seurat object or a SingleCellExperiment with
   covariates in metadata/colData.
2. Use `sn_simulate(method = "scdesign3")` with a `celltype`, `pseudotime`, or
   spatial design column.
3. Return `"seurat"` for a ready-to-use simulated object or `"result"` when
   scDesign3 model diagnostics are needed.

## Recipe: Find markers or run a contrast

1. Use `sn_find_de()`.
2. Prefer `return_object = TRUE` when the result should be stored.
3. Use `sn_annotate_de_features()` when the next step is to prioritize
   transcription factors, surface/plasma-membrane genes, cytokines, or
   chemokines from the marker/DE table.
4. Reuse later with `sn_get_de_result()` or plotting helpers.

## Recipe: Run pathway analysis

1. If enriching stored DE on a Seurat object, use `sn_enrich(x = object, source_de_name = ...)`.
2. If enriching a direct table, supply `gene_clusters = gene ~ cluster` for grouped ORA or `gene ~ log2fc` for ranked GSEA.
3. Store reusable results with `sn_store_enrichment()` when needed.

## Recipe: Analyze standalone bulk transcriptomics

1. Supply a feature-by-sample matrix plus row-named sample metadata, a list
   containing those objects, or a `SummarizedExperiment`.
2. Run `sn_assess_bulk_qc()` first and review `tables$samples`,
   `embeddings$pca`, and `tables$correlation`; do not remove a sample solely
   because one automatic outlier flag is true.
3. Use `sn_find_bulk_de(design = ..., contrast = c(variable, numerator,
   denominator))`. Keep `method = "auto"` unless the statistical backend is
   prespecified: integer counts choose edgeR, continuous expression chooses
   limma, and mixed-effects designs choose dream.
4. Use `sn_score_bulk_pathways()` for sample pathway scores and inspect
   `tables$coverage`. Use `sn_run_wgcna()` for modules/eigengenes and
   sample-level trait associations.
5. Use `sn_run_survival()` for adjusted Cox models and
   `sn_run_clinical_association()` for numeric or categorical phenotype tests.
   Retrieve evidence from each returned result's named `tables` rather than
   recomputing statistics inside plotting code.

## Recipe: Score and compare gene programs

1. Use `sn_score_programs(method = "ucell")` for sparse per-cell scoring.
2. Use `method = "gsva"` or `"ssgsea"` with `group_by = "sample"` for
   aggregated profiles; do not force large sparse cell matrices through a
   dense pathway-scoring backend.
3. Inspect `result$tables$coverage` before interpreting a program with poor
   feature overlap.
4. Use `sn_test_programs(sample_by = "patient", condition_by = ...)` so
   biological samples, not cells, are the inferential units.
5. Retrieve the score and comparison results with `sn_get_result()`.

## Recipe: Discover programs and infer regulons

1. Use `sn_discover_programs(method = "nmf")` for a lightweight local
   discovery run and inspect both `tables$fit_diagnostics` and
   `tables$stability` before naming programs.
2. Supply `group_by` only when each stratum contains enough cells and features
   for an independent factorization. cNMF and Hotspot results enter through an
   explicit `backend_control$runner` or `backend_control$result` adapter.
3. Use `sn_run_grn(method = "genie3", regulators = ...)` for direct R
   inference. Use pySCENIC through an explicit adapter when motif pruning and
   AUCell activity are required; keep motif-database provenance with the
   adapter result.
4. Retrieve `sn_get_result(object, "grn", name)` and inspect
   `tables$edges`, `tables$regulons`, `tables$activity`, and
   `tables$specificity` before plotting with `sn_plot_regulon()`.

## Recipe: Summarize composition shifts

1. Use `sn_calculate_composition(group_by = ..., variable = ..., measure = "both")` for grouped counts and percentages.
2. Use `sn_calculate_roe(group_by = ..., variable = ...)` when the question is categorical observed-over-expected enrichment; set `return_matrix = TRUE` for heatmap inputs.
3. Use `sn_compare_composition(sample_by = ..., group_by = ..., variable = ..., contrast = ...)` for replicate-aware condition comparisons.

## Recipe: Run communication or regulatory activity analysis

1. Use `sn_run_cell_communication(method = c("liana", "cellchat"),
   consensus = TRUE)` for cross-method ranking, `method = "nichenet"` for
   ligand-target inference, or `method = "multinichenet"` with `sample_by`,
   `condition_by`, and a two-level `contrast` for official sample-aware
   differential communication.
2. Supply `sample_by` for any backend when a condition claim should be based on
   biological-sample LR expression rather than cell counts. Inspect
   `result$tables$sample_evidence`, `condition_comparison`,
   `method_concordance`, and `ligand_targets`.
3. Use `sn_plot_communication()`, `sn_plot_ligand_target()`, and
   `sn_plot_communication_comparison()` for result-aware figures.
4. Use `sn_run_regulatory_activity(method = "dorothea")` for TF activity and
   `method = "progeny"` for pathway activity.
5. Retrieve stored outputs with `sn_get_cell_communication_result()` or
   `sn_get_regulatory_activity_result()`.

## Recipe: Infer CNV and compare metabolic activity

1. Define trustworthy normal cells explicitly, then run
   `sn_run_cnv(reference_cells = ...)` or use `reference_by` plus
   `reference_cat`; do not infer malignancy without an auditable reference.
2. Supply `sample_by` for multi-patient tumors. Review `tables$primary`,
   `chromosome`, `sample_summary`, and `expression_association` before using a
   malignant call or subclone in downstream figures.
3. Use `sn_plot_cnv()` for the chromosome, CNV UMAP, malignancy, sample, and
   expression-association views. Retrieve with
   `sn_get_result(object, "cnv", store_name)`.
4. Start metabolism with `sn_run_metabolism(scoring_method = "ucell")` and
   curated signatures. Supply `sample_by`, `condition_by`, and optionally
   `group_by` so differential activity is estimated from samples, not cells.
5. Inspect signature coverage and sample scores. Use explicit
   `backend_control$runner` or `backend_control$result` when standardizing
   scFEA/Compass output, then plot with `sn_plot_metabolism()`.

## Recipe: Build annotation evidence

1. Use `sn_run_annotation(method = "consensus")` for marker-only annotation,
   or supply `reference`, `reference_label_by`, and a reference backend when a
   biologically matched atlas is available.
2. Discover and retrieve the result with
   `sn_list_results(object, type = "annotation")` and
   `sn_get_result(object, "annotation", store_name)`.
3. Review `sn_review_annotation()` plus the confidence and marker plots before
   accepting low-margin labels.
4. Use `sn_prepare_annotation_evidence()` when stored DE/enrichment/QC evidence
   is needed for narrative interpretation.
5. `sn_interpret_annotation()` may explain or rank evidence, but must not
   overwrite the computational labels stored by `sn_run_annotation()`.

## Recipe: Infer a trajectory and test dynamic genes

1. Choose a biologically defensible reduction and cluster label, then call
   `sn_run_trajectory(start = ..., end = ...)`; do not infer orientation from
   an unlabeled embedding alone.
2. Review `sn_plot_trajectory()` and the stored `graphs$lineages`,
   `tables$terminal_states`, and warnings before interpreting pseudotime.
3. Inspect per-lineage pseudotime and probabilities with
   `sn_plot_pseudotime()` and `sn_plot_lineage_probability()` rather than using
   only the primary-lineage metadata shortcut.
4. For a formal tradeSeq analysis, supply `dynamic_features`, verify
   `tables$convergence`, and use BH-adjusted values from `dynamic_genes` and
   `branch_genes`.
5. Use `sn_plot_dynamic_heatmap()`, `sn_plot_gene_trend()`, and
   `sn_plot_branch_comparison()` for result-backed review; retrieve the durable
   result with `sn_get_result(object, "trajectory", store_name)`.

## Recipe: Estimate RNA velocity and fate

1. Keep raw spliced and unspliced counts in explicit Seurat layers and verify
   cell/feature overlap before calling `sn_run_velocity()`.
2. Review velocity confidence, projected vectors, and transition edges. A
   visually smooth arrow field is not sufficient evidence by itself.
3. Run `sn_run_fate(velocity_name = ...)` so CellRank consumes the retained
   scVelo H5AD transition evidence. Use the stability terminal-state rule by
   default; set `terminal_method = "top_n"` or `terminal_states` only with a
   documented biological rationale.
4. Retrieve `velocity` and `fate` results separately with `sn_get_result()`;
   plot them with `sn_plot_velocity()` and `sn_plot_fate()`.

## Recipe: Test abundance and prioritize states

1. Use `sn_test_abundance(sample_by = ..., condition_by = ...,
   cell_type_by = ...)`; every condition and design covariate must be constant
   within a biological sample.
2. Start with Propeller, validate important effects with sample-label
   permutation, and use Milo only when neighborhood-level resolution is needed.
3. Audit completed zero-count rows, sample proportions, contributions, and
   adjusted significance in the stored result before plotting.
4. Use `sn_prioritize_states(method = "augur")` for state-wise condition
   separability, with `sample_by` so entire samples are held out.
5. Use RareQ when topology-supported rare states are the target. Use Scissor
   only with explicit gene-by-bulk-sample expression and corresponding bulk
   phenotype; cell metadata is not a valid replacement.

## Recipe: Run a spatial analysis

1. Put two finite coordinate columns in metadata or pass `spatial_cols`.
2. Start with `sn_find_spatial_features(method = "morans_i")` and inspect the
   permutation null; use nnSVG when its Gaussian-process model is required.
3. Use `sn_find_spatial_domains(method = "banksy")` when BANKSY is installed,
   or pass an explicit result from stLearn, BayesSpace, or CellCharter.
4. Use `sn_run_spatial_neighborhood(group_by = ...)` to retain the graph,
   permutation enrichment, and distance-bin co-occurrence together.
5. Run ordinary communication first, then call
   `sn_run_spatial_communication()` to add proximity evidence. Do not interpret
   proximity as ligand-receptor evidence by itself.

## Recipe: Reuse stored results

1. Discover assets with `sn_list_results()`.
2. Retrieve with `sn_get_de_result()`, `sn_get_enrichment_result()`,
   `sn_get_milo_result()`, `sn_get_deconvolution_result()`,
   `sn_get_cell_communication_result()`,
   `sn_get_regulatory_activity_result()`, or `sn_get_interpretation_result()`.
3. Use `sn_download_zenodo(record_id = ..., files = ...)` when a reusable
   public Zenodo record should be pulled into a local cache. Public records do
   not need a token; pass `token = ...` only for restricted/private files.
4. Use `sn_list_datasets()` to inspect the current Shennong public collection,
   then `sn_load_data(dataset = <sample_id>)` or
   `sn_load_data(dataset = <study_id>, sample_id = <sample_id>)` to cache the
   study ZIP and extract one sample's filtered/raw H5 or metrics file.
5. Use `sn_upload_zenodo(files = ..., title = ..., version = ...)` when a
   reusable object, marker table, reference, or manifest should be released for
   cross-project reuse with stable checksums and Zenodo version metadata.
6. Avoid direct `object@misc` access in user-facing workflows.

## Recipe: LLM-backed interpretation

1. Build evidence with `sn_prepare_*_evidence()`.
2. Build prompt bundles with `sn_build_prompt()` when the caller wants a prompt only.
3. Create or resolve a provider with `sn_make_ellmer_provider()`.
4. Run `sn_interpret_annotation()`, `sn_interpret_de()`,
   `sn_interpret_enrichment()`, `sn_write_results()`,
   `sn_write_figure_legend()`, or `sn_write_presentation_summary()`.

## Recipe: Bulk deconvolution from single-cell reference

1. Prepare or choose reference labels in the Seurat object.
2. Run `sn_deconvolve_bulk()`.
3. Store results with `sn_store_deconvolution()`.
4. Retrieve later with `sn_get_deconvolution_result()`.

## Recipe: Initialize a governed project

1. Use `sn_initialize_project(with_agent = TRUE)`.
2. Follow the generated `AGENTS.md`, `memory/`, and `skills/`.
3. Keep package usage and project governance distinct.

## Common agent rules

- Prefer one Shennong entry point over reconstructing a workflow from raw Seurat calls.
- Prefer stored-result reuse over recomputation when the object already contains the needed analysis product.
- Keep the main object as Seurat unless a function explicitly returns a table, matrix, or prompt bundle.
- When a user asks for interpretation, separate evidence preparation from inference.
- Use current Shennong public argument names only. Do not revive retired
  compatibility aliases such as `group_col`, `sample_col`, `annotation_col`,
  `condition_col`, `cluster_col`, `label_col`, `labels_key`, `groupby`,
  `cnv_score_groupby`, `slot`, `angle`, `query`, `github_repo`,
  `github_ref`, or `local_path`.
- Use `sn_call_*()` for direct managed-Python command execution. Use
  object-level `sn_run_*()` Python wrappers only when a Seurat object is being
  exported to and imported back from the backend workflow.
