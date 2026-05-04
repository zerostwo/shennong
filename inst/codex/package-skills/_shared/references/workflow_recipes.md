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

## Recipe: Normalize and cluster a single dataset

1. `sn_normalize_data()` when explicit normalization control is needed.
2. `sn_run_cluster()` for PCA, neighbors, clustering, and embeddings. Add
   `hvg_features = c(...)` when known rare-population markers should be forced
   into the PCA feature set.
3. `sn_plot_dim()` for cluster and metadata visualization.

## Recipe: Integrate multiple samples

1. Ensure batch metadata are present.
2. Run `sn_run_cluster(batch = ...)`. Set
   `normalization_method = "sctransform"` for SCTransform followed by Harmony
   integration when that normalization is desired.
3. Evaluate with `sn_assess_integration()` and metric helpers.
4. Use `sn_identify_challenging_groups()` or rare-cell helpers when needed.
5. Use `sn_transfer_labels()` when labels should be projected from a reference
   Seurat object to a query object.

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
3. Reuse later with `sn_get_de_result()` or plotting helpers.

## Recipe: Run pathway analysis

1. If enriching stored DE on a Seurat object, use `sn_enrich(x = object, source_de_name = ...)`.
2. If enriching a direct table, supply `gene_clusters = gene ~ cluster` for grouped ORA or `gene ~ log2fc` for ranked GSEA.
3. Store reusable results with `sn_store_enrichment()` when needed.

## Recipe: Build annotation evidence

1. Start from stored marker results.
2. Use `sn_prepare_annotation_evidence()`.
3. Optionally add enrichment and QC evidence.
4. For automated annotation, use `sn_interpret_annotation()`.

## Recipe: Reuse stored results

1. Discover assets with `sn_list_results()`.
2. Retrieve with `sn_get_de_result()`, `sn_get_enrichment_result()`,
   `sn_get_milo_result()`, `sn_get_deconvolution_result()`, or
   `sn_get_interpretation_result()`.
3. Avoid direct `object@misc` access in user-facing workflows.

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
