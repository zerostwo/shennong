# Package index

## Data and IO

- [`marker_genes`](https://songqi.org/shennong/reference/marker_genes.md)
  : Pan-Immune CellTypist Metadata
- [`hom_genes`](https://songqi.org/shennong/reference/hom_genes.md) :
  Human-Mouse Homologous Genes Table
- [`shennong_gene_annotations`](https://songqi.org/shennong/reference/shennong_gene_annotations.md)
  : Bundled Shennong Gene Annotations
- [`shennong_signature_catalog`](https://songqi.org/shennong/reference/shennong_signature_catalog.md)
  : Bundled Shennong Signature Catalog
- [`sn_load_data()`](https://songqi.org/shennong/reference/sn_load_data.md)
  : Load example datasets from Zenodo
- [`sn_load_pbmc()`](https://songqi.org/shennong/reference/sn_load_pbmc.md)
  : Load PBMC example datasets
- [`sn_read()`](https://songqi.org/shennong/reference/sn_read.md)
  [`.import.rio_bpcells()`](https://songqi.org/shennong/reference/sn_read.md)
  [`.import.rio_10x()`](https://songqi.org/shennong/reference/sn_read.md)
  [`.import.rio_starsolo()`](https://songqi.org/shennong/reference/sn_read.md)
  [`.import.rio_h5ad()`](https://songqi.org/shennong/reference/sn_read.md)
  [`.import.rio_h5()`](https://songqi.org/shennong/reference/sn_read.md)
  [`.import.rio_gmt()`](https://songqi.org/shennong/reference/sn_read.md)
  : Read tabular and bioinformatics file formats
- [`sn_write()`](https://songqi.org/shennong/reference/sn_write.md)
  [`.export.rio_bpcells()`](https://songqi.org/shennong/reference/sn_write.md)
  [`.export.rio_h5ad()`](https://songqi.org/shennong/reference/sn_write.md)
  [`.export.rio_h5()`](https://songqi.org/shennong/reference/sn_write.md)
  : Write tabular and bioinformatics file formats
- [`sn_add_data_from_anndata()`](https://songqi.org/shennong/reference/sn_add_data_from_anndata.md)
  : Add metadata and embeddings exported from AnnData
- [`sn_check_version()`](https://songqi.org/shennong/reference/sn_check_version.md)
  : Check whether Shennong is up to date
- [`sn_install_shennong()`](https://songqi.org/shennong/reference/sn_install_shennong.md)
  : Install Shennong from CRAN or GitHub

## Preprocessing and QC

- [`sn_initialize_seurat_object()`](https://songqi.org/shennong/reference/sn_initialize_seurat_object.md)
  : Initialize a Seurat object with optional QC metrics
- [`sn_normalize_data()`](https://songqi.org/shennong/reference/sn_normalize_data.md)
  : Normalize data in a Seurat object
- [`sn_standardize_gene_symbols()`](https://songqi.org/shennong/reference/sn_standardize_gene_symbols.md)
  : Standardize gene symbols in a count matrix or Seurat object
- [`sn_score_cell_cycle()`](https://songqi.org/shennong/reference/sn_score_cell_cycle.md)
  : Score Cell Cycle Phases
- [`sn_filter_genes()`](https://songqi.org/shennong/reference/sn_filter_genes.md)
  : Filter genes based on the number of expressing cells
- [`sn_filter_cells()`](https://songqi.org/shennong/reference/sn_filter_cells.md)
  : Filter cells in a Seurat object based on QC metrics
- [`sn_find_doublets()`](https://songqi.org/shennong/reference/sn_find_doublets.md)
  : Find doublets using scDblFinder
- [`sn_remove_ambient_contamination()`](https://songqi.org/shennong/reference/sn_remove_ambient_contamination.md)
  : Remove ambient RNA contamination from counts

## Clustering and metrics

- [`sn_run_cluster()`](https://songqi.org/shennong/reference/sn_run_cluster.md)
  : Run clustering for a single dataset or Harmony integration workflow
- [`sn_find_de()`](https://songqi.org/shennong/reference/sn_find_de.md)
  : Run differential expression analysis on a Seurat object
- [`sn_run_celltypist()`](https://songqi.org/shennong/reference/sn_run_celltypist.md)
  : Run CellTypist for automated cell type annotation
- [`sn_calculate_lisi()`](https://songqi.org/shennong/reference/sn_calculate_lisi.md)
  : Calculate LISI score
- [`sn_calculate_rogue()`](https://songqi.org/shennong/reference/sn_calculate_rogue.md)
  : Calculate ROGUE score for Seurat Object
- [`sn_calculate_composition()`](https://songqi.org/shennong/reference/sn_calculate_composition.md)
  : Calculate Composition Proportions

## Visualization and signatures

- [`sn_plot_dim()`](https://songqi.org/shennong/reference/sn_plot_dim.md)
  : Create a dimensionality reduction plot for categorical data
- [`sn_plot_feature()`](https://songqi.org/shennong/reference/sn_plot_feature.md)
  : Plot feature expression in reduced dimensions
- [`sn_plot_violin()`](https://songqi.org/shennong/reference/sn_plot_violin.md)
  : Plot a violin plot with categorical groups
- [`sn_plot_dot()`](https://songqi.org/shennong/reference/sn_plot_dot.md)
  : Plot a dot plot with categorical groups
- [`sn_plot_boxplot()`](https://songqi.org/shennong/reference/sn_plot_boxplot.md)
  : Create a boxplot from a data frame
- [`sn_plot_barplot()`](https://songqi.org/shennong/reference/sn_plot_barplot.md)
  : Create a bar plot from a data frame
- [`show_all_palettes()`](https://songqi.org/shennong/reference/show_all_palettes.md)
  : List available color palettes
- [`sn_list_signatures()`](https://songqi.org/shennong/reference/sn_list_signatures.md)
  : List bundled Shennong signatures
- [`sn_get_signatures()`](https://songqi.org/shennong/reference/sn_get_signatures.md)
  : Retrieve bundled Shennong signature genes by category or tree path
- [`sn_add_signature()`](https://songqi.org/shennong/reference/sn_add_signature.md)
  : Add a signature to the editable source registry
- [`sn_update_signature()`](https://songqi.org/shennong/reference/sn_update_signature.md)
  : Update a signature in the editable source registry
- [`sn_delete_signature()`](https://songqi.org/shennong/reference/sn_delete_signature.md)
  : Delete a signature from the editable source registry
- [`sn_enrich()`](https://songqi.org/shennong/reference/sn_enrich.md) :
  Run gene set enrichment analysis

## Interpretation

- [`sn_list_results()`](https://songqi.org/shennong/reference/sn_list_results.md)
  : List stored Shennong analysis and interpretation results on a Seurat
  object
- [`sn_get_de_result()`](https://songqi.org/shennong/reference/sn_get_de_result.md)
  : Retrieve a stored DE result from a Seurat object
- [`sn_get_enrichment_result()`](https://songqi.org/shennong/reference/sn_get_enrichment_result.md)
  : Retrieve a stored enrichment result from a Seurat object
- [`sn_get_interpretation_result()`](https://songqi.org/shennong/reference/sn_get_interpretation_result.md)
  : Retrieve a stored interpretation result from a Seurat object
- [`sn_store_enrichment()`](https://songqi.org/shennong/reference/sn_store_enrichment.md)
  : Store an enrichment result on a Seurat object
- [`sn_prepare_annotation_evidence()`](https://songqi.org/shennong/reference/sn_prepare_annotation_evidence.md)
  : Prepare cluster-annotation evidence from a Seurat object
- [`sn_prepare_de_evidence()`](https://songqi.org/shennong/reference/sn_prepare_de_evidence.md)
  : Prepare differential-expression evidence from a stored DE result
- [`sn_prepare_enrichment_evidence()`](https://songqi.org/shennong/reference/sn_prepare_enrichment_evidence.md)
  : Prepare enrichment evidence
- [`sn_prepare_results_evidence()`](https://songqi.org/shennong/reference/sn_prepare_results_evidence.md)
  : Prepare manuscript-style results evidence
- [`sn_build_prompt()`](https://songqi.org/shennong/reference/sn_build_prompt.md)
  : Build an LLM prompt from structured Shennong evidence
- [`sn_run_llm()`](https://songqi.org/shennong/reference/sn_run_llm.md)
  : Run an LLM provider on a prepared message list
- [`sn_interpret_annotation()`](https://songqi.org/shennong/reference/sn_interpret_annotation.md)
  : Interpret cluster markers for cell-type annotation
- [`sn_interpret_de()`](https://songqi.org/shennong/reference/sn_interpret_de.md)
  : Interpret a stored differential-expression result
- [`sn_interpret_enrichment()`](https://songqi.org/shennong/reference/sn_interpret_enrichment.md)
  : Interpret a stored enrichment result
- [`sn_write_results()`](https://songqi.org/shennong/reference/sn_write_results.md)
  : Write a manuscript-style results summary from stored analysis
  outputs
- [`sn_write_figure_legend()`](https://songqi.org/shennong/reference/sn_write_figure_legend.md)
  : Write a figure legend from stored analysis outputs
- [`sn_write_presentation_summary()`](https://songqi.org/shennong/reference/sn_write_presentation_summary.md)
  : Write a presentation-style summary from stored analysis outputs

## Utilities

- [`sn_check_file()`](https://songqi.org/shennong/reference/sn_check_file.md)
  : Check if files exist
- [`sn_get_codex_skill_path()`](https://songqi.org/shennong/reference/sn_get_codex_skill_path.md)
  : Return installed Shennong Codex asset paths
- [`sn_get_species()`](https://songqi.org/shennong/reference/sn_get_species.md)
  : Retrieve or infer species information
- [`sn_initialize_project()`](https://songqi.org/shennong/reference/sn_initialize_project.md)
  : Initialize a Shennong analysis project
- [`sn_initialize_codex_project()`](https://songqi.org/shennong/reference/sn_initialize_codex_project.md)
  : Initialize Codex-style project guidance for a Shennong analysis
- [`sn_install_codex_skill()`](https://songqi.org/shennong/reference/sn_install_codex_skill.md)
  : Install the bundled Shennong Codex skill for end-user agents
- [`sn_set_path()`](https://songqi.org/shennong/reference/sn_set_path.md)
  : Create a directory path if needed
