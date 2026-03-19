# Function Map

## Data and Installation

- `sn_load_data()`: load packaged PBMC example data
- `sn_load_pbmc()`: deprecated compatibility wrapper
- `sn_read()`, `sn_write()`: IO helpers
- `sn_check_version()`, `sn_install_shennong()`: package maintenance helpers
- `sn_get_codex_skill_path()`, `sn_install_codex_skill()`: bundled skill discovery and installation

## Object Creation and Preprocessing

- `sn_initialize_seurat_object()`: create a Seurat object from counts
- `sn_standardize_gene_symbols()`: harmonize gene names
- `sn_score_cell_cycle()`: add cell-cycle scores
- `sn_normalize_data()`: normalization workflows

## Quality Control

- `sn_filter_cells()`: metadata-based cell QC
- `sn_filter_genes()`: gene filtering by detection frequency
- `sn_find_doublets()`: doublet detection
- `sn_remove_ambient_contamination()`: SoupX or decontX

## Clustering, Integration, and Annotation

- `sn_run_cluster()`: single-sample clustering or Harmony integration
- `sn_run_celltypist()`: CellTypist annotation
- `sn_find_de()`: markers, direct contrasts, and pseudobulk DE

## Metrics and Summaries

- `sn_calculate_composition()`: per-group composition summaries
- `sn_calculate_lisi()`: integration quality / mixing score
- `sn_calculate_rogue()`: purity / heterogeneity score

## Visualization

- `sn_plot_dim()`: dimensionality reduction plots
- `sn_plot_feature()`: continuous feature overlays
- `sn_plot_violin()`: violin plots
- `sn_plot_dot()`: marker dot plots, including stored top markers
- `sn_plot_boxplot()`, `sn_plot_barplot()`: table-driven plots

## Signatures and Interpretation

- `sn_get_signatures()`: built-in signature/blocklist retrieval through the optional `SignatuR` package
- `sn_enrich()`: ORA and GSEA with GO or KEGG

## Utility Helpers

- `sn_check_file()`: validate file paths
- `sn_get_species()`: resolve object species
- `sn_set_path()`: create output directories
