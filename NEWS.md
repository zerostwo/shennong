# Shennong 0.1.0 (2023-11-12)

## New Features

-   `sn_run_seurat()`: Implements a comprehensive workflow to run Seurat analyses. This function is a high-level wrapper that integrates various steps of Seurat analysis into a single, streamlined function.

-   `sn_run_seurat_standard()`: Implements a standard Seurat workflow.

-   `sn_run_seurat_sctransform()`: This function is designed to normalize single-cell RNA sequencing data using the SCTransform method, providing more accurate and reliable results.

-   `sn_run_seurat_postprocessing()`: Executes post-processing steps for Seurat workflow including neighbor finding, clustering, and UMAP/t-SNE reduction.

-   `sn_check_package()`: This function check if packages are installed and load them.
