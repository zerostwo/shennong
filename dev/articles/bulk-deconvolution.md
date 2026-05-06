# Bulk deconvolution from a PBMC3k reference

Bulk deconvolution starts with one clear contract: a single-cell
reference, labels on the reference cells, and a bulk matrix with
overlapping genes. Shennong wraps that contract for CIBERSORTx and
BayesPrism while keeping credentialed or containerized steps explicit.

This article uses PBMC3k clusters as a teaching reference.

## Build a small reference and mock bulk matrix

``` r

library(Shennong)
library(Seurat)
library(dplyr)

pbmc <- sn_load_data("pbmc3k")
#> INFO [2026-05-05 23:45:37] Initializing Seurat object for project: pbmc3k.
#> INFO [2026-05-05 23:45:37] Running QC metrics for human.
#> INFO [2026-05-05 23:45:37] Seurat object initialization complete.

pbmc <- sn_run_cluster(
  object = pbmc,
  normalization_method = "seurat",
  nfeatures = 1500,
  dims = 1:15,
  resolution = 0.6,
  species = "human",
  verbose = FALSE
)

pbmc$cell_type <- paste0("cluster_", pbmc$seurat_clusters)

counts <- SeuratObject::LayerData(pbmc, assay = "RNA", layer = "counts")

bulk <- cbind(
  bulk_a = Matrix::rowSums(counts[, pbmc$cell_type %in% c("cluster_0", "cluster_1"), drop = FALSE]),
  bulk_b = Matrix::rowSums(counts[, pbmc$cell_type %in% c("cluster_2", "cluster_3"), drop = FALSE])
)
```

The mock bulk samples are sums of selected single-cell clusters. Real
bulk data should be normalized and gene-aligned according to the backend
requirements.

## Prepare a CIBERSORTx run without launching the container

Use `cibersortx_dry_run = TRUE` to validate exported files and the
command bundle locally. This is the right mode for examples, tests, and
documentation.

``` r

bundle <- sn_deconvolve_bulk(
  x = pbmc,
  bulk = bulk,
  method = "cibersortx",
  cell_type_by = "cell_type",
  outdir = file.path(tempdir(), "pbmc3k-cibersortx"),
  prefix = "pbmc3k_demo",
  cibersortx_email = "demo@example.org",
  cibersortx_token = "fake-token",
  cibersortx_dry_run = TRUE,
  return_object = FALSE
)

names(bundle)
#> [1] "table"     "files"     "artifacts" "method"
bundle$command
#> NULL
bundle$files
#> $single_cell_reference
#> [1] "/tmp/RtmpSHYZeJ/pbmc3k-cibersortx/sample_file_for_cibersort.txt"
#> 
#> $mixture
#> [1] "/tmp/RtmpSHYZeJ/pbmc3k-cibersortx/mixture_file_for_cibersort.txt"
#> 
#> $signature_matrix
#> [1] "/tmp/RtmpSHYZeJ/pbmc3k-cibersortx/CIBERSORTx_sample_file_for_cibersort_inferred_phenoclasses.CIBERSORTx_sample_file_for_cibersort_inferred_refsample.bm.K999.txt"
#> 
#> $result
#> [1] "/tmp/RtmpSHYZeJ/pbmc3k-cibersortx/CIBERSORTx_pbmc3k_demo_Results.txt"
```

For an actual CIBERSORTx run, store credentials locally once and omit
the placeholder values.

``` r

sn_set_cibersortx_credentials(
  email = Sys.getenv("CIBERSORTX_EMAIL"),
  token = Sys.getenv("CIBERSORTX_TOKEN")
)

result <- sn_deconvolve_bulk(
  x = pbmc,
  bulk = bulk,
  method = "cibersortx",
  cell_type_by = "cell_type",
  outdir = "results/deconvolution/cibersortx",
  prefix = "pbmc3k"
)
```

## Import or store completed fractions

When CIBERSORTx has already produced a fraction table, import it through
`cibersortx_result` or store the table explicitly. Keeping results on
the reference object makes downstream retrieval predictable.

``` r

fractions <- data.frame(
  sample = c("bulk_a", "bulk_b"),
  cluster_0 = c(0.45, 0.10),
  cluster_1 = c(0.35, 0.05),
  cluster_2 = c(0.10, 0.50),
  cluster_3 = c(0.10, 0.35),
  check.names = FALSE
)

pbmc <- sn_store_deconvolution(
  object = pbmc,
  result = fractions,
  store_name = "pbmc3k_mock_bulk",
  method = "cibersortx",
  bulk_samples = fractions$sample,
  reference_label = "cell_type",
  artifacts = bundle$files,
  return_object = TRUE
)

sn_get_deconvolution_result(
  pbmc,
  deconvolution_name = "pbmc3k_mock_bulk"
)
#> # A tibble: 2 × 5
#>   sample cluster_0 cluster_1 cluster_2 cluster_3
#>   <chr>      <dbl>     <dbl>     <dbl>     <dbl>
#> 1 bulk_a      0.45      0.35       0.1      0.1 
#> 2 bulk_b      0.1       0.05       0.5      0.35
```

## BayesPrism is a local optional backend

BayesPrism runs inside R and is useful when you want a Bayesian local
workflow. It is optional because it is much heavier than the dry-run
CIBERSORTx export.

``` r

bayesprism_result <- sn_deconvolve_bulk(
  x = pbmc,
  bulk = bulk,
  method = "bayesprism",
  cell_type_by = "cell_type",
  key = NULL,
  n_cores = 2,
  store_name = "pbmc3k_bayesprism"
)
```

The rule of thumb is: use dry-run mode to verify file contracts, then
run the backend deliberately in the environment where its dependencies
and credentials are configured.
