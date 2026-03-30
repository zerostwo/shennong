# Bulk Deconvolution Workflow

This workflow shows how to combine a single-cell reference with bulk
RNA-seq mixtures for cell-composition inference.

Shennong currently supports two backends:

- `BayesPrism` for local Bayesian deconvolution
- `CIBERSORTx` for local container-based deconvolution

## 1. Build a single-cell reference

``` r
data("pbmc_small", package = "Shennong")
ref <- pbmc_small
ref$cell_type <- rep(c("T_like", "B_like", "Myeloid_like"), length.out = ncol(ref))
```

## 2. Create a mock bulk matrix

``` r
counts <- SeuratObject::LayerData(ref, assay = "RNA", layer = "counts")
bulk <- cbind(
  sample_a = rowSums(counts[, 1:40, drop = FALSE]),
  sample_b = rowSums(counts[, 41:ncol(counts), drop = FALSE])
)
bulk[1:5, 1:2]
#>                 sample_a sample_b
#> DDX11L16               0        0
#> WASH7P                 1        4
#> MIR1302-2HG            0        0
#> FAM138A                0        0
#> ENSG00000308361        0        0
```

## 3. Prepare a local CIBERSORTx run

``` r
bundle <- sn_deconvolve_bulk(
  ref,
  bulk = bulk,
  method = "cibersortx",
  cell_type_col = "cell_type",
  layer = "counts",
  outdir = tempdir(),
  prefix = "pbmc_demo",
  cibersortx_email = "user@example.org",
  cibersortx_token = "replace-with-real-token",
  cibersortx_dry_run = TRUE,
  return_object = FALSE
)

bundle$files
#> $single_cell_reference
#> [1] "/tmp/RtmpP2i3JS/sample_file_for_cibersort.txt"
#> 
#> $mixture
#> [1] "/tmp/RtmpP2i3JS/mixture_file_for_cibersort.txt"
#> 
#> $signature_matrix
#> [1] "/tmp/RtmpP2i3JS/CIBERSORTx_sample_file_for_cibersort_inferred_phenoclasses.CIBERSORTx_sample_file_for_cibersort_inferred_refsample.bm.K999.txt"
#> 
#> $result
#> [1] "/tmp/RtmpP2i3JS/CIBERSORTx_pbmc_demo_Results.txt"
bundle$artifacts$commands
#> $create_signature
#> [1] "docker run -v '/tmp/RtmpP2i3JS':/src/data:z -v '/tmp/RtmpP2i3JS':/src/outdir:z cibersortx/fractions --single_cell TRUE --username 'user@example.org' --token 'replace-with-real-token' --refsample 'sample_file_for_cibersort.txt' --G.min 300 --G.max 500 --q.value 0.01 --filter FALSE --k.max 999 --remake FALSE --replicates 5 --sampling 0.5 --fraction 0.75"
#> 
#> $deconvolve
#> [1] "docker run -v '/tmp/RtmpP2i3JS':/src/data:z -v '/tmp/RtmpP2i3JS':/src/outdir:z cibersortx/fractions --single_cell TRUE --username 'user@example.org' --token 'replace-with-real-token' --mixture 'mixture_file_for_cibersort.txt' --sigmatrix 'CIBERSORTx_sample_file_for_cibersort_inferred_phenoclasses.CIBERSORTx_sample_file_for_cibersort_inferred_refsample.bm.K999.txt' --perm 0 --label 'pbmc_demo' --rmbatchBmode FALSE --rmbatchSmode FALSE --sourceGEPs 'CIBERSORTx_sample_file_for_cibersort_inferred_phenoclasses.CIBERSORTx_sample_file_for_cibersort_inferred_refsample.bm.K999.txt' --QN FALSE --absolute FALSE --abs_method 'sig.score'"
```

In `dry_run` mode, Shennong writes the local input files and returns the
exact container commands. To run CIBERSORTx locally, provide real
credentials and set `cibersortx_dry_run = FALSE`.

``` r
ref <- sn_deconvolve_bulk(
  ref,
  bulk = bulk,
  method = "cibersortx",
  cell_type_col = "cell_type",
  layer = "counts",
  cibersortx_email = Sys.getenv("SHENNONG_CIBERSORTX_EMAIL"),
  cibersortx_token = Sys.getenv("SHENNONG_CIBERSORTX_TOKEN"),
  store_name = "pbmc_cibersortx"
)
```

If you already have a local `CIBERSORTx_*_Results.txt` file, you can
still import it directly through `cibersortx_result = ...`.

## 4. Run BayesPrism locally

BayesPrism is an optional GitHub dependency:

``` r
remotes::install_github("Danko-Lab/BayesPrism/BayesPrism")
```

Then run:

``` r
ref <- sn_deconvolve_bulk(
  ref,
  bulk = bulk,
  method = "bayesprism",
  cell_type_col = "cell_type",
  layer = "counts",
  store_name = "pbmc_bayesprism"
)
```

## 5. Retrieve stored deconvolution results

``` r
sn_list_results(ref)
sn_get_deconvolution_result(ref, "pbmc_cibersortx")
```
