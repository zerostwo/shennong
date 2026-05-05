# Data input, output, and project setup

Shennong is meant to make the first step of a single-cell project
boring: find the data, read it with one function, add explicit sample
metadata, and write reusable artifacts without changing APIs for every
file type.

This article uses the package PBMC3k example throughout. To keep package
checks and website builds fast, code chunks are shown by default and are
evaluated only when you render with `SHENNONG_RUN_VIGNETTES=true`.

## Start with PBMC3k

The package loader returns a Seurat object by default. If you need the
original 10x H5 path instead, set `return_object = FALSE`; this is
useful when you want to show exactly what was read from disk. You can
also request several example datasets at once; filtered matrices are
initialized separately and returned as one merged Seurat object with the
source dataset stored in the `sample` metadata column.

``` r

library(Shennong)
library(Seurat)
library(dplyr)

pbmc <- sn_load_data(dataset = "pbmc3k")
#> INFO [2026-05-05 20:32:36] Initializing Seurat object for project: pbmc3k.
#> INFO [2026-05-05 20:32:37] Running QC metrics for human.
#> INFO [2026-05-05 20:32:37] Seurat object initialization complete.

pbmc_h5 <- sn_load_data(
  dataset = "pbmc3k",
  return_object = FALSE
)

pbmc_h5
#> [1] "/home/runner/.shennong/data/pbmc3k_filtered_feature_bc_matrix.h5"

pbmc_merged <- sn_load_data(dataset = c("pbmc1k", "pbmc3k"))
#> INFO [2026-05-05 20:32:41] Initializing Seurat object for project: pbmc1k.
#> INFO [2026-05-05 20:32:41] Running QC metrics for human.
#> INFO [2026-05-05 20:32:42] Seurat object initialization complete.
#> INFO [2026-05-05 20:32:42] Initializing Seurat object for project: pbmc3k.
#> INFO [2026-05-05 20:32:42] Running QC metrics for human.
#> INFO [2026-05-05 20:32:42] Seurat object initialization complete.
table(pbmc_merged$sample)
#> 
#> pbmc1k pbmc3k 
#>   1215   2753
```

[`sn_load_data()`](https://songqi.org/shennong/dev/reference/sn_load_data.md)
uses the general Zenodo downloader underneath. Public Zenodo records do
not need a token; pass `token = ...` only for restricted or private
records that your account can access.

``` r

public_h5 <- sn_download_zenodo(
  record_id = "14884845",
  files = "pbmc3k_filtered_feature_bc_matrix.h5"
)

public_h5
#>                               pbmc3k_filtered_feature_bc_matrix.h5 
#> "/home/runner/.shennong/data/pbmc3k_filtered_feature_bc_matrix.h5"
```

## Read once, then initialize with metadata

[`sn_read()`](https://songqi.org/shennong/dev/reference/sn_read.md)
handles common tabular files, serialized objects, 10x matrices, H5/H5AD,
GMT files, and Shennong’s custom dispatchers. For a 10x H5 matrix, you
can read counts directly and then make the metadata decision explicit in
[`sn_initialize_seurat_object()`](https://songqi.org/shennong/dev/reference/sn_initialize_seurat_object.md).

``` r

counts <- sn_read(pbmc_h5)

pbmc <- sn_initialize_seurat_object(
  x = counts,
  project = "pbmc3k_demo",
  sample_name = "pbmc3k",
  study = "10x_pbmc",
  species = "human"
)
#> INFO [2026-05-05 20:32:44] Initializing Seurat object for project: pbmc3k_demo.
#> INFO [2026-05-05 20:32:44] Running QC metrics for human.
#> INFO [2026-05-05 20:32:44] Seurat object initialization complete.

pbmc
#> An object of class Seurat 
#> 54872 features across 2753 samples within 1 assay 
#> Active assay: RNA (54872 features, 0 variable features)
#>  1 layer present: counts
```

The important design choice is that provenance is not hidden in the
project name. `sample_name`, `study`, and `species` are separate inputs,
so downstream composition, integration, and interpretation steps can
reuse them.

``` r

head(pbmc[[]][, c("sample", "study", "nCount_RNA", "nFeature_RNA", "percent.mt")])
#>                  sample    study nCount_RNA nFeature_RNA percent.mt
#> AAACATACAACCAC-1 pbmc3k 10x_pbmc       2844         1046  2.7777778
#> AAACATTGAGCTAC-1 pbmc3k 10x_pbmc       5717         1722  3.2534546
#> AAACATTGATCAGC-1 pbmc3k 10x_pbmc       3766         1500  0.6903877
#> AAACCGTGCTTCCG-1 pbmc3k 10x_pbmc       2981         1153  1.5095606
#> AAACCGTGTATGCG-1 pbmc3k 10x_pbmc       1229          698  0.8950366
#> AAACGCACTGGTAC-1 pbmc3k 10x_pbmc       2536         1032  1.4984227
```

## Discover 10x folders before reading them

Real projects often start with several Cell Ranger outputs.
[`sn_list_10x_paths()`](https://songqi.org/shennong/dev/reference/sn_list_10x_paths.md)
scans a directory tree and returns the paths Shennong can feed back into
[`sn_initialize_seurat_object()`](https://songqi.org/shennong/dev/reference/sn_initialize_seurat_object.md).

``` r

tenx_paths <- sn_list_10x_paths(
  path = "data/cellranger",
  what = "filtered"
)

pbmc_samples <- sn_initialize_seurat_object(
  x = tenx_paths,
  species = "human",
  study = "pbmc_atlas"
)
```

This pattern keeps import code short while still making the biological
metadata explicit after discovery.

## Write durable intermediate files

Use
[`sn_write()`](https://songqi.org/shennong/dev/reference/sn_write.md)
and [`sn_read()`](https://songqi.org/shennong/dev/reference/sn_read.md)
together for analysis handoffs. The same API works for tables and
serialized R objects; Shennong chooses the writer from the file
extension.

``` r

outdir <- sn_set_path(file.path(tempdir(), "shennong-pbmc3k-io"))

metadata_path <- file.path(outdir, "pbmc3k_metadata.csv")
object_path <- file.path(outdir, "pbmc3k_initialized.qs2")

metadata_export <- cbind(cell = rownames(pbmc[[]]), pbmc[[]])
sn_write(metadata_export, metadata_path, row.names = FALSE)
sn_write(pbmc, object_path)

sn_check_file(c(metadata_path, object_path))

metadata <- sn_read(metadata_path, row_names = "cell")
pbmc_cached <- sn_read(object_path)

dim(metadata)
#> [1] 2753    8
pbmc_cached
#> An object of class Seurat 
#> 54872 features across 2753 samples within 1 assay 
#> Active assay: RNA (54872 features, 0 variable features)
#>  1 layer present: counts
```

`row_names = "cell"` is useful when a CSV carries cell IDs as a regular
column. This avoids ad hoc `read.csv(..., row.names = ...)` calls that
are easy to forget in later scripts.

## Prepare a reusable Zenodo upload

When an intermediate object or reference dataset should be reused by
another project, upload the exact files and a versioned manifest to
Zenodo with
[`sn_upload_zenodo()`](https://songqi.org/shennong/dev/reference/sn_upload_zenodo.md).
The default is draft-first: inspect the Zenodo draft, then rerun with
`publish = TRUE` only when the record is ready.

``` r

zenodo_plan <- sn_upload_zenodo(
  files = c(metadata_path, object_path),
  title = "PBMC3k initialized Shennong example",
  creators = "Duan, Songqi",
  version = "2026.05.04",
  sandbox = TRUE,
  dry_run = TRUE
)

zenodo_plan$files[, c("file", "size", "md5")]
#> # A tibble: 2 × 3
#>   file                      size md5                             
#>   <chr>                    <dbl> <chr>                           
#> 1 pbmc3k_metadata.csv     252144 30e88b334de6ed4d04a0393921f0f216
#> 2 pbmc3k_initialized.qs2 6229213 e9e72447e97d06d665224462532cf9a9
zenodo_plan$manifest_path
#> [1] "/tmp/RtmpUk5ndh/shennong_zenodo_manifest.json"
```

In a real upload, set `ZENODO_TOKEN` or `ZENODO_SANDBOX_TOKEN` and
remove `dry_run = TRUE`. The uploaded `shennong_zenodo_manifest.json`
records the dataset version, Shennong version, file sizes, and checksums
so future analyses can verify they are reusing the same data release.

## Add metadata exported from AnnData

When a Python workflow exports AnnData metadata or embeddings, use
[`sn_add_data_from_anndata()`](https://songqi.org/shennong/dev/reference/sn_add_data_from_anndata.md)
to merge those files back into the Seurat object with one explicit
import step.

``` r

pbmc <- sn_add_data_from_anndata(
  object = pbmc,
  metadata = "results/scanpy/obs.csv",
  reductions = c(umap = "results/scanpy/X_umap.csv")
)
```

This keeps cross-language handoffs visible in the analysis script
instead of leaving them as undocumented slot edits.

## Initialize a governed analysis project

For a new analysis repository,
[`sn_initialize_project()`](https://songqi.org/shennong/dev/reference/sn_initialize_project.md)
creates a lightweight project scaffold with Shennong/Codex-oriented
governance files. Use it when the analysis itself needs to be
reproducible, not just the R object.

``` r

sn_initialize_project(
  path = "~/projects/pbmc3k-shennong-demo",
  project_name = "pbmc3k-shennong-demo"
)
```

Before installing optional backends, inspect what Shennong considers
required, recommended, or optional for your workflow.

``` r

sn_check_version()

sn_install_shennong(source = "github", ref = "main")

deps <- sn_list_dependencies()
dplyr::count(deps, scope, source)

sn_install_dependencies(scope = "recommended")
```

For serialized project objects, prefer `.qs2` files in new workflows.
Legacy `.qs` files remain supported when the archived `qs` package is
already present in the active R library, but
[`sn_install_dependencies()`](https://songqi.org/shennong/dev/reference/sn_install_dependencies.md)
does not attempt to install `qs` on current R releases.

The package also ships Codex skill assets for Shennong-style analysis
projects. These helpers expose their installed locations without making
users search the package directory manually.

``` r

sn_get_codex_skill_path()
sn_install_codex_skill(destination = ".codex/skills")
```
