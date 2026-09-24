# Data input, output, and project setup

Shennong begins after data distribution: it reads materialized files,
adds explicit sample metadata, and writes reusable analysis artifacts
without changing APIs for every file type. Dataset discovery, download,
caching, and publication are handled by the separate `ShennongData`
package.

This article uses a small, provenance-tracked public PBMC subset
materialized under `SHENNONG_REAL_DATA_DIR`. The data are used locally
for validation and pkgdown output but are intentionally not stored in
the Shennong repository.

## Start with materialized public data

The local fixture is a Seurat object prepared from public data. Its
acquisition manifest and checksums live alongside the ignored local data
cache, while the analysis code remains reviewable here.

``` r

library(Shennong)
library(Seurat)
library(dplyr)

pbmc <- qs2::qs_read(file.path(
  real_data_root,
  "single-cell", "kotliarov_pbmc.qs2"
))

fixture_summary <- data.frame(
  cells = ncol(pbmc),
  features = nrow(pbmc),
  assays = paste(names(pbmc@assays), collapse = ", "),
  biological_samples = length(unique(pbmc$real_sample)),
  acquisition_batches = length(unique(pbmc$real_batch))
)
knitr::kable(fixture_summary)
```

## Use ShennongData for discovery and materialization

Shennong does not wrap or re-export the data client. Use qualified
`ShennongData::` calls so the ownership boundary is explicit, inspect
the lazy query, and materialize only the assay/layer required by an
analysis.

``` r

connection <- ShennongData::sn_connect(
  url = "https://data.example.org",
  set_default = FALSE
)

remote <- ShennongData::sn_load_data(
  resource = "toil",
  connection = connection
)

rna <- ShennongData::sn_assay(remote, assay = "rna", layer = "expression")
ShennongData::sn_show_query(rna)
materialized <- ShennongData::collect(rna, allow_large = TRUE)
```

Pass `materialized` to
[`sn_initialize_seurat_object()`](https://zerostwo.github.io/shennong/dev/reference/sn_initialize_seurat_object.md)
or another Shennong analysis entry point. Data publication and cache
policy remain in `ShennongData`, not the analysis package.

## Read once, then initialize with metadata

[`sn_read()`](https://zerostwo.github.io/shennong/dev/reference/sn_read.md)
handles common tabular files, serialized objects, 10x matrices, H5/H5AD,
GMT files, and Shennong’s custom dispatchers. For a 10x H5 matrix, you
can read counts directly and then make the metadata decision explicit in
[`sn_initialize_seurat_object()`](https://zerostwo.github.io/shennong/dev/reference/sn_initialize_seurat_object.md).

``` r

counts <- SeuratObject::LayerData(pbmc, assay = "RNA", layer = "counts")
metadata <- pbmc[[]]

pbmc <- sn_initialize_seurat_object(
  x = counts,
  metadata = metadata,
  project = "kotliarov_pbmc",
  species = "human"
)

pbmc
```

The important design choice is that provenance is not hidden in the
project name. The materialized metadata retain the real sample and batch
fields, while species is an explicit initialization input.

``` r

head(pbmc[[]][, c("real_sample", "real_batch", "nCount_RNA", "nFeature_RNA", "percent.mt")])
```

## Discover 10x folders before reading them

Real projects often start with several Cell Ranger outputs.
[`sn_list_10x_paths()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_10x_paths.md)
scans a directory tree and returns the paths Shennong can feed back into
[`sn_initialize_seurat_object()`](https://zerostwo.github.io/shennong/dev/reference/sn_initialize_seurat_object.md).

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
[`sn_write()`](https://zerostwo.github.io/shennong/dev/reference/sn_write.md)
and
[`sn_read()`](https://zerostwo.github.io/shennong/dev/reference/sn_read.md)
together for analysis handoffs. The same API works for tables and
serialized R objects; Shennong chooses the writer from the file
extension.

``` r

outdir <- sn_set_path(file.path(tempdir(), "shennong-kotliarov-io"))

metadata_path <- file.path(outdir, "kotliarov_metadata.csv")
object_path <- file.path(outdir, "kotliarov_initialized.qs2")

metadata_export <- cbind(cell = rownames(pbmc[[]]), pbmc[[]])
sn_write(metadata_export, metadata_path, row.names = FALSE)
sn_write(pbmc, object_path)

sn_check_file(c(metadata_path, object_path))

metadata <- sn_read(metadata_path, row_names = "cell")
pbmc_cached <- sn_read(object_path)

dim(metadata)
pbmc_cached
```

`row_names = "cell"` is useful when a CSV carries cell IDs as a regular
column. This avoids ad hoc `read.csv(..., row.names = ...)` calls that
are easy to forget in later scripts.

For large Seurat objects, switch bulky assay layers to the BPCells
backend before serializing the Seurat object. The returned object keeps
the same Seurat assay/layer interface, but the selected matrices live in
BPCells directories.

``` r

pbmc <- sn_set_layer_backend(
  pbmc,
  backend = "bpcells",
  directory = file.path(outdir, "pbmc3k_bpcells"),
  layers = c("counts", "data"),
  overwrite = TRUE
)

sn_write(pbmc, file.path(outdir, "pbmc3k_bpcells_bound.qs2"))
```

When an operation explicitly requires an in-memory sparse matrix, switch
only the required layer back. This materializes the complete selected
layer as a `dgCMatrix`; it does not delete the external BPCells matrix
directory.

``` r

pbmc <- sn_set_layer_backend(
  pbmc,
  backend = "memory",
  assays = "RNA",
  layers = "decontaminated_counts"
)
```

[`sn_convert_bpcells()`](https://zerostwo.github.io/shennong/dev/reference/sn_convert_bpcells.md)
remains a compatibility wrapper for
`sn_set_layer_backend(backend = "bpcells")`.

[`sn_initialize_seurat_object()`](https://zerostwo.github.io/shennong/dev/reference/sn_initialize_seurat_object.md)
also accepts a BPCells `IterableMatrix` directly and preserves that
on-disk backend instead of converting the complete matrix to an
in-memory sparse matrix.

``` r

counts <- BPCells::open_matrix_dir("data/processed/pbmc_counts.bpcells")
pbmc <- sn_initialize_seurat_object(
  x = counts,
  metadata = cell_metadata,
  project = "PBMC"
)
inherits(SeuratObject::LayerData(pbmc, layer = "counts"), "IterableMatrix")
```

## Publish reusable data outside the analysis package

When an intermediate object or reference dataset should be shared, hand
the exact files and their checksums to the `ShennongData` publication
workflow. Shennong intentionally stops at producing validated analysis
artifacts; it no longer owns repository discovery, download, caching,
credentials, or remote publication state.

## Add metadata exported from AnnData

When a Python workflow exports AnnData metadata or embeddings, use
[`sn_add_data_from_anndata()`](https://zerostwo.github.io/shennong/dev/reference/sn_add_data_from_anndata.md)
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
[`sn_initialize_project()`](https://zerostwo.github.io/shennong/dev/reference/sn_initialize_project.md)
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

sn_install_shennong(channel = "github", ref = "main")

# From a local Shennong source checkout:
sn_install_shennong(channel = "local", source = ".")

deps <- sn_list_dependencies()
dplyr::count(deps, scope, source)

sn_install_dependencies(scope = "recommended")
```

Serialized project objects use `.qs2` and the `qs2` package:

``` r

sn_write(object, "data/processed/object.qs2")
object <- sn_read("data/processed/object.qs2")
```

Legacy `.qs` import/export and its automatic installer have been
removed.
[`sn_read()`](https://zerostwo.github.io/shennong/dev/reference/sn_read.md)
and
[`sn_write()`](https://zerostwo.github.io/shennong/dev/reference/sn_write.md)
reject that format before IO or installation, even if an older
serializer is installed. Convert existing files in an older compatible
environment before using them here; changing the extension alone does
not convert the serialization format.

The package also ships Codex skill assets for Shennong-style analysis
projects. These helpers expose their installed locations without making
users search the package directory manually.

``` r

sn_get_codex_skill_path()
sn_install_codex_skill(destination = ".codex/skills")
```

## Preserve matrix objects during handoff

[`sn_write()`](https://zerostwo.github.io/shennong/dev/reference/sn_write.md)
preserves matrix classes and dimnames in RDS, RData, and QS2. CSV and
other tabular formats continue to write tables. Dense numeric and sparse
matrices can be written directly to BPCells directories or 10x HDF5:

``` r

counts <- SeuratObject::LayerData(pbmc, assay = "RNA", layer = "counts")
sn_write(counts, "counts.rds", auto_install = FALSE)
stopifnot(identical(readRDS("counts.rds"), counts))
sn_write(counts, "counts.h5", auto_install = FALSE)
sn_write(counts, "counts.bpcells", auto_install = FALSE)
```
