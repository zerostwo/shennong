
# Shennong

<!-- badges: start -->

[![R-CMD-check](https://github.com/zerostwo/shennong/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zerostwo/shennong/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/zerostwo/shennong/graph/badge.svg)](https://codecov.io/gh/zerostwo/shennong)
[![lifecycle](https://img.shields.io/badge/lifecycle-Experimental-important.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

`Shennong` is an experimental R package for single-cell and broader
omics workflows. It provides Seurat-oriented helpers for example-data
loading, preprocessing, quality control, ambient RNA correction,
clustering and integration, and common visualization tasks.

## Installation

Install the development version from GitHub:

``` r
install.packages("remotes")
remotes::install_github("zerostwo/shennong")
```

Once the package is installed, you can check whether your local copy is
current and install the preferred release channel directly from R:

``` r
sn_check_version()
sn_install_shennong(channel = "github")
```

## What the package covers

- Example data loading with `sn_load_data()`
- Seurat object initialization and QC helpers
- Ambient contamination correction with SoupX or decontX
- Single-dataset clustering and Harmony-based integration through
  `sn_run_cluster()`
- Plot helpers for embeddings, violin plots, dot plots, boxplots, and
  bar plots
- Signature scoring, composition analysis, and enrichment helpers

## Quick start

The package now uses `sn_load_data()` as the main example-data entry
point.

``` r
library(Shennong)

pbmc <- sn_load_data("pbmc1k")
pbmc <- sn_initialize_seurat_object(pbmc, project = "pbmc1k", species = "human")
pbmc <- sn_filter_cells(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  plot = FALSE
)
pbmc <- sn_filter_genes(pbmc, min_cells = 3, plot = FALSE)
pbmc <- sn_run_cluster(
  pbmc,
  normalization_method = "seurat",
  resolution = 0.6
)

sn_plot_dim(pbmc, group_by = "seurat_clusters", label = TRUE)
```

## Ambient contamination correction

`sn_remove_ambient_contamination()` exposes a unified interface for
multiple methods.

``` r
filtered_path <- sn_load_data(
  dataset = "pbmc3k",
  matrix_type = "filtered",
  return_object = FALSE
)

raw_path <- sn_load_data(
  dataset = "pbmc3k",
  matrix_type = "raw",
  return_object = FALSE
)

ambient_counts <- sn_remove_ambient_contamination(
  x = filtered_path,
  raw = raw_path,
  method = "soupx",
  return_object = FALSE
)
```

For decontX, switch the method:

``` r
corrected_counts <- sn_remove_ambient_contamination(
  x = ambient_counts,
  method = "decontx",
  return_object = FALSE
)
```

If you write decontX output back into a Seurat object, the corrected
counts are stored in a separate layer by default and corrected per-cell
totals are added to metadata.

## Integration workflow

`sn_run_cluster()` is the main clustering entry point for both single
datasets and batch-aware workflows.

``` r
pbmc1k <- sn_load_data("pbmc1k")
pbmc3k <- sn_load_data("pbmc3k")

obj1 <- sn_initialize_seurat_object(pbmc1k, project = "pbmc1k")
obj1$sample <- "pbmc1k"
obj2 <- sn_initialize_seurat_object(pbmc3k, project = "pbmc3k")
obj2$sample <- "pbmc3k"

merged <- merge(obj1, y = obj2, add.cell.ids = c("pbmc1k", "pbmc3k"))

merged <- sn_run_cluster(
  object = merged,
  batch = "sample",
  normalization_method = "seurat",
  resolution = 0.6
)
```

## Development

Common local commands:

``` r
devtools::document()
testthat::test_local(stop_on_failure = TRUE)
pkgdown::build_site()
```

Package builds and checks can be run from the shell:

``` sh
R CMD build .
R CMD check --no-manual Shennong_*.tar.gz
```

Contributor guidance, testing conventions, and commit rules are
documented in `AGENTS.md` and `CONTRIBUTING.md`.

## Project status

The package is being modernized incrementally. Current work focuses on:

- stabilizing core APIs with test coverage
- aligning roxygen documentation with exports
- improving CI and pkgdown readiness
- keeping user-facing behavior compatible unless a change is explicitly
  documented

## Code of Conduct

This project follows the Contributor Covenant code of conduct:
<https://www.contributor-covenant.org/version/2/1/code_of_conduct/>.
