# Infer transcription-factor or pathway activity

`sn_run_regulatory_activity()` runs fast footprint-style activity
inference using DoRothEA regulons or PROGENy pathway models through
decoupleR. DoRothEA returns transcription-factor activities; PROGENy
returns pathway activities.

## Usage

``` r
sn_run_regulatory_activity(
  object,
  method = c("dorothea", "progeny"),
  assay = NULL,
  layer = "data",
  group_by = NULL,
  species = NULL,
  network = NULL,
  confidence_levels = c("A", "B", "C"),
  progeny_top = 500,
  minsize = 5L,
  store_name = "default",
  return_object = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object.

- method:

  One of `"dorothea"` or `"progeny"`.

- assay, layer:

  Assay and layer used to retrieve expression.

- group_by:

  Optional metadata column. When supplied, expression is averaged by
  group before activity inference.

- species:

  Species used to choose built-in resources.

- network:

  Optional user-supplied regulatory network. If omitted, DoRothEA uses
  the dorothea package data and PROGENy uses progeny or decoupleR
  resources when installed.

- confidence_levels:

  DoRothEA confidence levels to keep.

- progeny_top:

  Number of target genes per PROGENy pathway.

- minsize:

  Minimum target-set size passed to `decoupleR::run_ulm()`.

- store_name:

  Name used under `object@misc$regulatory_activity_results`.

- return_object:

  If `TRUE`, return the updated Seurat object.

- ...:

  Additional arguments passed to `decoupleR::run_ulm()`.

## Value

A Seurat object or stored-result list.
