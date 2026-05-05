# Run cell-cell communication inference

`sn_run_cell_communication()` wraps established communication backends:
CellChat for global interaction networks, NicheNet for
sender-to-receiver ligand activity, and LIANA for consensus
ligand-receptor scoring when the optional package is installed.

## Usage

``` r
sn_run_cell_communication(
  object,
  method = c("cellchat", "nichenetr", "liana"),
  group_by,
  assay = NULL,
  layer = "data",
  species = NULL,
  sender = NULL,
  receiver = NULL,
  geneset = NULL,
  background_genes = NULL,
  condition_by = NULL,
  condition_oi = NULL,
  condition_reference = NULL,
  ligand_target_matrix = NULL,
  lr_network = NULL,
  expressed_pct = 0.1,
  top_n = 50,
  cellchat_db = NULL,
  min_cells = 10,
  population_size = FALSE,
  raw_use = TRUE,
  resource = NULL,
  store_name = "default",
  return_object = TRUE,
  condition_col = NULL,
  ...
)
```

## Arguments

- object:

  A Seurat object.

- method:

  One of `"cellchat"`, `"nichenetr"`, or `"liana"`.

- group_by:

  Metadata column defining sender/receiver cell groups.

- assay, layer:

  Assay and layer used to retrieve expression.

- species:

  Species passed to species-aware resources. Defaults to
  `sn_get_species(object)` when available.

- sender, receiver:

  Sender and receiver group labels required by `method = "nichenetr"`.

- geneset:

  Receiver gene set of interest for NicheNet. If omitted, receiver
  differential expression is computed from `condition_by`,
  `condition_oi`, and `condition_reference`.

- background_genes:

  Background expressed genes for NicheNet.

- condition_by:

  Metadata column containing receiver conditions.

- condition_oi, condition_reference:

  Receiver condition contrast used to derive `geneset` for NicheNet.

- ligand_target_matrix, lr_network:

  NicheNet prior matrices/networks. These must be supplied for
  `method = "nichenetr"`.

- expressed_pct:

  Minimum fraction of cells expressing a gene before it is considered
  expressed in NicheNet filtering.

- top_n:

  Number of top NicheNet ligands to keep.

- cellchat_db:

  Optional CellChat database object.

- min_cells:

  Minimum cells per group passed to CellChat filtering.

- population_size, raw_use:

  CellChat communication-probability controls.

- resource:

  Optional LIANA resource.

- store_name:

  Name used under `object@misc$cell_communication_results`.

- return_object:

  If `TRUE`, return the updated Seurat object; otherwise return the
  stored-result list.

- condition_col:

  Deprecated alias for `condition_by`.

- ...:

  Backend-specific arguments.

## Value

A Seurat object or stored-result list.
