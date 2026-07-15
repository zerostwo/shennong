# Run cell-cell communication inference

`sn_run_cell_communication()` wraps established communication backends
and stores a comparable ligand-receptor schema. Multiple backends can be
run together to calculate method concordance and a consensus rank.

## Usage

``` r
sn_run_cell_communication(
  object,
  method = c("liana", "cellchat", "cellphonedb", "nichenet", "multinichenet"),
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
  sample_by = NULL,
  consensus = TRUE,
  contrast = NULL,
  backend_control = list(),
  store_name = "default",
  return_object = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object.

- method:

  One or more of `"liana"`, `"cellchat"`, `"cellphonedb"`, `"nichenet"`,
  or `"multinichenet"`. The legacy alias `"nichenetr"` is accepted.

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
  `method = "nichenet"`.

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

- sample_by:

  Metadata column defining biological samples. When supplied, ligand and
  receptor expression is aggregated within each sample before condition
  comparison.

- consensus:

  If `TRUE` and multiple methods are requested, use the cross-method
  consensus ranking as the primary result table.

- contrast:

  Optional length-two condition contrast, with case first and reference
  second.

- backend_control:

  Named list of method-specific argument lists.

- store_name:

  Name used under `object@misc$cell_communication_results`.

- return_object:

  If `TRUE`, return the updated Seurat object; otherwise return the
  stored-result list.

- ...:

  Backend-specific arguments.

## Value

A Seurat object or stored-result list.
