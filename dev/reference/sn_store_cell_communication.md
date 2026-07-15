# Store a cell-cell communication result on a Seurat object

Store a cell-cell communication result on a Seurat object

## Usage

``` r
sn_store_cell_communication(
  object,
  result,
  store_name = "default",
  method = "cellchat",
  backend = method,
  group_by = NULL,
  sender = NULL,
  receiver = NULL,
  species = NULL,
  artifacts = NULL,
  raw_result = NULL,
  consensus_result = NULL,
  sample_evidence = NULL,
  comparison = NULL,
  concordance = NULL,
  ligand_targets = NULL,
  warnings = character(),
  sample_by = NULL,
  condition_by = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- result:

  Communication result table.

- store_name:

  Name used under `object@misc$cell_communication_results`.

- method:

  User-facing communication method label.

- backend:

  Canonical backend identifier. This differs from `method` only when
  preserving a legacy alias such as `"nichenetr"`.

- group_by:

  Metadata column used for groups.

- sender, receiver:

  Optional sender/receiver labels.

- species:

  Optional species label.

- artifacts:

  Optional backend-specific artifacts.

- raw_result, consensus_result, sample_evidence, comparison,
  concordance, ligand_targets:

  Optional standardized secondary tables stored in the unified result.

- warnings:

  Character vector of backend warnings retained for audit.

- sample_by, condition_by:

  Optional sample and condition metadata columns.

- return_object:

  If `TRUE`, return the updated object.

## Value

A Seurat object or stored-result list.
