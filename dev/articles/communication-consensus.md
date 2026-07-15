# Comparable cell-cell communication and consensus

Shennong keeps backend-native artifacts while mapping communication
evidence to one table contract. The common columns include sender,
receiver, ligand, receptor, score, p/q values, rank, method, condition,
sample, pathway, target-gene evidence, evidence source, and a reserved
spatial-distance field.

## Run one or several backends

Use a vector of methods when comparable backends should contribute to a
consensus. Method-specific arguments belong in `backend_control`.

``` r

object <- sn_run_cell_communication(
  object,
  method = c("liana", "cellchat"),
  group_by = "cell_type",
  sample_by = "patient",
  condition_by = "response",
  contrast = c("Responder", "Nonresponder"),
  consensus = TRUE,
  backend_control = list(
    liana = list(resource = "consensus"),
    cellchat = list(min_cells = 20)
  ),
  store_name = "tumor_communication"
)
```

NicheNet adds explicit ligand-target links when its prior network and
ligand-target matrix are supplied. MultiNicheNet requires biological
sample and condition columns and delegates the pseudobulk differential
analysis to the official backend.

``` r

object <- sn_run_cell_communication(
  object,
  method = "multinichenet",
  group_by = "cell_type",
  sample_by = "patient",
  condition_by = "response",
  contrast = c("Responder", "Nonresponder"),
  sender = c("Myeloid", "Fibroblast"),
  receiver = "Tumor",
  lr_network = lr_network,
  ligand_target_matrix = ligand_target_matrix,
  store_name = "differential_communication"
)
```

## Discover and retrieve stored evidence

``` r

sn_list_results(object, type = "cell_communication")

communication <- sn_get_result(
  object,
  type = "cell_communication",
  name = "tumor_communication"
)

communication$tables$primary
communication$tables$backend_raw
communication$tables$consensus
communication$tables$method_concordance
communication$tables$sample_evidence
communication$tables$condition_comparison
communication$tables$ligand_targets
communication$warnings
```

The sample-evidence table recomputes ligand and receptor expression
within each biological sample. Condition effects therefore use patients
or samples as replicates. They do not treat cells as independent
experimental units.

## Plot the stored result

``` r

sn_plot_communication(object, "tumor_communication", type = "bubble")
sn_plot_communication(object, "tumor_communication", type = "heatmap")
sn_plot_communication(object, "tumor_communication", type = "network")
sn_plot_ligand_target(object, "differential_communication")
sn_plot_communication_comparison(object, "tumor_communication")
```
