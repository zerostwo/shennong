---
name: manage-shennong-results
description: Use when storing, listing, retrieving, and reusing Shennong results on Seurat objects, especially for differential expression, enrichment, and interpretation outputs.
---

# manage-shennong-results

## Purpose

Teach the agent how to store, retrieve, and reuse results through the Shennong
package API.

## When To Use

- storing DE or enrichment outputs
- retrieving stored results for downstream plots or interpretation
- auditing what is available on a Seurat object

## Required Inputs

- Seurat object
- stored-result name
- intended reuse pattern

## Required Outputs

- stored or retrieved result
- explicit reuse path

## Rules

- use package result stores rather than unstructured `misc` access in user code
- prefer `sn_list_results()` for discovery
- prefer `sn_get_*_result()` helpers for retrieval
- remember that DE, enrichment, milo, deconvolution, and interpretation each
  have separate storage entry points
- when stored outputs become durable project assets, promote them according to
  the initialized-project governance rather than leaving them only in transient state

## References

- `../_shared/references/package_api_map.md`
- `../_shared/references/workflow_recipes.md`

## Procedure

1. Store outputs through Shennong APIs that support object return.
2. Discover result names with `sn_list_results()`.
3. Retrieve specific results with `sn_get_de_result()`,
   `sn_get_enrichment_result()`, `sn_get_milo_result()`,
   `sn_get_deconvolution_result()`,
   `sn_get_cell_communication_result()`,
   `sn_get_regulatory_activity_result()`, or
   `sn_get_interpretation_result()`.
4. Reuse stored outputs in visualization or interpretation helpers.

## Common Mistakes

- hard-coding untracked `object@misc` paths in analysis scripts
- storing outputs without clear names
- failing to preserve reusable provenance

## Examples

- `sn_list_results(object)`
- `sn_get_de_result(object, de_name = "cluster_markers")`
- `sn_get_enrichment_result(object, enrichment_name = "cluster_pathways")`
- `sn_get_milo_result(object, milo_name = "condition_da")`
- `sn_get_deconvolution_result(object, deconvolution_name = "bulk_mix")`
- `sn_get_cell_communication_result(object, communication_name = "cellchat")`
- `sn_get_regulatory_activity_result(object, activity_name = "dorothea")`
- `sn_get_interpretation_result(object, interpretation_name = "annotation_note")`
