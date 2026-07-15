# CNV, Malignancy, and Metabolism

CNV and metabolic activity answer different questions, but both require
an explicit evidence boundary. CNV malignancy calls need trusted normal
cells; condition-level metabolism needs biological samples rather than
treating cells as replicates. Shennong stores both analyses through the
common result contract.

## Reference-aware CNV

Use cell names directly or a metadata column to declare normal
references. inferCNVpy is the default pixi backend. CopyKAT is also
available and is run sample by sample when `sample_by` is supplied.

``` r

tumor <- sn_run_cnv(
  tumor,
  method = "infercnvpy",
  reference_by = "cell_type",
  reference_cat = c("T cell", "Myeloid"),
  sample_by = "patient",
  store_name = "tumor_cnv"
)
```

Discover and retrieve the durable result without reading `object@misc`
directly. The primary table contains per-cell CNV, malignancy, call,
reference, and subclone fields. The other tables retain chromosome,
sample, and CNV-expression evidence.

``` r

sn_list_results(tumor, type = "cnv")
cnv <- sn_get_result(tumor, "cnv", "tumor_cnv")
head(cnv$tables$primary)
cnv$tables$sample_summary

sn_plot_cnv(tumor, "tumor_cnv", type = "heatmap")
sn_plot_cnv(tumor, "tumor_cnv", type = "umap")
sn_plot_cnv(tumor, "tumor_cnv", type = "association")
```

## Curated metabolic activity

The lightweight default uses curated pathways with UCell. GSVA, ssGSEA,
and a sparse-aware mean score are alternatives. Always supply
`sample_by` before a condition comparison; `group_by` optionally
stratifies by cell type or state.

``` r

sn_metabolic_signatures("human") |> names()

tumor <- sn_run_metabolism(
  tumor,
  method = "geneset",
  scoring_method = "ucell",
  sample_by = "patient",
  condition_by = "response",
  group_by = "cell_type",
  contrast = c("responder", "nonresponder"),
  store_name = "metabolism_response"
)

sn_list_results(tumor, type = "metabolism")
metabolism <- sn_get_result(tumor, "metabolism", "metabolism_response")
metabolism$tables$coverage
metabolism$tables$sample_scores
metabolism$tables$differential

sn_plot_metabolism(tumor, "metabolism_response", type = "heatmap")
sn_plot_metabolism(tumor, "metabolism_response", type = "differential")
```

scMetabolism is an optional R backend. scFEA and Compass remain
heavyweight external runtimes: pass a pathway-by-cell result or an
explicit runner through `backend_control`, after which Shennong applies
the same sample-level storage, comparison, and plotting contract.
