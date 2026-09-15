# Plotting stored results and interpretation summaries

Use only the stages needed by the requested analysis. Each section is an API
reference, not a requirement to run all listed methods.

## 8

Prefer `sn_plot_result()` for stored/downstream analysis results; inspect
   valid analysis/view pairs with `sn_list_plot_methods()`. Use
   `sn_plot_distribution()` for numeric violin/box/histogram/density/ridge
   views and `sn_plot_association()` for cell- or sample-level correlations.
   Use `sn_plot_dim()`, `sn_plot_feature()`, `sn_plot_heatmap()`,
   `sn_plot_violin()`, `sn_plot_dot()`, `sn_plot_boxplot()`,
   `sn_plot_barplot()`, `sn_plot_composition()`, and `sn_plot_milo()` for
   package-style plots; resolve reusable colors with `sn_list_palettes()` and
   `sn_get_palette()`.
   `sn_plot_composition()` accepts a Seurat object directly. Use ordinary
   `type = "bar"` for descriptive cell counts/proportions; use
   `type = "sample_bar"` or `"sample_boxplot"` with an explicit `sample_by`
   when biological samples are the replicates; use `type = "alluvial"` with
   two or more `flow_by` metadata columns for source/annotation flows. Use
   `unit_by` for donor metadata repeated across cells so donor distributions do
   not count cells as independent donors.

## 9

Build prompts or stored-result summaries with the interpretation helpers
   when a narrative or report-ready output is needed.
