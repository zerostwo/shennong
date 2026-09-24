# Plot a stored or direct Shennong analysis result

`sn_plot_result()` is the canonical result-aware visualization entry
point. It uses one stable interface across analysis domains: `object`,
`analysis_type`, `result_id`, and `view`. Existing domain-specific
`sn_plot_*()` functions remain available as compatibility and advanced
interfaces.

## Usage

``` r
sn_plot_result(
  object,
  analysis_type = NULL,
  result_id = NULL,
  view = NULL,
  ...
)
```

## Arguments

- object:

  A Seurat object with stored results, a direct Shennong result, or a
  compatible table/legacy assessment supported by the chosen analysis
  type.

- analysis_type:

  Optional analysis type. Required for ambiguous Seurat stores and for
  table/legacy inputs. See
  [`sn_list_plot_methods()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_plot_methods.md).

- result_id:

  Optional stored result identifier. A unique result or one named
  `"default"` is selected automatically when omitted.

- view:

  Optional registered plot view. The analysis-specific default is used
  when omitted.

- ...:

  View-specific parameters forwarded to the specialized plotter.

## Value

Usually a `ggplot` object; heatmap or multi-panel backends may return a
compatible composed plot object.

## Details

A direct Shennong result supplies its own `analysis_type`. For Seurat
input, the function retrieves a stored result. If `analysis_type` or
`result_id` is omitted, it is selected only when unambiguous; a stored
result named `"default"` is preferred within a selected analysis type.

## Examples

``` r
de_table <- data.frame(
  gene = c("G1", "G2", "G3"),
  log2_fold_change = c(2, -1.5, 0.2),
  adjusted_p_value = c(0.01, 0.03, 0.8)
)
sn_plot_result(de_table, analysis_type = "de", view = "volcano")
if (FALSE) { # \dontrun{
sn_list_results(seurat_object, type = "de")
sn_plot_result(
  seurat_object,
  analysis_type = "de",
  result_id = "cluster_markers",
  view = "volcano"
)
} # }
```
