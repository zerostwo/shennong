# Add or refresh count-based QC percentages

Recalculate mitochondrial, ribosomal and hemoglobin percentages from a
selected count layer, including after ambient RNA correction.

## Usage

``` r
sn_add_qc_metrics(
  object,
  species = NULL,
  assay = NULL,
  layer = "counts",
  suffix = NULL
)
```

## Arguments

- object:

  A Seurat object with gene symbols as feature names.

- species:

  Either `"human"` or `"mouse"`. If NULL, use stored species metadata or
  infer species from the selected assay's feature names.

- assay:

  Assay to use; defaults to the object's default assay.

- layer:

  Count layer to use. An exact layer name takes precedence; otherwise
  split layers starting with `paste0(layer, ".")` are processed
  separately. Selected layers must not contain overlapping cells.

- suffix:

  String appended to `percent.mt`, `percent.ribo`, and `percent.hb`.
  NULL (default) automatically uses `"_corrected"` for
  `layer = "decontaminated_counts"` or its `decontaminated_counts.*`
  split layers, and `""` for other layers. An explicit string always
  takes precedence: `""` overwrites the original columns, while
  `"_corrected"` preserves them for comparison.

## Value

The Seurat object with three QC metadata columns and a command record.

## Details

Percentages are 100 times marker counts divided by total counts in the
selected layer, never by potentially stale `nCount_*` metadata.
Mitochondrial and ribosomal features use
[`sn_get_signatures()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_signatures.md);
hemoglobin matching retains the initialization patterns `^HB[^(P)]`
(human) and `^Hb[^(p)]` (mouse). No matching genes gives zero for
nonzero-total cells; zero-total cells give NaN. Cells absent from the
selected layers get NA. Supply non-negative counts, not normalized or
scaled expression. Counts, default assay, `nCount_*`, `nFeature_*`, and
stored species are unchanged. Sparse and BPCells matrices are not
converted to dense matrices.

## See also

[`sn_initialize_seurat_object()`](https://zerostwo.github.io/shennong/dev/reference/sn_initialize_seurat_object.md),
[`sn_assess_qc()`](https://zerostwo.github.io/shennong/dev/reference/sn_assess_qc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_add_qc_metrics(object, species = "human")
object <- sn_add_qc_metrics(object, layer = "decontaminated_counts")
head(object[[]][, c("percent.mt", "percent.mt_corrected")])
object <- sn_add_qc_metrics(object, assay = "RNA", layer = "counts.corrected",
                            suffix = "_corrected")
head(object[[]][, c("percent.mt", "percent.mt_corrected")])
} # }
```
