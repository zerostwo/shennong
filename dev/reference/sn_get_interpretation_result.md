# Retrieve a stored interpretation result from a Seurat object

Retrieve a stored interpretation result from a Seurat object

## Usage

``` r
sn_get_interpretation_result(object, result_id = NULL)
```

## Arguments

- object:

  A `Seurat` object.

- result_id:

  Identifier of the stored interpretation result.

## Value

The stored interpretation-result list.

## Examples

``` r
if (FALSE) { # \dontrun{
interpretation <- sn_get_interpretation_result(seurat_obj, "annotation_note")
} # }
```
