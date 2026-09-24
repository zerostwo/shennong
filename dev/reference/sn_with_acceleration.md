# Run an expression with selected acceleration patches active

Activates the requested ShennongOpt patches for the duration of `expr`
and restores the previous activation state afterwards. Patches without a
ShennongOpt counterpart run unaccelerated and are recorded as suppressed
in the acceleration provenance context.

## Usage

``` r
sn_with_acceleration(expr, name = NULL)
```

## Arguments

- expr:

  Expression to evaluate.

- name:

  Patch name(s) to activate.

## Value

The value of `expr`.

## Examples

``` r
if (FALSE) sn_with_acceleration(obj <- Seurat::RunPCA(obj), name = "seurat") # \dontrun{}
```
