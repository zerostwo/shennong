# Activate ShennongOpt acceleration patches

Rebinds accelerated implementations into their upstream namespaces.
Activation is process-local and idempotent; pass patch names to activate
a subset, or omit names to activate eligible installed patches. On
Seurat 5.6 and later, Seurat handles variable feature selection, PCA,
neighbor graphs and clustering. Requests for those four companion
patches return FALSE. The single-covariate ScaleData regression patch
remains eligible.

## Usage

``` r
sn_enable_acceleration(name = NULL)
```

## Arguments

- name:

  Optional patch name(s) to activate.

## Value

Invisible named logical vector of activation results.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_enable_acceleration()
sn_disable_acceleration()
} # }
```
