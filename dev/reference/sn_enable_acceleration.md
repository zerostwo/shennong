# Activate ShennongOpt acceleration patches

Rebinds accelerated implementations into their upstream namespaces.
Activation is process-local and idempotent; pass patch names to activate
a subset, or nothing to activate every registered patch whose upstream
package is installed. See `ShennongOpt::sn_list_accelerations()` for
available names.

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
