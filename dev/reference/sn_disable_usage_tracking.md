# Disable local workflow usage tracking

Restores the original Shennong function bindings and closes the logical
tracking session. There is no persistent database connection to close.

## Usage

``` r
sn_disable_usage_tracking()
```

## Value

Invisibly, `TRUE` when tracking was disabled and `FALSE` when it was
already disabled.
