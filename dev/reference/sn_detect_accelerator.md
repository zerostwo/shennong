# Detect local accelerator support for pixi-managed Python methods

This helper performs lightweight command-line checks for GPUs.
CUDA/NVIDIA is currently the only accelerator profile Shennong uses to
select a pixi GPU environment automatically; other detected accelerators
are reported for the user's information and fall back to CPU unless
explicitly handled later.

## Usage

``` r
sn_detect_accelerator(quiet = FALSE)
```

## Arguments

- quiet:

  Logical; suppress status messages.

## Value

A named list with `has_gpu`, `backend`, `devices`, and detected CUDA
version when available.

## Examples

``` r
accel <- sn_detect_accelerator(quiet = TRUE)
accel$backend
#> [1] "cpu"
```
