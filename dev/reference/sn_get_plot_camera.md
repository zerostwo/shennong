# Retrieve a reproducible 3D embedding camera

Read the JSON downloaded by an interactive embedding, a camera list, or
the initial camera attached to a static plot/widget. A standalone
browser cannot mutate an R session: use its Download camera or Copy R
camera button after adjusting the view, then pass the returned value to
`camera` in either plotter.

## Usage

``` r
sn_get_plot_camera(x = NULL)
```

## Arguments

- x:

  A JSON file, JSON string, camera list, or styled embedding
  plot/widget.

## Value

A list with azimuth, elevation and roll in degrees, positive zoom, and
two-element pan in normalized screen coordinates. Angles use an
orthographic projection shared by both renderers.

## Examples

``` r
camera <- sn_get_plot_camera(list(azimuth = 45, elevation = 20))
camera
```
