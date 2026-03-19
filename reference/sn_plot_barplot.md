# Create a bar plot from a data frame

Create a bar plot from a data frame

## Usage

``` r
sn_plot_barplot(data, x, y, fill, sort_by = NULL)
```

## Arguments

- data:

  A data frame.

- x, y:

  Columns mapped to the x- and y-axes.

- fill:

  Column mapped to the fill aesthetic.

- sort_by:

  Currently reserved for future sorting support.

## Value

A ggplot object.

## Examples

``` r
plot_data <- data.frame(group = c("A", "B"), value = c(10, 15), type = c("x", "y"))
sn_plot_barplot(plot_data, x = group, y = value, fill = type)

```
