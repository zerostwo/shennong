#' @export
sn_plot_boxplot <- function(data, x, y, sort = FALSE) {
  ggplot(data = data, mapping = aes(x = {{ x }}, y = {{ y }})) +
    geom_boxplot(
      outlier.shape = NA,
      staplewidth = 0.2,
      fatten = 1
    )
}

#' @export
sn_plot_barplot <- function(data, x, y, fill, sort_by = NULL) {
  ggplot(
    data = data,
    mapping = aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})
  ) +
    geom_col()
}
