# Cell compostition: Barplot, boxplot
#' @export
sn_plot_barplot <- function(x, group_by, variable) {
  metadata <- x@meta.data
  adata <- table(metadata[[group_by]], metadata[[variable]]) |>
    prop.table(margin = 1) |>
    as.data.frame() |>
    mutate(Freq = Freq * 100)
  colnames(x = adata) <- c(group_by, variable, "percentage")

  adata |>
    ggplot(mapping = aes(x = .data[[group_by]], y = percentage, fill = .data[[variable]])) +
    geom_col()
}
