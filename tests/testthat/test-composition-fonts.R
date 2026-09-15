.composition_text_sizes <- function(g) {
  own <- if (inherits(g, "text")) g$gp$fontsize else numeric()
  children <- c(g$grobs, as.list(g$children))
  c(own, unlist(lapply(children, .composition_text_sizes), use.names = FALSE))
}

test_that("all rendered composition text defaults to eight points", {
  d <- data.frame(sample = rep(paste0("s", 1:4), each = 4),
                  group = rep(c("case", "control"), each = 8),
                  cell = rep(c("B", "T"), 8), value = 1:16)
  plots <- list(
    sn_plot_bar(d, "group", "cell", title = "Title", facet_col_by = "group"),
    sn_plot_sample_bar(d, "group", "cell", "sample", title = "Title"),
    sn_plot_sample_boxplot(d, "group", "cell", "sample", title = "Title"),
    sn_plot_histogram(d, "value", "group", title = "Title")
  )
  for (style in c("A", "B", "C")) plots[[length(plots) + 1L]] <- sn_plot_sankey(
    d, c("group", "cell"), style = style, title = "Title", facet_col_by = "group")
  for (plot in plots) {
    sizes <- .composition_text_sizes(ggplot2::ggplotGrob(plot))
    expect_gt(length(sizes), 0)
    expect_equal(sizes, rep(8, length(sizes)), tolerance = 1e-8)
  }
})

test_that("C pie toggle only removes pie layers and preserves flow weights", {
  d <- data.frame(group = c("Young", "Aged", "Young"), cell = c("T", "T", "B"))
  yes <- sn_plot_sankey(d, c("group", "cell"), style = "C", show_pies = TRUE)
  no <- sn_plot_sankey(d, c("group", "cell"), style = "C", show_pies = FALSE)
  expect_true("show_pies" %in% names(formals(sn_plot_sankey)))
  pie_layer <- function(p) any(vapply(p$layers, function(l) "pie" %in% names(l$data), logical(1)))
  expect_true(pie_layer(yes))
  expect_false(pie_layer(no))
  expect_equal(yes$layers[[1]]$data, no$layers[[1]]$data)
  expect_no_error(ggplot2::ggplotGrob(no))
})
