.panel_size_pt <- function(plot) {
  g <- ggplot2::ggplotGrob(plot)
  panels <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  list(width = grid::convertWidth(g$widths[unique(panels$l)], "pt", valueOnly = TRUE),
       height = grid::convertHeight(g$heights[unique(panels$t)], "pt", valueOnly = TRUE))
}

test_that("every composition entry point renders requested panel dimensions", {
  d <- data.frame(sample = rep(paste0("s", 1:4), each = 4),
                  group = rep(c("case", "control"), each = 8),
                  cell = rep(c("B", "T"), 8), value = 1:16)
  plots <- list(
    sn_plot_bar(d, "group", "cell", panel_widths = 240, panel_heights = 180),
    sn_plot_sample_bar(d, "group", "cell", "sample", panel_widths = 240, panel_heights = 180),
    sn_plot_sample_boxplot(d, "group", "cell", "sample", panel_widths = 240, panel_heights = 180),
    sn_plot_histogram(d, "value", panel_widths = 240, panel_heights = 180)
  )
  for (style in c("A", "B", "C")) plots[[length(plots) + 1L]] <- sn_plot_sankey(
    d, flow_by = c("group", "cell"), style = style, panel_widths = 240, panel_heights = 180)
  for (p in plots) expect_equal(.panel_size_pt(p), list(width = 240, height = 180))
})

test_that("Sankey dimensions support derivation and facet sizes", {
  d <- data.frame(a = rep(c("X", "Y"), 2), b = rep(c("T", "B"), 2), cohort = rep(c("one", "two"), each = 2))
  override <- sn_plot_sankey(d, c("a", "b"), panel_widths = 240,
                             panel_heights = 180, aspect_ratio = 2)
  expect_equal(.panel_size_pt(override), list(width = 240, height = 180))
  p <- sn_plot_sankey(d, c("a", "b"), panel_widths = 240, aspect_ratio = 0.5)
  expect_equal(.panel_size_pt(p), list(width = 240, height = 120))
  p <- sn_plot_sankey(d, c("a", "b"), facet_col_by = "cohort", panel_widths = c(200, 300), panel_heights = 180)
  expect_equal(.panel_size_pt(p), list(width = c(200, 300), height = 180))
  for (size in list(0, -1, NA_real_, Inf, "large", numeric())) {
    expect_error(sn_plot_sankey(d, c("a", "b"), panel_widths = size), "positive finite")
  }
  expect_error(sn_plot_sankey(d, c("a", "b"), style = "C", panel_widths = c(200, 300)), "scalar")
})

test_that("C physical sizing keeps pie circles and derives a missing dimension", {
  d <- data.frame(a = "Young", b = c("T", "B"))
  p <- sn_plot_sankey(d, c("a", "b"), style = "C", panel_widths = 240, panel_heights = 180)
  built <- ggplot2::ggplot_build(p)
  panel <- built$layout$panel_params[[1]]
  pie <- p$layers[[4]]$data
  pie <- pie[pie$target == "T", ]
  width <- diff(range(pie$x)) / diff(panel$x.range) * 240
  height <- diff(range(pie$y)) / diff(panel$y.range) * 180
  expect_equal(width, height, tolerance = 0.002)
  p <- sn_plot_sankey(d, c("a", "b"), style = "C", panel_widths = 240)
  expect_equal(.panel_size_pt(p)$height, 240 * 0.45 * 1.83 / 1.06)
})
