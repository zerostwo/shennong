test_that("composition string selectors work with variables and preserve old calls", {
  d <- data.frame(group = c("b", "a", "b"), cell = c("T", "B", "B"))
  column <- "cell"
  p <- sn_plot_bar(d, x_by = "group", fill_by = column)
  old <- sn_plot_composition(d, x = group, fill = cell)
  expect_equal(p$data, old$data)
  expect_no_error(ggplot2::ggplot_build(p))
  expect_error(sn_plot_composition(d, x = group, x_by = "group", fill_by = column), "Supply only")
  expect_error(sn_plot_bar(d, x_by = c("group", "cell"), fill_by = column), "single non-empty")
  expect_error(sn_plot_bar(d, x_by = "group", fill_by = column, fill = "cell"), "Dedicated")
})

test_that("all composition dedicated entrypoints agree with the dispatcher", {
  d <- data.frame(sample = rep(paste0("s", 1:4), each = 4),
                  group = rep(c("case", "control"), each = 8),
                  cell = rep(c("B", "T"), 8), age = rep(c(20, 30, 40, 50), each = 4))
  for (type in c("sample_bar", "sample_boxplot")) {
    dedicated <- get(paste0("sn_plot_", type))(d, x_by = "group", fill_by = "cell", sample_by = "sample")
    unified <- sn_plot_composition(d, type = type, x_by = "group", fill_by = "cell", sample_by = "sample")
    expect_equal(dedicated$data, unified$data)
    expect_no_error(ggplot2::ggplot_build(dedicated))
  }
  p <- sn_plot_histogram(d, x_by = "age", unit_by = "sample")
  expect_equal(nrow(p$data), 4L)
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("Sankey B preserves source/target order and uses destination color", {
  d <- data.frame(a = factor(c("X", "Y", "X"), c("Y", "X")),
                  b = factor(c("X", "Y", "Z"), c("Z", "Y", "X")), count = c(8, 1, 1))
  p <- sn_plot_sankey(d, flow_by = c("a", "b"), style = "B")
  nodes <- p$layers[[2]]$data
  expect_identical(nodes$label[nodes$axis == 1], c("Y", "X"))
  expect_identical(nodes$label[nodes$axis == 2], c("Z", "Y", "X"))
  expect_equal(sum(vapply(p$layers, function(l) inherits(l$geom, "GeomSegment"), logical(1))), 3)
  expect_equal(levels(p$layers[[1]]$data$.sn_fill), levels(d$b))
  expect_no_error(ggplot2::ggplot_build(p))
  expect_error(sn_plot_sankey(d, flow_by = c("a", "b", "count"), style = "B"), "exactly two")
})

test_that("Sankey C pies use retained weights within each target and facet", {
  d <- data.frame(source = factor(c("Young", "Aged", "Young", "Aged", "Young", "Aged"),
                                  c("Young", "Aged")),
                  target = c("T", "T", "B", "B", "T", "T"),
                  tissue = c(rep("one", 4), "two", "two"), weight = c(3, 1, 1, 3, 1, 4))
  p <- sn_plot_sankey(d, flow_by = c("source", "target"), style = "C", y_by = "weight",
                      facet_row_by = "tissue", palette = c(Young = "grey50", Aged = "steelblue"))
  slices <- unique(p$layers[[4]]$data[, c("target", ".sn_fill", "proportion", "tissue")])
  expect_equal(slices$proportion[slices$target == "T" & slices$tissue == "one"], c(.75, .25))
  expect_equal(slices$proportion[slices$target == "T" & slices$tissue == "two"], c(.2, .8))
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(length(unique(ggplot2::ggplot_build(p)$data[[1]]$PANEL)), 2L)
  no_pies <- sn_plot_sankey(d, flow_by = c("source", "target"), style = "C", y_by = "weight", show_pies = FALSE)
  expect_length(no_pies$layers, 5L)
  expect_error(sn_plot_sankey(d, flow_by = c("source", "target"), style = "C", fill_by = "target"), "first flow axis")
  expect_error(sn_plot_sankey(d, flow_by = c("source", "target"), style = "A", show_pies = TRUE), "requires style C")
})


test_that("Sankey drops zero-weight paths and rejects empty weighted plots", {
  d <- data.frame(source = c("Young", "Aged", "Young"), target = c("T", "T", "absent"), count = c(3, 1, 0))
  p <- sn_plot_sankey(d, flow_by = c("source", "target"), style = "C")
  expect_equal(p$layers[[3]]$data$label, "T")
  expect_no_error(ggplot2::ggplot_build(p))
  d$count <- 0
  expect_error(sn_plot_sankey(d, flow_by = c("source", "target"), style = "C"), "positive-weight")
})

test_that("style C pie size is constant across unequal facet compositions", {
  d <- data.frame(source = "Young", target = rep(c("A", "B", "C"), 2),
                  cohort = rep(c("one", "two"), each = 3), count = c(1000, 1, 1, 1, 1, 1))
  p <- sn_plot_sankey(d, flow_by = c("source", "target"), style = "C", facet_row_by = "cohort")
  pies <- p$layers[[4]]$data
  diameters <- tapply(pies$x, interaction(pies$target, pies$cohort), function(x) diff(range(x)))
  expect_equal(as.numeric(diameters), rep(as.numeric(diameters[1]), 6), tolerance = 1e-10)
  expect_no_error(ggplot2::ggplot_build(p))
})
