test_that("Sankey preserves independent factor order and quantitative widths", {
  d <- data.frame(a = factor(c("A", "Z", "A"), c("Z", "A", "unused")),
                  b = factor(c("Z", "A", "B"), c("B", "A", "Z")), count = c(2, 4, 1))
  p <- sn_plot_composition(d, type = "sankey", flow_by = c("a", "b"),
                           fill = a, data_kind = "summary")
  nodes <- p$layers[[2]]$data
  expect_equal(nodes$label[nodes$axis == 1], c("Z", "A"))
  expect_equal(nodes$label[nodes$axis == 2], c("B", "A", "Z"))
  expect_equal(nodes$ymax - nodes$ymin, c(4, 3, 1, 4, 2))
  expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomRect"), logical(1))))
  expect_no_error(ggplot2::ggplot_build(p))
  expanded <- sn_plot_composition(d, type = "sankey", flow_by = c("a", "b"),
                                  sankey_layout = "expanded", show_stratum_boxes = TRUE)
  nodes <- expanded$layers[[2]]$data
  height <- tapply(nodes$ymin, nodes$axis, min)
  expect_lt(height[[2]], height[[1]])
  expect_equal(nodes$ymax - nodes$ymin, c(4, 3, 1, 4, 2))
  expect_no_error(ggplot2::ggplot_build(expanded))
})

test_that("Sankey validates weights and spacing", {
  d <- data.frame(a = "A", b = "B", count = -1)
  expect_error(sn_plot_composition(d, type = "sankey", flow_by = c("a", "b")), "non-negative")
  d$count <- 1
  expect_error(sn_plot_composition(d, type = "sankey", flow_by = c("a", "b"),
                                   stratum_gap = -1), "stratum_gap")
})

test_that("Sankey uses Seurat factor levels without modifying metadata", {
  skip_if_not_installed("SeuratObject")
  counts <- matrix(1:12, 3, dimnames = list(paste0("g", 1:3), paste0("c", 1:4)))
  obj <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  obj$a <- factor(c("A", "Z", "A", "Z"), levels = c("Z", "A"))
  obj$b <- factor(c("Z", "A", "B", "Z"), levels = c("B", "A", "Z"))
  before <- obj[[]]
  p <- sn_plot_composition(obj, type = "sankey", flow_by = c("a", "b"))
  nodes <- p$layers[[2]]$data
  expect_equal(nodes$label[nodes$axis == 1], c("Z", "A"))
  expect_equal(nodes$label[nodes$axis == 2], c("B", "A", "Z"))
  expect_identical(obj[[]], before)
  expect_equal(as.numeric(tapply(nodes$ymax - nodes$ymin, nodes$axis, sum)), c(4, 4))
})

test_that("Sankey supports facets, missing paths and independently ordered stages", {
  d <- data.frame(a = c("Z", "A", "Z", NA), b = c("X", "Y", "Y", "X"),
                  c = c("Q", "P", "P", "Q"), cohort = c("one", "one", "two", "two"))
  p <- sn_plot_composition(d, type = "sankey", flow_by = c("a", "b", "c"),
                           facet_col = cohort, sankey_layout = "expanded")
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$PANEL)), 2L)
  nodes <- p$layers[[2]]$data
  expect_equal(sum(nodes$ymax[nodes$axis == 1] - nodes$ymin[nodes$axis == 1]), 3)
  expect_equal(nodes$label[nodes$axis == 1 & nodes$cohort == "one"], c("Z", "A"))
  hidden <- sn_plot_composition(d, type = "sankey", flow_by = c("a", "b"),
                                show_stratum_labels = FALSE, stratum_gap = 0)
  expect_length(hidden$layers, 1L)
  expect_no_error(ggplot2::ggplot_build(hidden))
})
