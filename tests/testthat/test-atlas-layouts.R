test_that("triangle layouts mask geometry without changing numerical sources", {
  m <- matrix(1:16, 4, dimnames = list(letters[1:4], letters[1:4]))
  m <- m * 1
  p <- sn_plot_matrix(m, style = "upper_triangle")
  expect_equal(nrow(p$data), 6)
  expect_true(all(p$data$.row < p$data$.col))
  expect_equal(attr(p, "shennong_figure_data")$matrix, m)
  p <- sn_plot_matrix(m, style = "split_triangle", compare = m + 20)
  expect_equal(p$data$value[p$data$.row > p$data$.col], (m + 20)[lower.tri(m)])
  expect_equal(attr(p, "shennong_figure_data")$compare, m + 20)
  expect_error(sn_plot_matrix(m, style = "upper_triangle", cluster_rows = TRUE), "order|cluster")
  expect_error(sn_plot_matrix(m, style = "split_triangle"), "compare")
})

test_that("matrix composites retain annotation alignment and summaries explicitly", {
  m <- matrix(seq_len(12) / 10, 3, dimnames = list(LETTERS[1:3], letters[1:4]))
  ann <- data.frame(tissue = c("Skin", "Lung", "Skin", "Lung"), row.names = rev(letters[1:4]))
  p <- sn_plot_matrix(m, style = "annotated", column_annotation = ann)
  expect_s3_class(p, "patchwork")
  expect_equal(attr(p, "shennong_figure_data")$column_annotation, ann[colnames(m), , drop = FALSE])
  for (s in c("circular", "dendrogram")) {
    p <- sn_plot_matrix(m, style = s)
    expect_s3_class(if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplot2::ggplotGrob(p), "gtable")
    expect_equal(attr(p, "shennong_figure_data")$matrix, m)
  }
  expect_error(sn_plot_matrix(m, style = "heatmap_bar"), "row_summary")
  p <- sn_plot_matrix(m, style = "heatmap_bar", row_summary = c(C = 6, B = 4, A = 2))
  expect_equal(attr(p, "shennong_figure_data")$row_summary$value, c(2, 4, 6))
})

test_that("paired layouts require unique biological units and complete pairs", {
  d <- data.frame(donor = rep(letters[1:4], 2), group = rep(c("N", "T"), each = 4), value = 1:8)
  p <- sn_plot_distribution(d, "value", "group", style = "paired_box", pair_by = "donor")
  expect_equal(attr(p, "shennong_figure_data")$values$pair, d$donor)
  expect_equal(length(unique(ggplot2::ggplot_build(p)$data[[1]]$group)), 4)
  expect_error(sn_plot_distribution(d, "value", "group", style = "paired_slope"), "pair_by")
  expect_error(sn_plot_distribution(rbind(d, d[1, ]), "value", "group", style = "paired_box", pair_by = "donor"), "one|Duplicate|duplicate")
  expect_error(sn_plot_distribution(d[-1, ], "value", "group", style = "paired_slope", pair_by = "donor"), "complete|missing")
  cells <- d[rep(seq_len(nrow(d)), each = 3), ]; cells$sample <- rep(paste0("s", 1:8), each = 3)
  p <- sn_plot_distribution(cells, "value", "group", sample_by = "sample", pair_by = "donor", style = "paired_slope")
  expect_equal(nrow(attr(p, "shennong_figure_data")$values), 8)
})

test_that("distribution styles have distinct geometries and preserve observations", {
  d <- data.frame(group = rep(c("Pre", "Post"), each = 20), value = sin(1:40) + rep(c(0, 1), each = 20))
  for (style in c("violin_box", "raincloud", "split_violin", "stacked_violin")) {
    p <- sn_plot_distribution(d, "value", "group", style = style)
    expect_equal(attr(p, "shennong_figure_data")$values$value, d$value)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  }
  d$group <- rep(letters[1:4], 10)
  expect_error(sn_plot_distribution(d, "value", "group", style = "split_violin"), "two")
})

test_that("effect plots distinguish coefficients from ratios", {
  d <- data.frame(label = c("A", "B"), estimate = c(-1, 2), lower = c(-2, 1), upper = c(0, 3))
  p <- sn_plot_effect(d)
  expect_equal(attr(p, "shennong_figure_data")$effects, d)
  expect_equal(ggplot2::ggplot_build(p)$data[[1]]$xintercept, 0)
  expect_error(sn_plot_effect(d, scale = "ratio"), "positive")
  d$lower[1] <- 2
  expect_error(sn_plot_effect(d), "interval|lower")
})

test_that("UpSet counts exclusive membership and retains original set IDs", {
  sets <- list(B = c("s1", "s2", "s3"), T = c("s2", "s4"), Myeloid = c("s2", "s3"))
  p <- sn_plot_overlap(sets)
  tab <- attr(p, "shennong_figure_data")$intersections
  expect_equal(sum(tab$count), 4)
  expect_equal(tab$count[tab$intersection == "B&T&Myeloid"], 1)
  expect_s3_class(patchwork::patchworkGrob(p), "gtable")
  expect_error(sn_plot_overlap(list(B = c("s1", "s1"), T = "s1")), "duplicate|Duplicate")
})

test_that("network layouts preserve topology and reject misleading trees/chords", {
  e <- data.frame(source = c("A", "A", "B"), target = c("B", "C", "D"), weight = c(1, 2, 3))
  n <- data.frame(node = LETTERS[1:4], side = c("L", "R", "R", "L"))
  for (style in c("bipartite", "clone_tree", "chord", "rewiring")) {
    p <- sn_plot_network(e, n, group_by = "side", style = style,
      compare_edges = if (style == "rewiring") e[1:2, ] else NULL)
    expect_equal(attr(p, "shennong_figure_data")$edges, e)
    expect_s3_class(if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else ggplot2::ggplotGrob(p), "gtable")
  }
  cycle <- rbind(e, data.frame(source = "D", target = "A", weight = 1))
  expect_error(sn_plot_network(cycle, n, style = "clone_tree"), "tree|root|cycle")
  e$weight[1] <- -1
  expect_error(sn_plot_network(e, n, style = "chord"), "non-negative")
})

test_that("existing Seurat heatmap and violin interfaces dispatch structural styles", {
  x <- SeuratObject::pbmc_small
  p <- sn_plot_heatmap(x, c("MS4A1", "CD79A"), layer = "data", mode = "average", style = "dendrogram")
  expect_s3_class(p, "patchwork")
  expect_true(all(c("matrix", "averages") %in% names(attr(p, "shennong_figure_data"))))
  expect_s3_class(patchwork::patchworkGrob(p), "gtable")
  p <- sn_plot_violin(x, c("MS4A1", "CD79A"), style = "raincloud")
  expect_equal(nrow(attr(p, "shennong_figure_data")$values), 160)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})

test_that("chord widths and rewiring coordinates retain quantitative meaning", {
  e <- data.frame(source = c("A", "A"), target = c("B", "C"), weight = c(1, 3))
  p <- sn_plot_network(e, style = "chord")
  a <- attr(p, "shennong_figure_data")$ribbon_endpoints
  expect_equal(diff(a$source_end - a$source_start) / (a$source_end[1] - a$source_start[1]), 2)
  p <- sn_plot_network(e, style = "rewiring", compare_edges = e[2, ])
  expect_equal(attr(p, "shennong_figure_data")$compare_edges, e[2, ])
  first <- p$patches$plots[[1]]
  expect_equal(first$data[c("node", "x", "y")], p$data[c("node", "x", "y")])
  expect_equal(ggplot2::ggplot_build(first)$plot$scales$get_scales("colour")$limits,
               ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$limits)
})
