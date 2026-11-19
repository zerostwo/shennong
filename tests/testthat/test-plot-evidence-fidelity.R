test_that("raw probabilities stay distinct from adjusted probabilities in result plots", {
  de <- data.frame(gene = c("g1", "g2"), avg_log2FC = c(-1.8, 2.1), pvalue = c(.04, .09))
  volcano <- sn_plot_de(de)
  expect_match(volcano$labels$y, "raw p-value", fixed = TRUE)
  expect_equal(volcano$data$p_value, de$pvalue)
  expect_true(all(is.na(volcano$data$significant)))
  expect_false(any(vapply(volcano$layers, function(x) inherits(x$geom, "GeomHline"), logical(1))))
  expect_error(sn_plot_de(de, adjusted_p_value = 0.05), "adjusted p-value")
  de$padj <- p.adjust(de$pvalue, "BH")
  adjusted <- sn_plot_de(de)
  expect_equal(adjusted$data$adjusted_p_value, c(.08, .09))
  expect_equal(adjusted$data$significant, c(FALSE, FALSE))
  expect_s3_class(ggplot2::ggplotGrob(volcano), "gtable")

  enrichment <- data.frame(ID = c("P1", "P2"), pvalue = de$pvalue, Count = 1:2)
  for (plot in list(sn_plot_enrichment(enrichment), sn_plot_gsea(enrichment))) {
    expect_equal(plot$data$p_value, enrichment$pvalue)
    expect_match(tolower(plot$labels$colour), "raw p-value", fixed = TRUE)
    expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
  }
})

test_that("MA views require measured average expression and retain its scale", {
  de <- data.frame(gene = c("g1", "g2"), avg_log2FC = c(-1.8, 2.1), pvalue = c(.04, .09))
  expect_error(sn_plot_de(de, "ma"), "average expression")
  de$pct.1 <- c(.1, .9)
  expect_error(sn_plot_de(de, "ma"), "average expression")
  de$logCPM <- c(9.55, 2.18)
  plot <- sn_plot_de(de, "ma")
  expect_equal(ggplot2::ggplot_build(plot)$data[[1]]$x, de$logCPM)
  expect_match(plot$labels$x, "logCPM", fixed = TRUE)
})

test_that("different enrichment records with equal labels retain distinct positions", {
  e <- data.frame(ID = c("P1", "P2"), Description = c("glycolysis", "glycolysis"),
                  p.adjust = c(.02, .03), Count = c(5, 10))
  for (p in list(sn_plot_enrichment(e, "dot"), sn_plot_enrichment(e, "bar"), sn_plot_gsea(e))) {
    build <- ggplot2::ggplot_build(p)
    expect_equal(length(unique(build$data[[1]]$y)), 2L)
    expect_setequal(build$layout$panel_params[[1]]$y$get_labels(), c("glycolysis [P1]", "glycolysis [P2]"))
    expect_setequal(p$data$id, e$ID)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  }
})

test_that("heatmaps derive measurements from the selected assay and layer", {
  skip_if_not_installed("Seurat")
  rna <- matrix(1:16, 2L, dimnames = list(c("g1", "g2"), paste0("c", 1:8)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(rna, sparse = TRUE))
  protein <- matrix(1:8, 2L, dimnames = list(c("CD4", "CD8"), paste0("c", 3:6)))
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = Matrix::Matrix(protein, sparse = TRUE))
  object$CD4 <- 999
  object$group <- rep(c("x", "y"), each = 4)
  plot <- sn_plot_heatmap(object, features = "CD4", assay = "ADT", layer = "counts",
                          mode = "average", average_scale = FALSE, group_by = "group")
  expect_equal(as.character(plot$data$feature), c("CD4", "CD4"))
  expect_equal(plot$data$expression, c(2, 6))
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
})

test_that("communication views draw autocrine interactions and retain their scores", {
  links <- data.frame(source = c("A", "A"), target = c("A", "B"),
                      ligand = c("L1", "L2"), receptor = c("R1", "R2"), score = 1:2)
  for (kind in c("network", "chord")) {
    for (edges in list(links, links[1, ])) {
      plot <- sn_plot_communication(edges, type = kind)
      expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
      loop <- Filter(function(layer) inherits(layer$geom, "GeomPath"), plot$layers)
      expect_length(loop, 1L)
      expect_equal(unique(loop[[1]]$data$score), 1)
    }
  }
})

test_that("PDF export embeds fonts on request and a heatmap uses a single page", {
  skip_if_not(Shennong:::.sn_has_cairo_pdf())
  dir <- withr::local_tempdir()
  # Search PDF resource dictionaries in the uncompressed, printable byte runs.
  has_token <- function(path, pattern) {
    bytes <- readBin(path, "raw", file.size(path))
    bytes[bytes < as.raw(32L) | bytes > as.raw(126L)] <- charToRaw(" ")
    grepl(pattern, rawToChar(bytes))
  }
  plot <- ggplot2::ggplot(iris, ggplot2::aes(.data$Sepal.Length, .data$Sepal.Width)) +
    ggplot2::geom_point()
  embedded <- sn_save_figure(plot, file.path(dir, "embedded.pdf"))
  plain <- sn_save_figure(plot, file.path(dir, "portable.pdf"), embed_fonts = FALSE)
  expect_true(has_token(embedded, "/FontFile[23]? "))
  expect_false(has_token(plain, "/FontFile[23]? "))
  skip_if_not_installed("ComplexHeatmap")
  heatmap <- sn_save_figure(ComplexHeatmap::Heatmap(matrix(1:9, 3)), file.path(dir, "heatmap.pdf"))
  expect_true(has_token(heatmap, "/FontFile[23]? "))
  plain_heatmap <- sn_save_figure(ComplexHeatmap::Heatmap(matrix(1:9, 3)),
                                  file.path(dir, "portable-heatmap.pdf"), embed_fonts = FALSE)
  expect_true(has_token(plain_heatmap, "/Count +1\\b"))
})


test_that("headless PDF fallback reports the device actually used", {
  local_mocked_bindings(.sn_has_cairo_pdf = function() FALSE, .package = "Shennong")
  plot <- ggplot2::ggplot(iris, ggplot2::aes(.data$Sepal.Length, .data$Sepal.Width)) +
    ggplot2::geom_point()
  path <- tempfile(fileext = ".pdf")
  expect_warning(result <- sn_save_figure(plot, path), "without embedded fonts")
  expect_true(file.exists(path))
  expect_identical(attr(result, "embed_fonts"), FALSE)
})

test_that("incomplete communication endpoints do not hide available links", {
  edges <- data.frame(source = c("A", NA, "A"), target = c("B", "A", NA),
                      ligand = "L", receptor = "R", score = 1:3)
  for (kind in c("network", "chord")) {
    plot <- sn_plot_communication(edges, type = kind)
    # The geometry can omit incomplete coordinates when it draws.
    grob <- suppressWarnings(ggplot2::ggplotGrob(plot))
    expect_s3_class(grob, "gtable")
    edge_layer <- Filter(function(layer) inherits(layer$geom, c("GeomSegment", "GeomCurve")), plot$layers)
    expect_length(edge_layer, 1L)
    links <- edge_layer[[1]]$data
    valid <- links[stats::complete.cases(links[c("x_source", "y_source", "x_target", "y_target")]), ]
    expect_equal(as.character(valid$source), "A")
    expect_equal(as.character(valid$target), "B")
    expect_equal(valid$score, 1)
  }
})
