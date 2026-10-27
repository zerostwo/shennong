atlas_fixture <- function() {
  set.seed(44)
  m <- matrix(rpois(8 * 36, 3), 8, dimnames = list(paste0("g", 1:8), paste0("c", 1:36)))
  x <- Seurat::CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE))
  x$celltype <- factor(rep(c("NK", "B", "Myeloid"), each = 12))
  x$sample <- rep(paste0("s", 1:6), each = 6)
  x$condition <- rep(c("Normal", "Tumour"), 18)
  x <- Seurat::NormalizeData(x, verbose = FALSE)
  x[["umap"]] <- SeuratObject::CreateDimReducObject(
    matrix(rnorm(72), ncol = 2, dimnames = list(colnames(x), c("UMAP_1", "UMAP_2"))),
    key = "UMAP_", assay = "RNA")
  x
}

test_that("atlas styles preserve marker statistics and expression values", {
  x <- atlas_fixture()
  base <- sn_plot_dot(x, features = c("g1", "g2"), group_by = "celltype", scale = FALSE)
  for (style in c("atlas", "atlas_outline", "atlas_module")) {
    p <- sn_plot_dot(x, features = c("g1", "g2"), group_by = "celltype", scale = FALSE, style = style)
    expect_equal(p$data, base$data)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
    expect_equal(attr(p, "shennong_plot_style")$name, style)
  }
  p <- sn_plot_dot(x, features = c("g1", "g2"), group_by = "celltype", scale = FALSE, style = "atlas_outline")
  expect_equal(ggplot2::ggplot_build(p)$data[[1]]$shape, rep(21, 6))
  plain <- sn_plot_heatmap(x, c("g1", "g2"), group_by = "celltype", layer = "data", mode = "average")
  styled <- sn_plot_heatmap(x, c("g1", "g2"), group_by = "celltype", layer = "data", mode = "average", style = "atlas")
  expect_equal(plain$data, styled$data)
  expect_s3_class(ggplot2::ggplotGrob(styled), "gtable")
})

test_that("presentation styles do not change samples, denominators or correlations", {
  d <- data.frame(sample = rep(letters[1:4], each = 8),
    condition = rep(c("N", "T"), each = 16), type = rep(c("A", "A", "A", "B"), 8))
  a <- sn_plot_sample_boxplot(d, "condition", "type", "sample")
  b <- sn_plot_sample_boxplot(d, "condition", "type", "sample", style = "atlas_outline")
  expect_equal(attr(a, "shennong_figure_data"), attr(b, "shennong_figure_data"))
  expect_s3_class(ggplot2::ggplotGrob(b), "gtable")
  d <- data.frame(x = 1:12, y = c(2, 4, 3, 5, 7, 6, 9, 8, 10, 12, 11, 13), group = rep(c("A", "B"), 6))
  a <- sn_plot_association(d, "x", "y", group_by = "group")
  b <- sn_plot_association(d, "x", "y", group_by = "group", style = "atlas")
  expect_equal(attr(a, "shennong_figure_data"), attr(b, "shennong_figure_data"))
  expect_identical(a$labels$subtitle, b$labels$subtitle)
  expect_error(sn_plot_association(d, "x", "y", style = "typo"), "style")
})

test_that("matrix plots retain missing values and require aligned bubble sizes", {
  m <- matrix(c(-1, NA, 2, 0, 1, 3), 2, dimnames = list(c("A", "B"), c("Lung", "Liver", "Skin")))
  p <- sn_plot_matrix(m, style = "atlas_module", midpoint = 0, show_values = TRUE)
  expect_equal(attr(p, "shennong_figure_data")$matrix, m)
  expect_equal(nrow(p$data), length(m))
  expect_equal(sum(is.na(p$data$value)), 1)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  sizes <- m; sizes[] <- seq(0, 100, length.out = 6)
  p <- sn_plot_matrix(m, type = "bubble", size = sizes, style = "atlas_outline")
  expect_equal(attr(p, "shennong_figure_data")$size, sizes)
  expect_equal(ggplot2::ggplot_build(p)$data[[1]]$size[1], 0)
  expect_equal(ggplot2::ggplot_build(p)$data[[1]]$stroke[1], 0)
  expect_error(sn_plot_matrix(m, type = "bubble", size = sizes[2:1, ]), "dimnames")
  expect_error(sn_plot_matrix(m, cluster_rows = TRUE), "missing|finite")
})

test_that("network renderer preserves supplied nodes and edge weights", {
  edges <- data.frame(source = c("A", "B"), target = c("B", "C"), weight = c(-0.4, 0.8))
  nodes <- data.frame(node = c("A", "B", "C", "isolated"), lineage = c("T", "B", "Myeloid", "T"))
  p <- sn_plot_network(edges, nodes, group_by = "lineage", style = "atlas_module")
  expect_equal(attr(p, "shennong_figure_data")$edges, edges)
  expect_equal(nrow(p$data), 4)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  expect_error(sn_plot_network(edges, nodes[-1, ]), "missing|Missing")
  expect_error(sn_plot_network(rbind(edges, edges[1, ]), nodes), "Duplicate|duplicate")
  expect_error(sn_plot_network(edges, nodes, layout = "supplied"), "x.*y")
})

test_that("atlas embedding styles preserve coordinates and explicit colors", {
  x <- atlas_fixture()
  palette <- c(NK = "#FF0000", B = "#00FF00", Myeloid = "#0000FF")
  plain <- sn_plot_dim(x, group_by = "celltype", raster = FALSE)
  styled <- sn_plot_dim(x, group_by = "celltype", raster = FALSE,
    style = "atlas", cols = palette)
  expect_equal(styled$data, plain$data)
  expect_equal(ggplot2::ggplot_build(styled)$data[[1]]$colour,
    unname(palette[as.character(styled$data$celltype)]))
  expect_true(inherits(styled$theme$axis.text, "element_blank"))
  split <- sn_plot_dim(x, group_by = "celltype", split_by = "condition", combine = FALSE,
    raster = FALSE, style = "atlas_outline")
  expect_type(split, "list")
  expect_true(all(vapply(split, inherits, logical(1), "ggplot")))
  plain <- sn_plot_feature(x, "g1", raster = FALSE)
  styled <- sn_plot_feature(x, "g1", raster = FALSE, style = "atlas_module")
  expect_equal(plain$data, styled$data)
  expect_s3_class(ggplot2::ggplotGrob(styled), "gtable")
  expect_error(sn_plot_dot(x, "g1", split_by = "condition", style = "atlas"), "split_by")
  expect_error(sn_plot_heatmap(x, "g1", mode = "average", backend = "pheatmap",
    style = "atlas"), "ggplot2")
})

test_that("new renderers export original matrices and network tables", {
  m <- matrix(c(-1, NA, 2, 0), 2, dimnames = list(c("A", "B"), c("T1", "T2")))
  p <- sn_plot_matrix(m, style = "atlas", limits = c(0, 1))
  directory <- file.path(withr::local_tempdir(), "matrix")
  sn_export_figure_bundle(p, directory, formats = "pdf")
  readback <- as.matrix(utils::read.csv(file.path(directory, "matrix_data.csv"),
    check.names = FALSE, row.names = 1))
  expect_equal(readback, m)
  expect_equal(Shennong:::.sn_figure_source_data(p)$matrix$.shennong_row, rownames(m))
  edges <- data.frame(source = "A", target = "B", weight = 0.5)
  nodes <- data.frame(node = c("B", "A", "isolated"), x = c(1, 2, 3), y = c(0, 4, 2))
  p <- sn_plot_network(edges, nodes, layout = "supplied", directed = TRUE, style = "atlas")
  tables <- Shennong:::.sn_figure_source_data(p)
  expect_equal(tables$edges, edges)
  expect_equal(tables$nodes, nodes)
  expect_equal(tables$coordinates, nodes)
  expect_equal(sn_get_figure_spec(p)$data_summary$n_edges, 1)
  endpoint <- ggplot2::ggplot_build(p)$data[[1]]
  expect_gt(endpoint$xend, nodes$x[nodes$node == "B"])
  expect_lt(endpoint$xend, nodes$x[nodes$node == "A"])
})

test_that("atlas style retains DE, enrichment and communication statistics", {
  de <- data.frame(gene = paste0("g", 1:4), log2FoldChange = c(-2, 1, 0.1, 3),
    padj = c(0.02, 0.04, 0.8, 0.001), baseMean = c(10, 20, 30, 40))
  enrichment <- data.frame(ID = c("P1", "P2"), Description = c("Pathway 1", "Pathway 2"),
    p.adjust = c(0.01, 0.04), NES = c(-2, 1.4), core_enrichment = c("A/B", "B/C"))
  for (type in c("volcano", "ma", "effect", "heatmap")) {
    p <- sn_plot_de(de, type = type, style = "atlas")
    expect_equal(p$data, sn_plot_de(de, type = type)$data)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  }
  for (type in c("dot", "bar", "network", "emap")) {
    p <- sn_plot_enrichment(enrichment, type = type, style = "atlas_outline")
    expect_equal(p$data, sn_plot_enrichment(enrichment, type = type)$data)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  }
  interactions <- data.frame(source = c("A", "B"), target = c("B", "C"),
    ligand = c("L1", "L2"), receptor = c("R1", "R2"), score = c(0.4, 0.8))
  for (type in c("bubble", "heatmap", "network", "chord", "river")) {
    p <- sn_plot_communication(interactions, type = type, style = "atlas_outline")
    expect_equal(attr(p, "shennong_figure_data"),
      attr(sn_plot_communication(interactions, type = type), "shennong_figure_data"))
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  }
})

test_that("spatial atlas styles preserve sections and raw expression", {
  x <- atlas_fixture()
  x$x <- rep(1:6, 6); x$y <- rep(1:6, each = 6)
  plain <- sn_plot_spatial_feature(x, "g1", split_by = "condition", layer = "counts")
  p <- sn_plot_spatial_feature(x, "g1", split_by = "condition", layer = "counts", style = "atlas")
  expect_equal(attr(p, "shennong_figure_data"), attr(plain, "shennong_figure_data"))
  expect_equal(p$coordinates$ratio, plain$coordinates$ratio)
  expect_equal(p$theme$text$size, 8)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})


test_that("summary bars retain sample means and SD across styles", {
  d <- data.frame(group = rep(c("A", "B"), each = 5), value = c(1:5, 4:8))
  plain <- sn_plot_barplot(d, group, value, stat = "summary", errorbar = "sd")
  p <- sn_plot_barplot(d, group, value, stat = "summary", errorbar = "sd", style = "atlas")
  expect_equal(attr(p, "shennong_figure_data"), attr(plain, "shennong_figure_data"))
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  p <- sn_plot_boxplot(d, group, value, style = "atlas_outline")
  expect_equal(attr(p, "shennong_figure_data"), d)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
})
