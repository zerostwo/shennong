library(testthat)

test_that("sn_write and sn_read round-trip a csv file", {
  input <- data.frame(a = 1:3, b = c("x", "y", "z"))
  path <- tempfile(fileext = ".csv")

  sn_write(input, path)
  output <- sn_read(path)

  expect_s3_class(output, "data.frame")
  expect_equal(output$a, input$a)
  expect_equal(output$b, input$b)
})

test_that("plot helpers return ggplot objects", {
  data <- data.frame(
    group = c("A", "B", "A", "B"),
    value = c(1, 2, 3, 4),
    fill_group = c("x", "x", "y", "y")
  )

  expect_s3_class(sn_plot_boxplot(data, x = group, y = value), "ggplot")
  expect_s3_class(sn_plot_barplot(data, x = group, y = value, fill = fill_group), "ggplot")
})

test_that("show_all_palettes prints without error", {
  expect_no_error(show_all_palettes())
})

test_that("sn_add_data_from_anndata adds metadata and embeddings", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(40, lambda = 3), nrow = 10, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(10))
  colnames(counts) <- paste0("cell", seq_len(4))
  object <- sn_initialize_seurat_object(x = counts, project = "ann")

  metadata <- data.frame(cell_type = c("T", "T", "B", "B"), row.names = colnames(object))
  meta_path <- tempfile(fileext = ".csv")
  write.csv(metadata, meta_path)

  umap <- data.frame(UMAP_1 = c(1, 2, 3, 4), UMAP_2 = c(4, 3, 2, 1), row.names = colnames(object))
  umap_path <- tempfile(fileext = ".csv")
  write.csv(umap, umap_path)

  updated <- sn_add_data_from_anndata(
    object = object,
    metadata_path = meta_path,
    umap_path = umap_path
  )

  expect_true("cell_type" %in% colnames(updated[[]]))
  expect_true("umap" %in% names(updated@reductions))
})
