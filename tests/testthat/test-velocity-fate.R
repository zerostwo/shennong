library(testthat)

make_velocity_test_object <- function() {
  set.seed(23)
  counts <- matrix(rpois(10 * 30, 3), nrow = 10)
  rownames(counts) <- paste0("G", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  embedding <- cbind(seq(-2, 2, length.out = 30), sin(seq(0, 2 * pi, length.out = 30)))
  rownames(embedding) <- colnames(object)
  colnames(embedding) <- c("UMAP_1", "UMAP_2")
  object[["umap"]] <- SeuratObject::CreateDimReducObject(embedding, key = "UMAP_", assay = "RNA")
  object$cell_type <- rep(c("early", "middle", "late"), each = 10)
  object
}

velocity_adapter_result <- function(object) {
  embedding <- SeuratObject::Embeddings(object[["umap"]])
  list(
    cells = tibble::tibble(
      cell = colnames(object), velocity_x = rep(0.1, ncol(object)),
      velocity_y = c(diff(embedding[, 2]), 0), velocity_pseudotime = seq(0, 1, length.out = ncol(object)),
      velocity_confidence = seq(0.5, 0.9, length.out = ncol(object))
    ),
    graph = tibble::tibble(
      source = head(colnames(object), -1), target = tail(colnames(object), -1),
      probability = rep(0.8, ncol(object) - 1)
    ),
    artifacts = list(output_h5ad = "synthetic.h5ad")
  )
}

test_that("scVelo adapter stores projected vectors and diagnostics", {
  object <- make_velocity_test_object()
  updated <- sn_run_velocity(object, backend_control = list(result = velocity_adapter_result(object)))
  result <- sn_get_result(updated, "velocity", "velocity")
  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_equal(nrow(result$tables$cells), ncol(object))
  expect_equal(nrow(result$tables$transition_edges), ncol(object) - 1L)
  expect_equal(result$diagnostics$finite_vectors, ncol(object))
  expect_true(all(c("velocity_pseudotime", "velocity_confidence") %in% colnames(updated[[]])))
})

test_that("CellRank adapter stores probabilities, terminals, and drivers", {
  object <- make_velocity_test_object()
  probabilities <- expand.grid(cell = colnames(object), state = c("A", "B"), stringsAsFactors = FALSE)
  probabilities$probability <- rep(c(0.7, 0.3), each = ncol(object))
  output <- list(
    probabilities = probabilities,
    terminal_states = tibble::tibble(cell = c("cell1", "cell30"), state = c("A", "B"), probability = 1),
    drivers = tibble::tibble(feature = c("G1", "G2"), state = c("A", "B"), correlation = c(0.8, 0.7))
  )
  updated <- sn_run_fate(object, backend_control = list(result = output))
  result <- sn_get_result(updated, "fate", "fate")
  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_equal(result$diagnostics$states, 2L)
  expect_equal(nrow(result$tables$probabilities), 2L * ncol(object))
  expect_true(all(c("fate_fate_A", "fate_fate_B") %in% colnames(updated[[]])))
})

test_that("velocity and fate runners receive backend context", {
  object <- make_velocity_test_object()
  velocity_runner <- function(object, method, spliced_assay, spliced_layer, unspliced_assay, unspliced_layer, reduction, dims, backend_control) {
    expect_equal(method, "scvelo")
    expect_equal(reduction, "umap")
    velocity_adapter_result(object)
  }
  velocity <- sn_run_velocity(object, backend_control = list(runner = velocity_runner), return_object = FALSE)
  expect_equal(velocity$method, "scvelo")

  fate_runner <- function(object, method, velocity_result, reduction, dims, backend_control) {
    list(probabilities = tibble::tibble(cell = colnames(object), state = "terminal", probability = 1))
  }
  fate <- sn_run_fate(object, backend_control = list(runner = fate_runner), return_object = FALSE)
  expect_equal(fate$method, "cellrank")
})

test_that("velocity and fate plots render", {
  object <- make_velocity_test_object()
  velocity <- sn_run_velocity(object, backend_control = list(result = velocity_adapter_result(object)), return_object = FALSE)
  fate <- sn_run_fate(
    object,
    backend_control = list(result = list(probabilities = tibble::tibble(
      cell = rep(colnames(object), 2), state = rep(c("A", "B"), each = ncol(object)),
      probability = rep(seq(0, 1, length.out = ncol(object)), 2)
    ))), return_object = FALSE
  )
  for (plot in list(sn_plot_velocity(velocity), sn_plot_fate(fate))) {
    expect_s3_class(plot, "ggplot")
    expect_silent(ggplot2::ggplotGrob(plot))
  }
})

test_that("trajectory pixi environment and dynamics methods are registered", {
  expect_true("trajectory" %in% sn_list_pixi_environments())
  expect_true(sn_method_status("scvelo", "velocity")$implemented)
  expect_true(sn_method_status("cellrank", "fate")$implemented)
})
