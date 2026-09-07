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
  expect_true(result$parameters$enforce_normalization)
  expect_true(result$parameters$log1p_transform)
  expect_identical(result$backend, "scvelo-provided-result")
  expect_true(all(c("velocity_pseudotime", "velocity_confidence") %in% colnames(updated[[]])))
})

test_that("managed velocity runs expose explicit retention and remove raw exports", {
  object <- make_velocity_test_object()
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  SeuratObject::LayerData(object, assay = "RNA", layer = "spliced") <- counts
  SeuratObject::LayerData(object, assay = "RNA", layer = "unspliced") <- counts
  embedding <- Shennong:::.sn_velocity_embedding(object, "umap", 1:2)
  observed_roots <- character()
  runner <- function(environment, command, args, ...) {
    input_dir <- args[[match("--input-dir", args) + 1L]]
    output_dir <- args[[match("--output-dir", args) + 1L]]
    observed_roots <<- c(observed_roots, dirname(input_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cells <- colnames(object)
    utils::write.csv(
      data.frame(
        cell = cells,
        velocity_1 = 0.1,
        velocity_2 = 0.2,
        pseudotime = seq(0, 1, length.out = length(cells)),
        confidence = 0.8
      ),
      file.path(output_dir, "velocity_cells.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      data.frame(source = head(cells, -1L), target = tail(cells, -1L), probability = 0.9),
      file.path(output_dir, "velocity_graph.csv"),
      row.names = FALSE
    )
    h5ad <- file.path(output_dir, "velocity.h5ad")
    file.create(h5ad)
    jsonlite::write_json(
      list(method = "scvelo", output_h5ad = h5ad, n_cells = length(cells)),
      file.path(output_dir, "manifest.json"),
      auto_unbox = TRUE
    )
  }
  runtime_dir <- withr::local_tempdir(pattern = "velocity-runtime-")

  retained <- testthat::with_mocked_bindings(
    Shennong:::.sn_run_velocity_pixi(
      object, "scvelo", "RNA", "spliced", "RNA", "unspliced", embedding,
      list(runtime_dir = runtime_dir)
    ),
    sn_call_pixi_environment = runner,
    .package = "Shennong"
  )
  retained_root <- retained$artifacts$run_dir
  withr::defer({
    if (dir.exists(retained_root)) Shennong:::.sn_cleanup_owned_run_dir(retained_root)
  })
  expect_true(retained$artifacts$run_dir_retained)
  expect_true(file.exists(retained$artifacts$output_h5ad))
  expect_false(dir.exists(file.path(retained_root, "input")))
  expect_false(file.exists(file.path(retained_root, "config.json")))

  cleanup_parent <- withr::local_tempdir(pattern = "velocity-clean-parent-")
  sentinel <- file.path(cleanup_parent, "user-sentinel.txt")
  writeLines("keep", sentinel)
  temporary <- testthat::with_mocked_bindings(
    Shennong:::.sn_run_velocity_pixi(
      object, "scvelo", "RNA", "spliced", "RNA", "unspliced", embedding,
      list(run_dir = cleanup_parent, keep_run_dir = FALSE)
    ),
    sn_call_pixi_environment = runner,
    .package = "Shennong"
  )
  expect_true(file.exists(sentinel))
  expect_false(dir.exists(utils::tail(observed_roots, 1L)))
  expect_false(temporary$artifacts$output_retained)
  expect_null(temporary$artifacts$output_h5ad)
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

test_that("RegVelo reuses the velocity result contract", {
  object <- make_velocity_test_object()
  result <- sn_run_velocity(
    object,
    method = "regvelo",
    backend_control = list(result = velocity_adapter_result(object)),
    return_object = FALSE
  )

  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_identical(result$method, "regvelo")
  expect_identical(result$backend, "regvelo-provided-result")
  expect_true(all(result$tables$cells$method == "regvelo"))
  expect_true(all(c("soft_constraint", "lam", "max_epochs", "enforce_normalization", "log1p_transform") %in% names(result$parameters)))
})

test_that("managed velocity preprocessing enforces and records normalization", {
  script <- paste(
    readLines(Shennong:::.sn_trajectory_pixi_script(), warn = FALSE),
    collapse = "\n"
  )

  expect_match(script, 'config.get("enforce_normalization", True)', fixed = TRUE)
  expect_match(script, "enforce=enforce_normalization", fixed = TRUE)
  expect_match(script, "sc.pp.log1p(adata)", fixed = TRUE)
  expect_match(script, '"velocity_mode": config.get', fixed = TRUE)
  expect_match(script, '"preprocessing": {', fixed = TRUE)
  expect_match(script, "estimator.compute_eigendecomposition()", fixed = TRUE)
})

test_that("velocity inference dimensions are independent of the 2D plotting embedding", {
  exported <- list(features = paste0("G", seq_len(80)), cells = paste0("C", seq_len(50)))
  expect_identical(Shennong:::.sn_velocity_inference_n_pcs(exported), 30L)
  expect_identical(Shennong:::.sn_velocity_inference_n_pcs(exported, 40L), 40L)
  expect_error(Shennong:::.sn_velocity_inference_n_pcs(exported, 1L), "at least 2")
  expect_error(
    Shennong:::.sn_velocity_inference_n_pcs(
      list(features = c("G1", "G2"), cells = paste0("C", 1:5)),
      30L
    ),
    "at least three exported features"
  )
})

test_that("fate adapter validates probability contracts", {
  object <- make_velocity_test_object()
  valid <- expand.grid(cell = colnames(object), state = c("A", "B"), stringsAsFactors = FALSE)
  valid$probability <- rep(c(0.7, 0.3), each = ncol(object))

  outside <- valid
  outside$probability[[1]] <- 1.1
  expect_error(
    Shennong:::.sn_standardize_fate(list(probabilities = outside), object),
    "between 0 and 1"
  )
  duplicated <- rbind(valid, valid[1, , drop = FALSE])
  expect_error(
    Shennong:::.sn_standardize_fate(list(probabilities = duplicated), object),
    "unique row"
  )
  unnormalized <- valid
  unnormalized$probability[unnormalized$state == "B"] <- 0.2
  expect_error(
    Shennong:::.sn_standardize_fate(list(probabilities = unnormalized), object),
    "sum to 1"
  )
})

test_that("velocity adapter rejects ambiguous cells and invalid transitions", {
  object <- make_velocity_test_object()
  embedding <- Shennong:::.sn_velocity_embedding(object, "umap", 1:2)
  valid <- velocity_adapter_result(object)

  duplicate_cells <- valid
  duplicate_cells$cells$cell[[2]] <- duplicate_cells$cells$cell[[1]]
  expect_error(
    Shennong:::.sn_standardize_velocity(duplicate_cells, object, embedding, "scvelo"),
    "unique row"
  )

  unknown_endpoint <- valid
  unknown_endpoint$graph$target[[1]] <- "not-a-cell"
  expect_error(
    Shennong:::.sn_standardize_velocity(unknown_endpoint, object, embedding, "scvelo"),
    "endpoint"
  )

  invalid_weight <- valid
  invalid_weight$graph$probability[[1]] <- -0.1
  expect_error(
    Shennong:::.sn_standardize_velocity(invalid_weight, object, embedding, "scvelo"),
    "non-negative"
  )

  no_vectors <- valid
  no_vectors$cells$velocity_x <- NA_real_
  no_vectors$cells$velocity_y <- NA_real_
  expect_error(
    Shennong:::.sn_standardize_velocity(no_vectors, object, embedding, "scvelo"),
    "no finite"
  )
})

test_that("fate terminal states agree with probability states", {
  object <- make_velocity_test_object()
  probabilities <- expand.grid(
    cell = colnames(object), state = c("A", "B"), stringsAsFactors = FALSE
  )
  probabilities$probability <- 0.5
  expect_error(
    Shennong:::.sn_standardize_fate(
      list(
        probabilities = probabilities,
        terminal_states = tibble::tibble(cell = "cell1", state = "C")
      ),
      object
    ),
    "absent from probability states"
  )
})

test_that("temporary CellRank run directories are cleaned after import", {
  object <- make_velocity_test_object()
  h5ad <- tempfile(fileext = ".h5ad")
  file.create(h5ad)
  root <- withr::local_tempdir(pattern = "fate-parent-")
  sentinel <- file.path(root, "user-owned.txt")
  writeLines("preserve", sentinel)
  owned_run_dir <- NULL
  velocity_result <- list(models = list(artifacts = list(output_h5ad = h5ad)))

  output <- testthat::with_mocked_bindings(
    Shennong:::.sn_run_fate_pixi(
      velocity_result,
      list(run_dir = root, keep_run_dir = FALSE)
    ),
    sn_call_pixi_environment = function(environment, command, args, ...) {
      output_dir <- args[[match("--output-dir", args) + 1L]]
      owned_run_dir <<- dirname(output_dir)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(
        data.frame(cell = colnames(object), state = "A", probability = 1),
        file.path(output_dir, "fate_probabilities.csv"), row.names = FALSE
      )
      utils::write.csv(
        data.frame(cell = "cell1", state = "A", probability = 1),
        file.path(output_dir, "terminal_states.csv"), row.names = FALSE
      )
      jsonlite::write_json(
        list(mode = "fate"), file.path(output_dir, "manifest.json"), auto_unbox = TRUE
      )
      invisible(character())
    },
    .package = "Shennong"
  )

  expect_true(dir.exists(root))
  expect_true(file.exists(sentinel))
  expect_false(dir.exists(owned_run_dir))
  expect_false(output$artifacts$run_dir_retained)
})

test_that("fate metadata ownership is stable across reruns", {
  object <- make_velocity_test_object()
  probabilities <- expand.grid(
    cell = colnames(object), state = c("A-B", "A B"), stringsAsFactors = FALSE
  )
  probabilities$probability <- 0.5
  updated <- sn_run_fate(
    object, backend_control = list(result = list(probabilities = probabilities))
  )
  result <- sn_get_result(updated, "fate", "fate")
  mapping <- result$tables$state_metadata

  expect_equal(length(unique(mapping$metadata_column)), 2L)
  expect_true(all(mapping$metadata_column %in% colnames(updated[[]])))
  rerun <- sn_run_fate(
    updated, backend_control = list(result = list(probabilities = probabilities))
  )
  rerun_mapping <- sn_get_result(rerun, "fate", "fate")$tables$state_metadata
  expect_identical(rerun_mapping, mapping)
  expect_equal(
    sum(grepl("^fate_fate_A_B", colnames(rerun[[]]))),
    2L
  )

  modified <- updated
  modified[[mapping$metadata_column[[1L]]]] <- 0
  expect_error(
    sn_run_fate(
      modified, backend_control = list(result = list(probabilities = probabilities))
    ),
    "user-modified"
  )
})

test_that("fate metadata never overwrites an unowned user column", {
  object <- make_velocity_test_object()
  object[["fate_fate_A_B"]] <- seq_len(ncol(object))
  probabilities <- expand.grid(
    cell = colnames(object), state = c("A-B", "A B"), stringsAsFactors = FALSE
  )
  probabilities$probability <- 0.5

  expect_error(
    sn_run_fate(
      object, backend_control = list(result = list(probabilities = probabilities))
    ),
    "not owned"
  )
})

test_that("RegVelo prior GRNs are normalized to regulator-target edges", {
  path <- tempfile(fileext = ".csv")
  prior <- matrix(
    c(0, 1.5, -0.5, 0),
    nrow = 2,
    dimnames = list(target = c("G1", "G2"), regulator = c("TF1", "TF2"))
  )
  Shennong:::.sn_write_regvelo_prior_grn(prior, path)
  edges <- utils::read.csv(path)

  expect_setequal(names(edges), c("regulator", "target", "weight"))
  expect_equal(nrow(edges), 2L)
  expect_setequal(edges$regulator, c("TF1", "TF2"))
  expect_error(
    Shennong:::.sn_write_regvelo_prior_grn(NULL, path),
    "prior_grn"
  )
})

test_that("velocity and fate plots render", {
  object <- make_velocity_test_object()
  velocity <- sn_run_velocity(object, backend_control = list(result = velocity_adapter_result(object)), return_object = FALSE)
  fate <- sn_run_fate(
    object,
    backend_control = list(result = list(probabilities = tibble::tibble(
      cell = rep(colnames(object), 2), state = rep(c("A", "B"), each = ncol(object)),
      probability = c(
        seq(0, 1, length.out = ncol(object)),
        1 - seq(0, 1, length.out = ncol(object))
      )
    ))), return_object = FALSE
  )
  for (plot in list(sn_plot_velocity(velocity), sn_plot_fate(fate))) {
    expect_s3_class(plot, "ggplot")
    expect_silent(ggplot2::ggplotGrob(plot))
  }
})

test_that("trajectory pixi environment and dynamics methods are registered", {
  expect_true("trajectory" %in% sn_list_pixi_environments())
  expect_true(sn_get_method_status("scvelo", "velocity")$implemented)
  expect_true(sn_get_method_status("regvelo", "velocity")$implemented)
  expect_true(sn_get_method_status("cellrank", "fate")$implemented)
  manifest <- readLines(sn_get_pixi_config_path("trajectory"), warn = FALSE)
  expect_true(any(grepl("regvelo", manifest, fixed = TRUE)))
})
