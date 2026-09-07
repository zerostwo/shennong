library(testthat)

make_grn_test_object <- function() {
  skip_if_not_installed("SeuratObject")
  set.seed(11)
  counts <- matrix(rpois(8 * 24, lambda = 2), nrow = 8)
  rownames(counts) <- paste0("G", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  counts[c("G1", "G2", "G3"), 1:12] <- counts[c("G1", "G2", "G3"), 1:12] + 8
  counts[c("G4", "G5", "G6"), 13:24] <- counts[c("G4", "G5", "G6"), 13:24] + 8
  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE), project = "grn-test"
  )
  object$sample <- rep(paste0("S", 1:4), each = 6)
  object$condition <- rep(c("control", "control", "treated", "treated"), each = 6)
  object$cell_type <- rep(c("type_a", "type_b"), 12)
  Seurat::NormalizeData(object, verbose = FALSE)
}

grn_adapter_output <- function(object) {
  list(edges = tibble::tibble(
    regulator = rep(c("G1", "G4"), each = 3),
    target_gene = c("G2", "G3", "G7", "G5", "G6", "G8"),
    importance = c(0.9, 0.8, 0.2, 0.95, 0.75, 0.1)
  ))
}

test_that("GRN adapters derive regulons, activity, and group specificity", {
  object <- make_grn_test_object()
  updated <- sn_run_grn(
    object, method = "genie3", result_id = "test_grn", group_by = "condition",
    backend_control = list(result = grn_adapter_output(object), top_targets = 2)
  )
  result <- sn_get_result(updated, "grn", "test_grn")

  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_equal(nrow(result$tables$edges), 6L)
  expect_equal(nrow(result$tables$regulons), 4L)
  expect_equal(nrow(result$tables$activity), 2L * ncol(object))
  expect_equal(unique(result$tables$specificity$metric), "eta_squared")
  expect_true(all(result$tables$specificity$specificity_score >= 0))
  expect_true(all(result$tables$specificity$specificity_score <= 1))
  expect_true(all(c("test_grn_G1", "test_grn_G4") %in% colnames(updated[[]])))
})

test_that("external GRN activity is standardized and retained", {
  object <- make_grn_test_object()
  output <- grn_adapter_output(object)
  output$activity <- transform(
    expand.grid(cell_id = colnames(object), tf = c("G1", "G4"), stringsAsFactors = FALSE),
    auc = seq_len(2L * ncol(object)) / 100
  )
  for (method in c("pyscenic", "scenic", "grnboost2")) {
    result <- sn_run_grn(
      object, method = method, group_by = "cell_type",
      backend_control = list(result = output), return_object = FALSE
    )
    expect_equal(result$method, method)
    expect_equal(nrow(result$tables$activity), 2L * ncol(object))
    expect_true(all(result$tables$activity$method == method))
  }
})

test_that("fallback GRN activity preserves signed edge weights", {
  expression <- Matrix::Matrix(
    matrix(
      c(10, 0, 0, 10),
      nrow = 2,
      dimnames = list(c("G2", "G3"), c("cell1", "cell2"))
    ),
    sparse = TRUE
  )
  regulons <- tibble::tibble(
    regulon = "G1", regulator = "G1", target = c("G2", "G3"),
    weight = c(1, -1), rank = 1:2
  )

  scores <- Shennong:::.sn_score_grn_weighted_mean(expression, regulons)

  expect_equal(as.numeric(scores["G1", ]), c(5, -5))
})

test_that("GRN backend seeds do not change caller RNG state", {
  set.seed(102)
  before <- .Random.seed
  first <- Shennong:::.sn_with_seed(7, runif(4))
  second <- Shennong:::.sn_with_seed(7, runif(4))
  expect_equal(first, second)
  expect_identical(.Random.seed, before)
})

test_that("signed GRNs fail clearly for every unweighted scoring fallback", {
  object <- make_grn_test_object()
  output <- list(edges = tibble::tibble(
    regulator = "G1", target_gene = c("G2", "G3"), importance = c(1, -1)
  ))
  for (activity_method in c("mean", "ucell")) {
    expect_error(
      sn_run_grn(
        object, backend_control = list(result = output, activity_method = activity_method),
        return_object = FALSE
      ),
      "Signed GRN edges require.*weighted_mean"
    )
  }
})

test_that("GRN metadata names are unique and retain a regulon mapping", {
  object <- make_grn_test_object()
  object$collision_A_B <- seq_len(ncol(object))
  output <- list(
    edges = tibble::tibble(
      regulator = c("A-B", "A B"), target_gene = c("G2", "G3"),
      importance = c(1, 1)
    ),
    activity = transform(
      expand.grid(cell_id = colnames(object), tf = c("A-B", "A B"), stringsAsFactors = FALSE),
      auc = seq_len(2L * ncol(object)) / 100
    )
  )
  updated <- sn_run_grn(
    object, method = "pyscenic", result_id = "collision",
    backend_control = list(result = output)
  )
  result <- sn_get_result(updated, "grn", "collision")
  mapping <- result$tables$regulon_metadata_columns

  expect_setequal(mapping$regulon, c("A-B", "A B"))
  expect_false(anyDuplicated(mapping$metadata_column) > 0L)
  expect_false("collision_A_B" %in% mapping$metadata_column)
  expect_true(all(mapping$metadata_column %in% colnames(updated[[]])))
})

test_that("GRN runners receive explicit backend context", {
  object <- make_grn_test_object()
  runner <- function(object, method, assay, layer, regulators, group_by, backend_control) {
    expect_equal(method, "pyscenic")
    expect_equal(regulators, c("G1", "G4"))
    grn_adapter_output(object)
  }
  result <- sn_run_grn(
    object, method = "pyscenic", regulators = c("G1", "G4"),
    backend_control = list(runner = runner), return_object = FALSE
  )
  expect_equal(result$diagnostics$regulons, 2L)
})

test_that("GRN edge contracts reject malformed adapter output", {
  object <- make_grn_test_object()
  expect_error(
    sn_run_grn(object, backend_control = list(result = list(edges = data.frame(a = 1, b = 2)))),
    "require source/regulator"
  )
})

test_that("regulon network, activity, and specificity plots render", {
  object <- make_grn_test_object()
  result <- sn_run_grn(
    object, group_by = "condition", backend_control = list(result = grn_adapter_output(object)),
    return_object = FALSE
  )
  for (type in c("network", "activity", "specificity")) {
    plot <- sn_plot_regulon(result, type = type)
    expect_s3_class(plot, "ggplot")
    expect_silent(ggplot2::ggplotGrob(plot))
  }
})

test_that("GRN methods are registered", {
  methods <- sn_list_methods("grn")
  expect_true(all(c("genie3", "pyscenic", "scenic", "grnboost2") %in% methods$name))
  expect_true(methods$default[methods$name == "genie3"])
})
