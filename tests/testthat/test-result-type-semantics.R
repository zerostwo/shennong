library(testthat)

make_type_semantic_result <- function(analysis_type, primary) {
  Shennong:::.sn_new_analysis_result(
    analysis_type = analysis_type,
    result_id = paste0(analysis_type, "_semantic_test"),
    method = "test",
    tables = list(primary = primary)
  )
}

expect_type_semantic_error <- function(result, pattern) {
  report <- sn_validate_result(result, error = FALSE)
  expect_false(report$valid)
  expect_match(paste(report$errors, collapse = "\n"), pattern)
  invisible(report)
}

test_that("every built-in analysis result type has a semantic registry entry", {
  expect_setequal(
    names(Shennong:::.sn_analysis_result_type_specs()),
    Shennong:::.sn_analysis_result_builtin_types()
  )
})

test_that("previously generic built-in result types reject unstructured primary tables", {
  types <- c(
    "deconvolution", "enrichment", "differential_abundance", "cnv",
    "cell_communication", "bulk_de", "spatial_domains"
  )
  for (type in types) {
    result <- Shennong:::.sn_new_analysis_result(
      analysis_type = "custom_extension",
      result_id = paste0(type, "_malformed"),
      method = "test",
      tables = list(primary = data.frame(foo = 1))
    )
    result$analysis_type <- type
    result$provenance$analysis_type <- type
    expect_type_semantic_error(result, "tables\\$primary|primary.*missing|must contain")
  }
})

test_that("trajectory semantics allow coherent unassigned cells but reject invalid values", {
  result <- make_type_semantic_result(
    "trajectory",
    data.frame(
      cell = c("cell1", "cell2"),
      primary_lineage = c("Lineage1", "Lineage1"),
      primary_pseudotime = c(0, 1),
      pseudotime_Lineage1 = c(0, 1),
      weight_Lineage1 = c(1, 0.5)
    )
  )

  duplicate_cell <- result
  duplicate_cell$tables$primary$cell[[2]] <- "cell1"
  expect_type_semantic_error(duplicate_cell, "unique rows.*`cell`")

  non_numeric <- result
  non_numeric$tables$primary$primary_pseudotime <- c("early", "late")
  expect_type_semantic_error(non_numeric, "`primary_pseudotime` must be numeric")

  for (value in c(Inf, -Inf, NaN)) {
    non_finite <- result
    non_finite$tables$primary$primary_pseudotime[[2]] <- value
    expect_type_semantic_error(non_finite, "`primary_pseudotime` may contain missing values")
  }

  unassigned <- result
  unassigned$tables$primary$primary_lineage[[2]] <- NA_character_
  unassigned$tables$primary$primary_pseudotime[[2]] <- NA_real_
  unassigned$tables$primary$pseudotime_Lineage1[[2]] <- NA_real_
  expect_true(sn_validate_result(unassigned)$valid)

  mismatched_missing <- unassigned
  mismatched_missing$tables$primary$primary_pseudotime[[2]] <- 0.5
  expect_type_semantic_error(mismatched_missing, "must be missing for the same rows")

  invalid_weight <- result
  invalid_weight$tables$primary$weight_Lineage1[[2]] <- 1.2
  expect_type_semantic_error(invalid_weight, "`weight_Lineage1` must lie in \\[0, 1\\]")
})

test_that("annotation primary tables require unique cells and present predictions", {
  result <- make_type_semantic_result(
    "annotation",
    data.frame(cell = c("cell1", "cell2"), prediction = c("T", "B"))
  )

  duplicate_cell <- result
  duplicate_cell$tables$primary$cell[[2]] <- "cell1"
  expect_type_semantic_error(duplicate_cell, "unique rows.*`cell`")

  missing_prediction <- result
  missing_prediction$tables$primary$prediction[[2]] <- NA_character_
  expect_type_semantic_error(missing_prediction, "`prediction` must not contain missing")

  non_character <- result
  non_character$tables$primary$prediction <- factor(non_character$tables$primary$prediction)
  expect_type_semantic_error(non_character, "`prediction` must be character")
})

test_that("DE primary tables require explicit character gene identifiers", {
  result <- make_type_semantic_result(
    "de",
    data.frame(gene = c("gene1", "gene2"))
  )

  for (value in c(NA_character_, "", "  ")) {
    missing_gene <- result
    missing_gene$tables$primary$gene[[2]] <- value
    expect_type_semantic_error(missing_gene, "`gene` must not contain missing or blank")
  }

  non_character <- result
  non_character$tables$primary$gene <- seq_len(nrow(non_character$tables$primary))
  expect_type_semantic_error(non_character, "`gene` must be character")
})

test_that("program scoring uses an entity-program key and numeric scores", {
  result <- make_type_semantic_result(
    "program_scoring",
    data.frame(
      entity = c("cell1", "cell1", "cell2"),
      program = c("program1", "program2", "program1"),
      score = c(0.1, 0.2, 0.3)
    )
  )
  expect_true(sn_validate_result(result, error = FALSE)$valid)

  duplicate_pair <- result
  duplicate_pair$tables$primary[2, c("entity", "program")] <-
    duplicate_pair$tables$primary[1, c("entity", "program")]
  expect_type_semantic_error(duplicate_pair, "unique rows.*`entity`, `program`")

  non_numeric <- result
  non_numeric$tables$primary$score <- as.character(non_numeric$tables$primary$score)
  expect_type_semantic_error(non_numeric, "`score` must be numeric")

  for (value in c(NA_real_, Inf)) {
    non_finite <- result
    non_finite$tables$primary$score[[2]] <- value
    expect_type_semantic_error(non_finite, "`score` must contain only finite")
  }
})

test_that("state priority tables use unique states and finite numeric priorities", {
  result <- make_type_semantic_result(
    "state_priority",
    data.frame(state = c("state1", "state2"), priority_score = c(0.8, 0.6))
  )

  duplicate_state <- result
  duplicate_state$tables$primary$state[[2]] <- "state1"
  expect_type_semantic_error(duplicate_state, "unique rows.*`state`")

  non_numeric <- result
  non_numeric$tables$primary$priority_score <- c("high", "low")
  expect_type_semantic_error(non_numeric, "`priority_score` must be numeric")

  non_finite <- result
  non_finite$tables$primary$priority_score[[2]] <- NA_real_
  expect_type_semantic_error(non_finite, "`priority_score` must contain only finite")
})

test_that("scissor primary tables use a typed cell-selection contract", {
  result <- make_type_semantic_result(
    "scissor",
    data.frame(
      cell = c("cell1", "cell2"),
      coefficient = c(-0.5, 0.5),
      selection = c("Scissor-", "Scissor+")
    )
  )

  duplicate_cell <- result
  duplicate_cell$tables$primary$cell[[2]] <- "cell1"
  expect_type_semantic_error(duplicate_cell, "unique rows.*`cell`")

  missing_cell <- result
  missing_cell$tables$primary$cell[[2]] <- NA_character_
  expect_type_semantic_error(missing_cell, "key column.*`cell`.*must not contain missing")

  non_numeric <- result
  non_numeric$tables$primary$coefficient <- as.character(
    non_numeric$tables$primary$coefficient
  )
  expect_type_semantic_error(non_numeric, "`coefficient` must be numeric")

  non_finite <- result
  non_finite$tables$primary$coefficient[[2]] <- Inf
  expect_type_semantic_error(non_finite, "`coefficient` must contain only finite")

  missing_selection <- result
  missing_selection$tables$primary$selection[[2]] <- NA_character_
  expect_type_semantic_error(missing_selection, "`selection` must not contain missing")

  invalid_selection <- result
  invalid_selection$tables$primary$selection[[2]] <- "positive"
  expect_type_semantic_error(invalid_selection, "`selection` must contain only: Scissor\\+")
})

test_that("bulk survival estimates are finite, positive, and probabilistic", {
  result <- make_type_semantic_result(
    "bulk_survival",
    data.frame(
      feature = c("gene1", "gene2"),
      hazard_ratio = c(1.2, 0.8),
      conf_low = c(1.1, 0.7),
      conf_high = c(1.3, 0.9),
      p_value = c(0, 1)
    )
  )

  for (column in c("hazard_ratio", "conf_low", "conf_high")) {
    non_finite <- result
    non_finite$tables$primary[[column]][[1]] <- Inf
    expect_type_semantic_error(non_finite, paste0("`", column, "` must contain only finite"))

    non_positive <- result
    non_positive$tables$primary[[column]][[1]] <- 0
    expect_type_semantic_error(non_positive, paste0("`", column, "` must contain only positive"))
  }

  non_numeric <- result
  non_numeric$tables$primary$hazard_ratio <- as.character(
    non_numeric$tables$primary$hazard_ratio
  )
  expect_type_semantic_error(non_numeric, "`hazard_ratio` must be numeric")

  non_finite_p <- result
  non_finite_p$tables$primary$p_value[[1]] <- NA_real_
  expect_type_semantic_error(non_finite_p, "`p_value` must contain only finite")

  for (value in c(-0.01, 1.01)) {
    out_of_range <- result
    out_of_range$tables$primary$p_value[[1]] <- value
    expect_type_semantic_error(out_of_range, "`p_value` must lie in \\[0, 1\\]")
  }

  duplicate_feature <- result
  duplicate_feature$tables$primary$feature[[2]] <- "gene1"
  expect_type_semantic_error(duplicate_feature, "unique rows.*`feature`")

  non_character_feature <- result
  non_character_feature$tables$primary$feature <- c(1, 2)
  expect_type_semantic_error(non_character_feature, "`feature` must be character")
})

test_that("cross-row and cross-column scientific invariants fail closed", {
  fate <- make_type_semantic_result(
    "fate",
    data.frame(
      cell = rep(c("cell1", "cell2"), each = 2L),
      state = rep(c("A", "B"), 2L),
      probability = c(0.7, 0.3, 0.4, 0.6)
    )
  )
  invalid_fate <- fate
  invalid_fate$tables$primary$probability[[4L]] <- 0.5
  expect_type_semantic_error(invalid_fate, "probabilities must sum to 1")

  spatial <- make_type_semantic_result(
    "spatial_integration",
    data.frame(cell = c("cell1", "cell2"), dim1 = c(0, 1), dim2 = c(1, 0))
  )
  one_dimension <- spatial
  one_dimension$tables$primary$dim2 <- NULL
  expect_type_semantic_error(one_dimension, "at least two numeric embedding dimensions")

  survival <- make_type_semantic_result(
    "bulk_survival",
    data.frame(
      feature = "gene1", hazard_ratio = 1.2, conf_low = 0.9,
      conf_high = 1.5, p_value = 0.05
    )
  )
  reversed_interval <- survival
  reversed_interval$tables$primary$conf_low <- 1.3
  expect_type_semantic_error(reversed_interval, "conf_low <= hazard_ratio <= conf_high")
})
