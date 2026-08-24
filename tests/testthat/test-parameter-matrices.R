.parameter_default_text <- function(value) {
  if (identical(value, quote(expr = ))) return("<required>")
  paste(deparse(value, width.cutoff = 500L), collapse = "")
}

.parameter_conformance_path <- function(...) {
  local <- testthat::test_path("..", "..", "inst", "conformance", ...)
  if (file.exists(local)) return(local)
  installed <- system.file("conformance", ..., package = "Shennong")
  if (!nzchar(installed)) {
    stop("Installed Shennong conformance artifact is missing.", call. = FALSE)
  }
  installed
}

test_that("public API parameter inventory is complete and fresh", {
  inventory <- jsonlite::read_json(
    .parameter_conformance_path("public-api-parameters.json"),
    simplifyVector = FALSE
  )
  exports <- sort(
    grep("^sn_", getNamespaceExports("Shennong"), value = TRUE),
    method = "radix"
  )
  recorded <- vapply(inventory$functions, `[[`, character(1), "function_name")

  expect_identical(inventory$schema_version, "shennong.public-api-parameters/v1")
  expect_identical(
    inventory$package_version,
    as.character(utils::packageVersion("Shennong"))
  )
  expect_equal(inventory$function_count, length(exports))
  expect_identical(recorded, exports)
  expect_identical(recorded, sort(unique(recorded), method = "radix"))
  expect_true(all(vapply(
    inventory$functions,
    `[[`,
    logical(1),
    "runtime_classified"
  )))

  for (subject in inventory$functions) {
    fun <- get(subject$function_name, envir = asNamespace("Shennong"), inherits = FALSE)
    defaults <- formals(fun)
    digest <- digest::digest(
      vapply(defaults, .parameter_default_text, character(1)),
      algo = "sha256",
      serialize = TRUE
    )
    expect_identical(subject$formal_digest, digest, info = subject$function_name)
    expect_identical(
      vapply(subject$parameters, `[[`, character(1), "name"),
      if (is.null(names(defaults))) character() else names(defaults),
      info = subject$function_name
    )
    selector_rows <- Filter(
      function(parameter) identical(parameter$class, "selector"),
      subject$parameters
    )
    for (parameter in selector_rows) {
      expect_equal(
        parameter$minimum_cases,
        length(parameter$selector_values),
        info = paste(subject$function_name, parameter$name)
      )
    }
  }
})

test_that("clustering selector matrix covers every current backend", {
  matrix <- jsonlite::read_json(
    .parameter_conformance_path("matrices", "clustering-methods.json"),
    simplifyVector = FALSE
  )
  methods <- vapply(matrix$integration_methods, `[[`, character(1), "method")
  expected_normalization <- eval(formals(Shennong::sn_run_cluster)$normalization_method)

  expect_identical(methods, Shennong:::.sn_supported_integration_methods())
  expect_identical(unlist(matrix$normalization_methods), expected_normalization)
  expect_true(all(vapply(
    matrix$integration_methods,
    function(subject) {
      nzchar(subject$contract_kind) && nzchar(subject$status) &&
        length(subject$required_cases) > 0L
    },
    logical(1)
  )))
  expect_true(all(c(
    "normalization_method", "integration_method", "batch"
  ) %in% unlist(matrix$pairwise_axes)))
  expect_true(all(c(
    "split_layers", "bpcells", "checkpoint_resume", "multi_method_grid"
  ) %in% unlist(matrix$high_risk_cases)))
})

test_that("enrichment matrix inventories dispatch without overstating evidence", {
  matrix <- jsonlite::read_json(
    .parameter_conformance_path("matrices", "enrichment-methods.json"),
    simplifyVector = FALSE
  )
  analyses <- unlist(matrix$analysis_modes)
  databases <- unlist(matrix$database_labels)
  input_modes <- unlist(matrix$input_modes)
  cases <- matrix$cases
  keys <- vapply(cases, function(case) {
    paste(case$analysis, case$database, case$input_mode, sep = "\r")
  }, character(1))
  expected <- as.vector(outer(
    as.vector(outer(analyses, databases, paste, sep = "\r")),
    input_modes,
    paste,
    sep = "\r"
  ))

  expect_identical(matrix$schema_version, "shennong.parameter-matrix/v1")
  expect_identical(matrix[["function"]], "sn_enrich")
  expect_setequal(keys, expected)
  expect_identical(anyDuplicated(keys), 0L)
  expect_setequal(databases, c("GO", "GOBP", "GOMF", "GOCC", "KEGG", "MSIGDB_COLLECTION"))
  expect_true(all(vapply(cases, function(case) {
    case$status %in% c("pilot", "pending", "unsupported")
  }, logical(1))))

  grouped_gsea <- Filter(function(case) {
    identical(case$analysis, "gsea") && identical(case$input_mode, "grouped")
  }, cases)
  expect_length(grouped_gsea, length(databases))
  expect_true(all(vapply(grouped_gsea, function(case) {
    identical(case$status, "unsupported") && nzchar(case$reason)
  }, logical(1))))

  pilots <- Filter(function(case) identical(case$status, "pilot"), cases)
  expect_setequal(
    vapply(pilots, `[[`, character(1), "contract_id"),
    unlist(matrix$pilot_contracts)
  )
  contracts <- .conformance_read_contracts()
  for (id in unlist(matrix$pilot_contracts)) {
    expect_true(id %in% names(contracts), info = id)
    expect_identical(contracts[[id]]$status, "pilot", info = id)
    expect_identical(
      contracts[[id]]$shennong[["function"]],
      "sn_enrich",
      info = id
    )
  }

  axes <- vapply(matrix$parameter_axes, `[[`, character(1), "axis")
  expect_setequal(
    axes,
    c("universe", "cutoffs", "gene_set_size", "duplicate_gene_method", "seed", "backend_version")
  )
  expect_false(any(vapply(matrix$parameter_axes, function(axis) {
    length(axis$cases) == 0L
  }, logical(1))))
})
