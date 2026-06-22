library(testthat)

make_registry_test_object <- function() {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 1, 0, 3, 2, 1), nrow = 2),
    sparse = TRUE
  )
  rownames(counts) <- c("gene1", "gene2")
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  SeuratObject::CreateSeuratObject(counts = counts, project = "registry-test")
}

test_that("stored result registry describes listable and direct misc collections", {
  registry <- Shennong:::.sn_misc_result_registry()

  expect_s3_class(registry, "data.frame")
  expect_true(all(c(
    "collection", "type", "schema_version", "required_fields",
    "listable", "table_required", "reader", "writer"
  ) %in% colnames(registry)))
  expect_true(any(registry$collection == "de_results" & registry$listable))
  expect_true(any(registry$collection == "enrichment_results" & registry$listable))
  expect_true(any(registry$collection == "interpretation_results" & registry$listable))
  expect_true(any(registry$collection == "deconvolution_results" & registry$listable))
  expect_true(any(registry$collection == "sn_run_cluster" & !registry$listable))
  expect_true(any(registry$collection == "bpcells_layers" & !registry$listable))
  expect_equal(anyDuplicated(registry$collection), 0L)
})

test_that("misc result storage validates known collection schema", {
  object <- make_registry_test_object()
  good <- list(
    schema_version = "1.0.0",
    package_version = "0.0.0",
    created_at = "2026-06-22 00:00:00 UTC",
    table = tibble::tibble(gene = "gene1", score = 1),
    analysis = "markers",
    method = "seurat"
  )

  object <- Shennong:::.sn_store_misc_result(
    object = object,
    collection = "de_results",
    store_name = "markers",
    result = good
  )
  expect_equal(object@misc$de_results$markers$table$gene, "gene1")

  expect_error(
    Shennong:::.sn_store_misc_result(
      object = object,
      collection = "unknown_results",
      store_name = "x",
      result = good
    ),
    "Unknown Shennong misc collection"
  )

  expect_error(
    Shennong:::.sn_store_misc_result(
      object = object,
      collection = "de_results",
      store_name = "bad",
      result = list(schema_version = "1.0.0", table = tibble::tibble())
    ),
    "missing required field"
  )
})

test_that("misc result retrieval rejects malformed stored entries", {
  object <- make_registry_test_object()
  object@misc$de_results <- list(
    legacy = list(
      schema_version = "1.0.0",
      package_version = "0.0.0",
      created_at = "2026-06-22 00:00:00 UTC",
      analysis = "markers"
    )
  )

  expect_error(
    Shennong:::.sn_get_misc_result(object, "de_results", "legacy"),
    "Stored result 'legacy' in collection 'de_results' does not match"
  )
})

test_that("sn_run_llm filters provider arguments and normalizes provider outputs", {
  captured <- NULL
  provider_without_dots <- function(messages, model = NULL) {
    captured <<- list(messages = messages, model = model)
    "plain response"
  }

  plain <- sn_run_llm(
    messages = list(list(role = "user", content = "hello")),
    provider = provider_without_dots,
    model = "demo-model",
    structured_type = list(name = "ignored"),
    tools = list("ignored"),
    extra = "ignored"
  )

  expect_equal(plain$text, "plain response")
  expect_equal(plain$model, "demo-model")
  expect_equal(names(captured), c("messages", "model"))

  structured <- sn_run_llm(
    messages = list(list(role = "user", content = "hello")),
    provider = function(messages, structured_type = NULL, ...) {
      list(structured = list(answer = "ok", n = 1))
    },
    structured_type = list(name = "demo")
  )
  expect_equal(structured$structured$answer, "ok")
  expect_match(structured$text, '"answer":"ok"', fixed = TRUE)

  expect_error(
    sn_run_llm(
      messages = list(list(role = "user", content = "hello")),
      provider = function(...) list(unexpected = TRUE)
    ),
    "provider.*single string.*text.*structured"
  )
})
