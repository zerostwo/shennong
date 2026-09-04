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

test_that("misc registry contains artifacts only", {
  registry <- Shennong:::.sn_misc_result_registry()

  expect_s3_class(registry, "data.frame")
  expect_true(all(c(
    "collection", "type", "schema_version", "required_fields",
    "contract_scope", "listable", "table_required", "reader", "writer"
  ) %in% colnames(registry)))
  expect_true(any(registry$collection == "sn_run_cluster" & !registry$listable))
  expect_true(any(registry$collection == "bpcells_layers" & !registry$listable))
  expect_true(any(registry$collection == "label_transfer" & registry$contract_scope == "artifact"))
  expect_true(any(registry$collection == "scdesign3" & registry$contract_scope == "artifact"))
  expect_true(all(registry$contract_scope == "artifact"))
  expect_equal(anyDuplicated(registry$collection), 0L)
})

test_that("canonical result storage validates the v2 contract", {
  object <- make_registry_test_object()
  good <- list(
    tables = list(primary = tibble::tibble(gene = "gene1", score = 1)),
    analysis = "markers",
    method = "seurat",
    backend = "seurat"
  )

  object <- sn_store_result(
    object = object,
    type = "de",
    result_id = "markers",
    result = good
  )
  expect_equal(sn_get_result(object, "de", "markers")$tables$primary$gene, "gene1")
  expect_true(is.null(object@misc$de_results))

  expect_error(
    sn_store_result(object, "de", "", good),
    "result_id"
  )
})

test_that("canonical retrieval rejects malformed stored entries", {
  object <- make_registry_test_object()
  object@misc$shennong <- list(results = list(de = list(
    malformed = list(schema_version = "2.0.0", result_id = "malformed")
  )))

  expect_error(
    sn_get_result(object, "de", "malformed"),
    "Invalid Shennong analysis result"
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
