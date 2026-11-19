bundle_consumer_reference <- function(value, path) {
  dimensions <- dim(value)
  axes <- dimnames(value)
  coordinates <- arrayInd(seq_along(value), dimensions)
  list(
    path = path,
    storage_type = typeof(value),
    dim = I(dimensions),
    dimnames = if (is.null(axes)) NULL else unname(lapply(axes, function(x) if (is.null(x)) NULL else I(x))),
    axis_names = if (is.null(names(axes))) NULL else I(names(axes)),
    coordinates = lapply(seq_along(value), function(index) {
      positions <- as.integer(coordinates[index, ])
      list(index = I(positions), value = unname(as.vector(value)[index]),
        labels = lapply(seq_along(dimensions), function(axis) {
          if (is.null(axes[[axis]])) NULL else unname(axes[[axis]][positions[[axis]]])
        }))
    })
  )
}

make_bundle_consumer_result <- function() {
  counts <- matrix(c(10, 20, 30, 8, 11, 19, 27, 12, 8, 25, 24, 14, 14, 16, 31, 9), 4,
    dimnames = list(c("g/1", "g~2", "g3", "g4"), c("sample/B", "sample~A", "sample C", "sampleD")))
  sn_assess_bulk_qc(counts, top_variable = 4L)
}

test_that("an independent Python consumer restores bundle array coordinates", {
  python <- Sys.which("python3")
  skip_if(!nzchar(python), "python3 is needed for the independent JSON consumer")
  result <- make_bundle_consumer_result()
  expected <- list(
    bundle_consumer_reference(result$embeddings$pca, list("result", "embeddings", "pca")),
    bundle_consumer_reference(result$models$pca$rotation, list("result", "models", "pca", "rotation")),
    bundle_consumer_reference(result$models$pca$x, list("result", "models", "pca", "x"))
  )
  examples <- list(
    "single/~" = matrix(7L, 1L, 1L, dimnames = list(sample = "S/~", component = "PC1")),
    empty_rows = matrix(numeric(), 0L, 3L, dimnames = list(NULL, c("a", "b", "c"))),
    empty_columns = matrix(numeric(), 2L, 0L, dimnames = list(c("a", "b"), NULL)),
    "cube/~" = array(seq_len(12L), c(2L, 3L, 2L),
      dimnames = list(sample = c("s2", "s1"), component = c("c/~", "c2", "c3"), draw = c("d1", "d2"))),
    empty_cube = array(numeric(), c(2L, 0L, 3L)),
    logical_values = matrix(c(TRUE, NA, FALSE, TRUE), 2L),
    characters = matrix(c("x/y", "~"), 1L)
  )
  result$models$examples <- examples
  for (name in names(examples)) expected[[length(expected) + 1L]] <-
    bundle_consumer_reference(examples[[name]], list("result", "models", "examples", name))

  result$models$mixed <- list(list(matrix(11:14, 2L)), named = matrix(9L, 1L, 1L))
  expected[[length(expected) + 1L]] <- bundle_consumer_reference(result$models$mixed[[1]][[1]], list("result", "models", "mixed", 1L, 1L))
  expected[[length(expected) + 1L]] <- bundle_consumer_reference(result$models$mixed$named, list("result", "models", "mixed", "named"))

  frame <- data.frame(sample = c("s2", "s1"), check.names = FALSE)
  frame[["coords/~"]] <- I(matrix(21:24, 2L, dimnames = list(c("s2", "s1"), c("PC2", "PC1"))))
  frame[["nested/~"]] <- I(list(matrix(31L, 1L, 1L), matrix(41:44, 2L)))
  result$tables[["frame/~"]] <- frame
  expected[[length(expected) + 1L]] <- bundle_consumer_reference(frame[["coords/~"]], list("result", "tables", "frame/~", "coords/~"))
  for (index in seq_len(2L)) expected[[length(expected) + 1L]] <- bundle_consumer_reference(
    frame[["nested/~"]][[index]], list("result", "tables", "frame/~", "nested/~", index))

  directory <- tempfile("bundle-consumer-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  exported <- file.path(directory, "bundle.json")
  expectation <- file.path(directory, "expected.json")
  sn_export_result_bundle(sn_build_result_bundle(result), exported)
  jsonlite::write_json(list(arrays = expected), expectation, auto_unbox = TRUE, na = "null", null = "null", digits = NA)
  consumer <- testthat::test_path("..", "fixtures", "result-bundle-array-consumer.py")
  consume <- function(path) suppressWarnings(system2(python,
    c(shQuote(consumer), shQuote(path), shQuote(expectation)), stdout = TRUE, stderr = TRUE))
  output <- consume(exported)
  expect_identical(attr(output, "status") %||% 0L, 0L, info = paste(output, collapse = "\n"))
  expect_match(paste(output, collapse = "\n"), "Validated 15 arrays and [0-9]+ labeled scalar coordinates")

  # The independent consumer must reject incorrect numbers and lost identities.
  original <- jsonlite::read_json(exported, simplifyVector = FALSE)
  corrupted_path <- file.path(directory, "corrupted.json")
  corrupted <- original
  corrupted$result$embeddings$pca[[1L]][[1L]] <- corrupted$result$embeddings$pca[[1L]][[1L]] + 1
  jsonlite::write_json(corrupted, corrupted_path, auto_unbox = TRUE, null = "null", digits = NA)
  expect_gt(attr(consume(corrupted_path), "status") %||% 0L, 0L)
  corrupted <- original
  corrupted$serialization$arrays[[1L]]$dimnames[[1L]][[1L]] <- "wrong sample"
  jsonlite::write_json(corrupted, corrupted_path, auto_unbox = TRUE, null = "null", digits = NA)
  expect_gt(attr(consume(corrupted_path), "status") %||% 0L, 0L)
})

test_that("bundle export rejects generated JSON key collisions", {
  result <- make_bundle_consumer_result()
  result$models$ambiguous <- structure(list(matrix(1L, 1L), matrix(2L, 1L)), names = c("", "1"))
  destination <- tempfile(fileext = ".json")
  on.exit(unlink(destination), add = TRUE)
  expect_error(sn_export_result_bundle(sn_build_result_bundle(result), destination), "ambiguous|unique|duplicate|colli")
  expect_false(file.exists(destination))
})

test_that("bundle export rejects list arrays whose nested identities are not supported", {
  result <- make_bundle_consumer_result()
  nested <- array(list(matrix(8L, 1L, 1L,
    dimnames = list("sample-X", "PC-X")), 9L), dim = 2L)
  result$models$nested_array <- nested
  expect_error(sn_build_result_bundle(result), "list.*array|array.*list|atomic")

  # Validate the export boundary too when a caller mutates an already-built bundle.
  bundle <- sn_build_result_bundle(make_bundle_consumer_result())
  bundle$result$models$nested_array <- nested
  destination <- tempfile(fileext = ".json")
  on.exit(unlink(destination), add = TRUE)
  expect_error(sn_export_result_bundle(bundle, destination), "list.*array|array.*list|atomic")
  expect_false(file.exists(destination))
})

test_that("bundle build and export reject unnamed data-frame columns", {
  result <- make_bundle_consumer_result()
  frame <- data.frame(sample = c("s1", "s2"))
  frame[[2L]] <- I(matrix(1:4, 2L,
    dimnames = list(c("s1", "s2"), c("PC1", "PC2"))))
  names(frame)[2L] <- ""
  result$tables$unnamed_column <- frame
  expect_error(sn_build_result_bundle(result), "non-empty.*data-frame.*column names")

  # An already-built bundle can be mutated, so export must enforce the same rule.
  bundle <- sn_build_result_bundle(make_bundle_consumer_result())
  bundle$result$tables$unnamed_column <- frame
  destination <- tempfile(fileext = ".json")
  on.exit(unlink(destination), add = TRUE)
  expect_error(sn_export_result_bundle(bundle, destination), "non-empty.*data-frame.*column names")
  expect_false(file.exists(destination))
})
