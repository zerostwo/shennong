library(testthat)

.make_python_contract_object <- function() {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(c(1, 0, 3, 2, 1, 0, 4, 2, 1), nrow = 3),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  SeuratObject::CreateSeuratObject(counts = counts)
}

.python_contract_query <- function(object) {
  list(
    layer = "counts",
    features = rownames(object),
    cells = colnames(object),
    n_features = nrow(object),
    n_cells = ncol(object)
  )
}

.write_backend_manifest <- function(output_dir, manifest) {
  jsonlite::write_json(
    manifest,
    file.path(output_dir, "manifest.json"),
    auto_unbox = TRUE,
    null = "null"
  )
}

.python_contract_pixi_root <- function() {
  candidates <- c(
    file.path("inst", "pixi"),
    file.path("..", "..", "inst", "pixi"),
    system.file("pixi", package = "Shennong")
  )
  available <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (length(available) == 0L) {
    stop("Could not locate the packaged Pixi resources.", call. = FALSE)
  }
  normalizePath(available[[1L]], winslash = "/", mustWork = TRUE)
}

test_that("method registry separates adapter implementation, run readiness, and conformance", {
  cnmf <- sn_get_method_status("cnmf", task = "program_discovery")
  expect_true(cnmf$implemented)
  expect_false(cnmf$runnable)
  expect_false(cnmf$available)
  expect_match(cnmf$reason, "requires an explicit upstream runner or result")

  moran <- sn_get_method_status("morans_i", task = "spatial_svg")
  expect_identical(moran$runtime, "r")
  expect_true(moran$runnable)

  stlearn <- sn_get_method_status("stlearn", task = "spatial_domain")
  expect_identical(stlearn$runtime, "external")
  expect_false(stlearn$implemented)
  expect_false(stlearn$runnable)
  expect_false(stlearn$available)
  expect_match(stlearn$reason, "not implemented yet")

  banksy <- sn_get_method_status("banksy", task = "spatial_domain")
  expect_identical(banksy$package, "Banksy")
  expect_identical(sn_get_method_status("popv", "annotation")$conformance_status, "admitted")
  expect_identical(sn_get_method_status("edger", "bulk_de")$conformance_status, "pilot")
})

test_that("object-level Python wrappers expose cleanup and import budgets", {
  wrappers <- c(
    "sn_run_scarches", "sn_run_scpoli", "sn_run_cellphonedb",
    "sn_run_cell2location", "sn_run_tangram", "sn_run_squidpy",
    "sn_run_spatialdata", "sn_run_stlearn", "sn_run_infercnvpy"
  )
  for (wrapper in wrappers) {
    expect_true(
      all(c("keep_run_dir", "max_artifact_import_gb") %in%
            names(formals(getExportedValue("Shennong", wrapper)))),
      info = wrapper
    )
  }
  expect_identical(formals(sn_run_cellphonedb)$layer, "data")
})

test_that("CellPhoneDB and infercnvpy require declared log-normalized layers", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  object$cell_type <- rep("type", ncol(object))

  expect_error(
    sn_run_cellphonedb(object, layer = "counts", group_by = "cell_type"),
    "requires normalized, log-transformed expression"
  )
  expect_error(
    sn_run_infercnvpy(object, layer = "counts"),
    "requires normalized, log-transformed expression"
  )
  expect_error(
    sn_run_infercnvpy(object),
    "has no `data`/`data\\.\\*` layer"
  )

  normalized <- log1p(SeuratObject::LayerData(object, assay = "RNA", layer = "counts"))
  SeuratObject::LayerData(object, assay = "RNA", layer = "data") <- normalized
  expect_identical(Shennong:::.sn_select_infercnvpy_layer(object, "RNA"), "data")
  expect_invisible(Shennong:::.sn_validate_log_normalized_layer(
    object, "RNA", "data", "test backend"
  ))
})

test_that("pixi manifests use exact direct dependency pins and matching lockfiles", {
  pixi_root <- .python_contract_pixi_root()
  paths <- list.files(pixi_root, pattern = "^pixi\\.toml$", recursive = TRUE, full.names = TRUE)
  expect_gt(length(paths), 0L)
  lock_paths <- list.files(pixi_root, pattern = "^pixi\\.lock$", recursive = TRUE, full.names = TRUE)
  expect_setequal(dirname(lock_paths), dirname(paths))
  pixi <- unname(Sys.which("pixi"))
  pixi_version <- if (nzchar(pixi)) {
    suppressWarnings(system2(pixi, "--version", stdout = TRUE, stderr = TRUE))
  } else {
    character()
  }
  can_verify_with_pixi <- identical(Sys.info()[["sysname"]], "Linux") &&
    identical(pixi_version, "pixi 0.69.0")

  for (path in paths) {
    lines <- readLines(path, warn = FALSE)
    dependency_section <- FALSE
    dependency_lines <- character()
    for (line in lines) {
      section <- regmatches(line, regexec("^\\s*\\[([^]]+)\\]\\s*$", line, perl = TRUE))[[1]]
      if (length(section) > 0L) {
        dependency_section <- grepl("(^|\\.)(pypi-)?dependencies$", section[[2]])
      } else if (dependency_section && grepl("^\\s*[A-Za-z0-9_.-]+\\s*=", line)) {
        dependency_lines <- c(dependency_lines, line)
      }
    }
    floating <- grepl(
      "=\\s*\"(?:\\*|[<>~^!]|=[^=])|version\\s*=\\s*\"(?:\\*|[<>~^!]|=[^=])",
      dependency_lines,
      perl = TRUE
    )
    expect_false(any(floating), info = paste("floating direct dependency in", path))
    versioned <- dependency_lines[!grepl("git\\s*=", dependency_lines)]
    expect_true(
      all(grepl("=\\s*(?:\"==[^\"]+\"|\\{[^}]*version\\s*=\\s*\"==[^\"]+\"[^}]*\\})\\s*$", versioned, perl = TRUE)),
      info = paste("non-exact direct dependency in", path)
    )
    git_lines <- grep("git\\s*=", lines, value = TRUE)
    if (length(git_lines) > 0L) {
      expect_true(all(grepl("rev\\s*=\\s*\"[0-9a-f]{40}\"", git_lines, perl = TRUE)))
    }

    lock_path <- file.path(dirname(path), "pixi.lock")
    lock <- readLines(lock_path, warn = FALSE)
    expect_match(lock[[1]], "^version: [0-9]+$")
    expect_true(any(grepl("^platforms:$", lock)))
    workspace_platforms <- jsonlite::fromJSON(sub(
      "^platforms\\s*=\\s*",
      "",
      grep("^platforms\\s*=", lines, value = TRUE)[[1]]
    ))
    lock_platforms <- sub("^- name: ", "", grep("^- name: ", lock, value = TRUE))
    expect_setequal(lock_platforms, workspace_platforms)
    environment_name <- basename(dirname(path))
    expected_platforms <- if (environment_name %in% c(
      "cell2location", "scib-metrics", "trajectory"
    )) {
      c("linux-64", "osx-64", "osx-arm64")
    } else {
      c("linux-64", "osx-64", "osx-arm64", "win-64")
    }
    expect_setequal(workspace_platforms, expected_platforms)
    gpu_section <- which(lines == "[feature.gpu]")
    if (length(gpu_section) > 0L) {
      gpu_platforms <- jsonlite::fromJSON(sub(
        "^platforms\\s*=\\s*",
        "",
        lines[[gpu_section[[1]] + 1L]]
      ))
      expect_identical(gpu_platforms, "linux-64")
    }
    expect_true(any(grepl("^environments:$", lock)))
    expect_true(any(grepl("^      - (conda|pypi): ", lock)))

    if (can_verify_with_pixi) {
      project_dir <- tempfile(paste0("pixi-lock-check-", basename(dirname(path)), "-"))
      dir.create(project_dir)
      rendered <- gsub("{{ platform }}", "linux-64", lines, fixed = TRUE)
      rendered <- gsub("{{ cuda_major }}", "12", rendered, fixed = TRUE)
      rendered <- gsub("{{ cuda_version }}", "12.6", rendered, fixed = TRUE)
      rendered_path <- file.path(project_dir, "pixi.toml")
      writeLines(rendered, rendered_path, useBytes = TRUE)
      expect_true(file.copy(lock_path, file.path(project_dir, "pixi.lock")))
      check <- suppressWarnings(system2(
        pixi,
        c("lock", "--check", "--no-progress", "--manifest-path", shQuote(rendered_path)),
        stdout = TRUE,
        stderr = TRUE
      ))
      status <- attr(check, "status")
      if (is.null(status)) status <- 0L
      expect_identical(
        as.integer(status),
        0L,
        info = paste("pixi rejected lock for", path, paste(check, collapse = "\n"))
      )
    }
  }
})

test_that("Python run directories are unique and stale runs are rejected", {
  runtime_dir <- tempfile("python-runtime-")
  first <- Shennong:::.sn_default_python_run_dir("squidpy", runtime_dir)
  second <- Shennong:::.sn_default_python_run_dir("squidpy", runtime_dir)
  expect_false(identical(first, second))

  stale <- tempfile("python-stale-")
  dir.create(stale)
  writeLines("old", file.path(stale, "manifest.json"))
  expect_error(
    Shennong:::.sn_prepare_python_run_directory(stale),
    "refusing to mix a new run with stale outputs"
  )
})

test_that("integration Python runs clean successes and sanitize failures by default", {
  success_parent <- tempfile("integration-python-success-parent-")
  dir.create(success_parent)
  writeLines("keep", file.path(success_parent, "user-sentinel.txt"))
  withr::defer(unlink(success_parent, recursive = TRUE, force = TRUE))
  successful_run <- NULL
  value <- Shennong:::.sn_with_integration_python_run(
    method = "scvi",
    integration_control = list(run_dir = success_parent, keep_run_dir = FALSE),
    code = function(control) {
      successful_run <<- control$run_dir
      dir.create(file.path(control$run_dir, "input"), recursive = TRUE)
      writeLines("private", file.path(control$run_dir, "input", "obs.csv"))
      list(value = 1L)
    }
  )
  expect_identical(value$value, 1L)
  expect_false(dir.exists(successful_run))
  expect_true(dir.exists(success_parent))
  expect_identical(readLines(file.path(success_parent, "user-sentinel.txt")), "keep")

  failure_parent <- tempfile("integration-python-failure-parent-")
  dir.create(failure_parent)
  writeLines("keep", file.path(failure_parent, "user-sentinel.txt"))
  withr::defer(unlink(failure_parent, recursive = TRUE, force = TRUE))
  failed_run <- NULL
  expect_error(
    Shennong:::.sn_with_integration_python_run(
      method = "scvi",
      integration_control = list(run_dir = failure_parent, keep_run_dir = FALSE),
      code = function(control) {
        failed_run <<- control$run_dir
        dir.create(file.path(control$run_dir, "input"), recursive = TRUE)
        dir.create(file.path(control$run_dir, "output"), recursive = TRUE)
        writeLines("secret", file.path(control$run_dir, "input", "counts.mtx"))
        writeLines("secret", file.path(control$run_dir, "config.json"))
        writeLines("secret", file.path(control$run_dir, "output", "latent.csv"))
        stop("integration failed")
      }
    ),
    "integration failed.*Sanitized diagnostics remain"
  )
  expect_true(file.exists(file.path(failed_run, "failure.json")))
  expect_false(dir.exists(file.path(failed_run, "input")))
  expect_false(file.exists(file.path(failed_run, "config.json")))
  expect_false(file.exists(file.path(failed_run, "output", "latent.csv")))
  expect_true(dir.exists(failure_parent))
  expect_identical(readLines(file.path(failure_parent, "user-sentinel.txt")), "keep")
})

test_that("explicit integration run directories are retained", {
  run_dir <- tempfile("retained-integration-run-")
  value <- Shennong:::.sn_with_integration_python_run(
    method = "bbknn",
    integration_control = list(run_dir = run_dir),
    code = function(control) {
      writeLines("retained", file.path(control$run_dir, "marker.txt"))
      list(value = TRUE)
    }
  )
  expect_true(value$value)
  expect_true(file.exists(file.path(run_dir, "marker.txt")))
})

test_that("rare-cell dense backends reject oversized sparse inputs before conversion", {
  huge_sparse <- Matrix::sparseMatrix(
    i = integer(), j = integer(),
    dims = c(30000L, 100000L)
  )
  expect_error(
    Shennong:::.sn_assert_rare_dense_budget(
      huge_sparse,
      max_dense_gb = 2,
      name = "rare test"
    ),
    "exceeding `max_dense_gb"
  )
  expect_silent(
    Shennong:::.sn_assert_rare_dense_budget(
      Matrix::Matrix(matrix(0, nrow = 10, ncol = 10), sparse = TRUE),
      max_dense_gb = 0.01
    )
  )
})

test_that("generic Python imports fail closed on missing or empty output", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  run_dir <- tempfile("python-import-")
  output_dir <- file.path(run_dir, "output")
  dir.create(output_dir, recursive = TRUE)
  query <- .python_contract_query(object)

  invoke <- function() {
    Shennong:::.sn_import_python_object_results(
      object = object,
      method = "demo",
      result_name = "demo",
      output_dir = output_dir,
      run_dir = run_dir,
      assay = "RNA",
      query = query,
      return_object = FALSE
    )
  }

  expect_error(invoke(), "missing required `manifest.json`")
  .write_backend_manifest(output_dir, list(method = "other"))
  expect_error(invoke(), "manifest method mismatch")
  .write_backend_manifest(output_dir, list(method = "demo"))
  expect_error(invoke(), "without any importable")

  malformed <- data.frame(axis = 1:2, row.names = colnames(object)[1:2])
  utils::write.csv(malformed, file.path(output_dir, "latent.csv"))
  expect_error(invoke(), "Invalid Python embedding output")
})

test_that("cell2location reports exported and model feature counts separately", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  query <- .python_contract_query(object)
  run_dir <- tempfile("cell2location-import-")
  output_dir <- file.path(run_dir, "output")
  dir.create(output_dir, recursive = TRUE)
  abundance <- data.frame(
    cell_type_a = c(0.8, 0.2, 0.5),
    row.names = colnames(object)
  )
  utils::write.csv(abundance, file.path(output_dir, "obs.csv"))
  .write_backend_manifest(output_dir, list(
    method = "cell2location",
    n_cells = query$n_cells,
    n_features = query$n_features,
    n_shared_features = query$n_features - 1L
  ))

  imported <- Shennong:::.sn_import_python_object_results(
    object = object,
    method = "cell2location",
    result_name = "cell2location",
    output_dir = output_dir,
    run_dir = run_dir,
    assay = "RNA",
    query = query,
    return_object = FALSE
  )
  expect_identical(imported$n_features, query$n_features)
  expect_equal(imported$n_shared_features, query$n_features - 1L)
  expect_false(anyDuplicated(names(imported)) > 0L)

  script <- file.path(
    .python_contract_pixi_root(), "cell2location", "scripts", "cell2location_run.py"
  )
  source <- readLines(script, warn = FALSE)
  expect_true(any(grepl("input_n_features = int(adata.n_vars)", source, fixed = TRUE)))
  expect_true(any(grepl('"n_features": input_n_features', source, fixed = TRUE)))
  expect_true(any(grepl('"n_shared_features": int(adata.n_vars)', source, fixed = TRUE)))
  expect_true(any(grepl("_validate_raw_counts(matrix)", source, fixed = TRUE)))
  expect_true(any(grepl("_read_reference_signatures", source, fixed = TRUE)))
})

test_that("cell2location requires raw counts and validated reference signatures", {
  skip_if_not_installed("SeuratObject")
  object <- .make_python_contract_object()
  object$x <- seq_len(ncol(object))
  object$y <- rev(seq_len(ncol(object)))
  signatures <- matrix(
    1,
    nrow = nrow(object),
    ncol = 2L,
    dimnames = list(rownames(object), c("state_a", "state_b"))
  )

  expect_error(
    sn_run_cell2location(
      object,
      layer = "data",
      reference_signatures = signatures,
      spatial_cols = c("x", "y"),
      install_pixi = FALSE
    ),
    "raw/count-like input layer"
  )

  fractional <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  fractional@x[[1L]] <- fractional@x[[1L]] + 0.25
  SeuratObject::LayerData(object, assay = "RNA", layer = "fractional_counts") <- fractional
  expect_error(
    sn_run_cell2location(
      object,
      layer = "fractional_counts",
      reference_signatures = signatures,
      spatial_cols = c("x", "y"),
      install_pixi = FALSE
    ),
    "integer-like raw counts"
  )

  expect_s3_class(
    Shennong:::.sn_validate_cell2location_reference_signatures(signatures),
    "data.frame"
  )
  negative <- signatures
  negative[[1L]] <- -1
  expect_error(
    Shennong:::.sn_validate_cell2location_reference_signatures(negative),
    "finite and non-negative"
  )
  non_finite <- signatures
  non_finite[[1L]] <- Inf
  expect_error(
    Shennong:::.sn_validate_cell2location_reference_signatures(non_finite),
    "finite and non-negative"
  )
  duplicate_features <- signatures
  rownames(duplicate_features)[[2L]] <- rownames(duplicate_features)[[1L]]
  expect_error(
    Shennong:::.sn_validate_cell2location_reference_signatures(duplicate_features),
    "unique, non-empty"
  )
  duplicate_states <- signatures
  colnames(duplicate_states) <- rep("state", ncol(duplicate_states))
  expect_error(
    Shennong:::.sn_validate_cell2location_reference_signatures(duplicate_states),
    "unique, non-empty"
  )
  character_values <- as.data.frame(signatures)
  character_values[[1L]] <- as.character(character_values[[1L]])
  expect_error(
    Shennong:::.sn_validate_cell2location_reference_signatures(character_values),
    "must be numeric"
  )
})

test_that("direct scPoli requires integer-like raw counts before Python execution", {
  skip_if_not_installed("SeuratObject")
  object <- .make_python_contract_object()
  object$batch <- rep(c("a", "b", "a"), length.out = ncol(object))

  expect_error(
    sn_run_scpoli(
      object,
      layer = "data",
      batch_by = "batch",
      install_pixi = FALSE
    ),
    "raw/count-like input layer"
  )

  fractional <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  fractional@x[[1L]] <- fractional@x[[1L]] + 0.25
  SeuratObject::LayerData(object, assay = "RNA", layer = "fractional_counts") <- fractional
  expect_error(
    sn_run_scpoli(
      object,
      layer = "fractional_counts",
      batch_by = "batch",
      install_pixi = FALSE
    ),
    "integer-like raw counts"
  )
})

test_that("Tangram and Squidpy imports retain and require their core artifacts", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  query <- .python_contract_query(object)

  exercise <- function(method, field, filename) {
    run_dir <- tempfile(paste0(method, "-import-"))
    output_dir <- file.path(run_dir, "output")
    dir.create(output_dir, recursive = TRUE)
    .write_backend_manifest(output_dir, list(method = method))
    expect_error(
      Shennong:::.sn_import_python_object_results(
        object, method, method, output_dir, run_dir, "RNA", query,
        return_object = FALSE
      ),
      "without required"
    )

    artifact <- file.path(output_dir, filename)
    if (identical(method, "tangram")) {
      mapping <- matrix(
        1 / query$n_cells,
        nrow = query$n_cells,
        ncol = query$n_cells,
        dimnames = list(query$cells, query$cells)
      )
      utils::write.csv(mapping, artifact)
    } else {
      utils::write.csv(
        data.frame(
          source = query$cells[[1L]],
          target = query$cells[[2L]],
          weight = 1
        ),
        artifact,
        row.names = FALSE
      )
    }
    manifest <- list(method = method)
    manifest[[field]] <- artifact
    .write_backend_manifest(output_dir, manifest)
    imported <- Shennong:::.sn_import_python_object_results(
      object, method, method, output_dir, run_dir, "RNA", query,
      reference = if (identical(method, "tangram")) query else NULL,
      return_object = FALSE
    )
    expect_identical(unname(imported$imported_artifacts[[field]]), normalizePath(artifact))
  }

  exercise("tangram", "mapping_path", "mapping.csv")
  exercise("squidpy", "spatial_graph_path", "spatial_graph.csv")
})

test_that("retained Python artifacts are validated before success", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  query <- .python_contract_query(object)
  run_dir <- tempfile("retained-tangram-invalid-")
  output_dir <- file.path(run_dir, "output")
  dir.create(output_dir, recursive = TRUE)
  mapping_path <- file.path(output_dir, "mapping.csv")
  invalid <- matrix(
    1,
    nrow = query$n_cells,
    ncol = query$n_cells,
    dimnames = list(query$cells, c(query$cells[-1L], "unknown-cell"))
  )
  utils::write.csv(invalid, mapping_path)
  .write_backend_manifest(output_dir, list(method = "tangram", mapping_path = mapping_path))
  expect_error(
    Shennong:::.sn_import_python_object_results(
      object, "tangram", "tangram", output_dir, run_dir, "RNA", query,
      reference = query, retain_run_dir = TRUE, return_object = FALSE
    ),
    "do not exactly match"
  )

  invalid_probability <- matrix(
    0.2,
    nrow = query$n_cells,
    ncol = query$n_cells,
    dimnames = list(query$cells, query$cells)
  )
  utils::write.csv(invalid_probability, mapping_path)
  expect_error(
    Shennong:::.sn_import_python_object_results(
      object, "tangram", "tangram", output_dir, run_dir, "RNA", query,
      reference = query, retain_run_dir = TRUE, return_object = FALSE
    ),
    "probability distribution.*sums to 1"
  )
})

test_that("retained large-table artifacts stream validation without import budgets", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  query <- .python_contract_query(object)
  run_dir <- tempfile("retained-streamed-artifacts-")
  output_dir <- file.path(run_dir, "output")
  dir.create(output_dir, recursive = TRUE)
  withr::defer(unlink(run_dir, recursive = TRUE, force = TRUE))

  mapping_path <- file.path(output_dir, "mapping.csv")
  mapping <- matrix(
    1 / query$n_cells,
    nrow = query$n_cells,
    ncol = query$n_cells,
    dimnames = list(query$cells, query$cells)
  )
  utils::write.csv(mapping, mapping_path)
  tangram <- Shennong:::.sn_import_python_artifact_tables(
    method = "tangram",
    artifact_paths = c(mapping_path = mapping_path),
    query = query,
    reference = query,
    import = FALSE,
    max_import_gb = 1e-12
  )
  expect_identical(tangram, list())

  cellphonedb_path <- file.path(output_dir, "significant_means.txt")
  writeLines(c("id\tvalue", "interaction-1\t1"), cellphonedb_path)
  cellphonedb <- Shennong:::.sn_import_python_artifact_tables(
    method = "cellphonedb",
    artifact_paths = c(result_files = cellphonedb_path),
    query = query,
    reference = NULL,
    import = FALSE,
    max_import_gb = 1e-12
  )
  expect_identical(cellphonedb, list())
})

test_that("streamed artifact parsing chunks logical records instead of physical lines", {
  path <- tempfile(fileext = ".csv")
  withr::defer(unlink(path))
  writeLines(c(
    "id,value",
    "one,\"first line",
    "second line\"",
    "two,plain"
  ), path)
  observed <- list()

  Shennong:::.sn_walk_csv_rows(
    path = path,
    header = c("id", "value"),
    label = "multiline CSV fixture",
    chunk_rows = 1L,
    callback = function(chunk) observed[[length(observed) + 1L]] <<- chunk
  )
  combined <- do.call(rbind, observed)

  expect_identical(combined$id, c("one", "two"))
  expect_identical(combined$value, c("first line\nsecond line", "plain"))

  tab_path <- tempfile(fileext = ".txt")
  withr::defer(unlink(tab_path))
  writeLines(c(
    "id\tvalue",
    "interaction-1\t\"first line",
    "second line\""
  ), tab_path)
  expect_invisible(Shennong:::.sn_validate_cellphonedb_table_header(tab_path))
})

test_that("retained Squidpy enrichment artifacts are validated", {
  object <- .make_python_contract_object()
  query <- .python_contract_query(object)
  output_dir <- tempfile("retained-squidpy-invalid-")
  dir.create(output_dir)
  withr::defer(unlink(output_dir, recursive = TRUE, force = TRUE))
  graph_path <- file.path(output_dir, "spatial_graph.csv")
  enrichment_path <- file.path(output_dir, "neighborhood_enrichment.csv")
  utils::write.csv(
    data.frame(source = query$cells[[1L]], target = query$cells[[2L]], weight = 1),
    graph_path,
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(group_1 = "A", group_2 = "B", zscore = Inf),
    enrichment_path,
    row.names = FALSE
  )
  expect_error(
    Shennong:::.sn_import_python_artifact_tables(
      method = "squidpy",
      artifact_paths = c(
        spatial_graph_path = graph_path,
        neighborhood_enrichment_path = enrichment_path
      ),
      query = query,
      reference = NULL,
      import = FALSE
    ),
    "invalid `zscore`"
  )
})

test_that("disabled placeholder backends fail before launching pixi", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  expect_error(sn_run_scarches(object, install_pixi = FALSE), "faithful upstream scarches workflow")
  expect_error(sn_run_stlearn(object, install_pixi = FALSE), "faithful upstream stlearn workflow")
})

test_that("default Python runs export only required metadata and clean sensitive files", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  object$secret_patient_note <- paste0("private-", seq_len(ncol(object)))
  output_parent <- tempfile("generic-python-success-parent-")
  dir.create(output_parent)
  writeLines("keep", file.path(output_parent, "user-sentinel.txt"))
  withr::defer(unlink(output_parent, recursive = TRUE, force = TRUE))
  captured <- NULL

  updated <- testthat::with_mocked_bindings(
    Shennong:::.sn_run_python_object_method(
      object = object,
      environment = "demo",
      script_name = "unused.py",
      method = "demo",
      output_dir = output_parent,
      keep_run_dir = FALSE
    ),
    .sn_pixi_script_path = function(...) "unused.py",
    .sn_execute_python_object_pixi = function(input_dir, output_dir, ...) {
      exported <- utils::read.csv(file.path(input_dir, "query", "obs.csv"), check.names = FALSE)
      captured <<- list(run_dir = dirname(input_dir), columns = colnames(exported))
      metadata <- data.frame(score = seq_len(nrow(exported)), row.names = exported$cell_id)
      utils::write.csv(metadata, file.path(output_dir, "obs.csv"))
      .write_backend_manifest(output_dir, list(method = "demo", n_cells = nrow(exported)))
    },
    .package = "Shennong"
  )

  expect_identical(captured$columns, "cell_id")
  expect_false(dir.exists(captured$run_dir))
  expect_true(dir.exists(output_parent))
  expect_identical(readLines(file.path(output_parent, "user-sentinel.txt")), "keep")
  expect_false(updated@misc$demo$demo$run_dir_retained)
  expect_null(updated@misc$demo$demo$run_dir)
  expect_true("demo_score" %in% colnames(updated[[]]))
})

test_that("failed temporary Python runs retain sanitized diagnostics only", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  output_parent <- tempfile("generic-python-failure-parent-")
  dir.create(output_parent)
  writeLines("keep", file.path(output_parent, "user-sentinel.txt"))
  withr::defer(unlink(output_parent, recursive = TRUE, force = TRUE))
  captured <- NULL
  expect_error(
    testthat::with_mocked_bindings(
      Shennong:::.sn_run_python_object_method(
        object = object,
        environment = "demo",
        script_name = "unused.py",
        method = "demo",
        output_dir = output_parent,
        keep_run_dir = FALSE
      ),
      .sn_pixi_script_path = function(...) "unused.py",
      .sn_execute_python_object_pixi = function(input_dir, output_dir, ...) {
        captured <<- list(run_dir = dirname(input_dir), output_dir = output_dir)
        file.create(file.path(output_dir, "partial.h5ad"))
        writeLines("private", file.path(output_dir, "obs.csv"))
        stop("backend failed")
      },
      .package = "Shennong"
    ),
    "backend failed"
  )
  expect_true(file.exists(file.path(captured$run_dir, "failure.json")))
  expect_false(dir.exists(file.path(captured$run_dir, "input")))
  expect_false(file.exists(file.path(captured$output_dir, "partial.h5ad")))
  expect_false(file.exists(file.path(captured$output_dir, "obs.csv")))
  expect_length(list.files(captured$run_dir, pattern = "_config\\.json$"), 0L)
  expect_true(dir.exists(output_parent))
  expect_identical(readLines(file.path(output_parent, "user-sentinel.txt")), "keep")
})

test_that("scCAD runner serializes rare sets without a generated-script NameError", {
  python <- Sys.which("python3")
  skip_if(!nzchar(python), "python3 is required for generated scCAD runner test")
  module_dir <- tempfile("fake-sccad-")
  dir.create(module_dir)
  # Keep this contract test independent of the developer machine's optional
  # Python packages. The generated bridge only needs this tiny pandas surface.
  writeLines(character(), file.path(module_dir, "numpy.py"))
  writeLines(
    c(
      "import csv",
      "class _Axis(list):",
      "    def to_numpy(self): return list(self)",
      "class _Frame:",
      "    def __init__(self, rows, index, columns):",
      "        self._rows = rows",
      "        self.index = _Axis(index)",
      "        self.columns = _Axis(columns)",
      "    def to_numpy(self, dtype=float):",
      "        return [[dtype(value) for value in row] for row in self._rows]",
      "def read_csv(path, index_col=0):",
      "    with open(path, newline='') as handle: rows = list(csv.reader(handle))",
      "    return _Frame([row[1:] for row in rows[1:]], [row[0] for row in rows[1:]], rows[0][1:])"
    ),
    file.path(module_dir, "pandas.py")
  )
  module <- file.path(module_dir, "scCAD.py")
  writeLines(
    c(
      "def scCAD(**kwargs):",
      "    return [['cell1', 'cell2']], [0.9], ['rare_1', 'rare_1', 'major_1'], [['gene1']]"
    ),
    module
  )
  withr::local_envvar(PYTHONPATH = module_dir)
  expr <- Matrix::Matrix(matrix(1:9, nrow = 3), sparse = TRUE)
  output <- Shennong:::.sn_run_sccad(
    expr = expr,
    cell_ids = paste0("cell", 1:3),
    gene_ids = paste0("gene", 1:3),
    python = python,
    script = module
  )
  expect_equal(unname(unlist(output$rare_sets, use.names = FALSE)), c("cell1", "cell2"))
  expect_equal(output$scores, 0.9)
})

test_that("scVI-family and direct scPoli export use target-assay raw counts", {
  skip_if_not_installed("SeuratObject")
  object <- .make_python_contract_object()
  alt <- Matrix::Matrix(
    matrix(seq_len(12), nrow = 4, dimnames = list(paste0("alt", 1:4), colnames(object))),
    sparse = TRUE
  )
  object[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = alt)
  input_dir <- tempfile("scvi-target-assay-")
  withr::defer(unlink(input_dir, recursive = TRUE, force = TRUE))

  exported <- Shennong:::.sn_write_scvi_input(
    object = object,
    input_dir = input_dir,
    features = c("alt2", "alt4", "gene1"),
    assay = "ALT",
    layer = "counts",
    batch = NULL
  )
  expect_identical(exported$features, c("alt2", "alt4"))
  expect_identical(
    utils::read.csv(file.path(input_dir, "features.csv"), stringsAsFactors = FALSE)$feature_id,
    c("alt2", "alt4")
  )

  SeuratObject::LayerData(object, assay = "ALT", layer = "fractional_counts") <- alt + 0.25
  expect_error(
    Shennong:::.sn_write_scvi_input(
      object, tempfile("scvi-fractional-"), rownames(alt), "ALT",
      layer = "fractional_counts", batch = NULL
    ),
    "integer-like raw counts"
  )
  expect_error(
    Shennong:::.sn_write_scvi_input(
      object, tempfile("scvi-normalized-"), rownames(alt), "ALT",
      layer = "data", batch = NULL
    ),
    "raw/count-like RNA layer"
  )

  scpoli_script <- file.path(
    .python_contract_pixi_root(), "scarches", "scripts", "scpoli_integration.py"
  )
  scpoli_source <- readLines(scpoli_script, warn = FALSE)
  expect_true(any(grepl("_validate_raw_counts(matrix)", scpoli_source, fixed = TRUE)))
  expect_true(any(grepl("integer-like raw counts", scpoli_source, fixed = TRUE)))
})

test_that("Python result imports enforce budgets and exact cell identity", {
  skip_if_not_installed("Seurat")
  object <- .make_python_contract_object()
  query <- .python_contract_query(object)
  query$cells <- colnames(object)
  run_dir <- tempfile("strict-python-import-")
  output_dir <- file.path(run_dir, "output")
  dir.create(output_dir, recursive = TRUE)
  withr::defer(unlink(run_dir, recursive = TRUE, force = TRUE))
  .write_backend_manifest(output_dir, list(method = "demo", n_cells = ncol(object)))

  extra <- data.frame(score = seq_len(ncol(object) + 1L))
  rownames(extra) <- c(colnames(object), "extra-cell")
  utils::write.csv(extra, file.path(output_dir, "obs.csv"))
  expect_error(
    Shennong:::.sn_import_python_object_results(
      object, "demo", "demo", output_dir, run_dir, "RNA", query,
      return_object = FALSE
    ),
    "one-per-input|do not exactly match"
  )

  exact <- data.frame(score = seq_len(ncol(object)), row.names = colnames(object))
  utils::write.csv(exact, file.path(output_dir, "obs.csv"))
  expect_error(
    Shennong:::.sn_import_python_object_results(
      object, "demo", "demo", output_dir, run_dir, "RNA", query,
      return_object = FALSE, max_artifact_import_gb = 1e-12
    ),
    "exceeding `max_artifact_import_gb"
  )

  latent <- data.frame(axis = c(1, 2, Inf), row.names = colnames(object))
  utils::write.csv(latent, file.path(output_dir, "latent.csv"))
  expect_error(
    Shennong:::.sn_import_python_object_results(
      object, "demo", "demo", output_dir, run_dir, "RNA", query,
      return_object = FALSE
    ),
    "finite numeric values"
  )
})

test_that("stored Python manifests use a scalar allowlist without raw paths", {
  safe <- Shennong:::.sn_sanitize_python_backend_manifest(list(
    method = "demo",
    n_cells = 3L,
    backend_version = "1.2.3",
    output_h5ad = "/private/output.h5ad",
    result_files = c("a.csv", "b.csv"),
    versions = list(python = "3.12"),
    arbitrary = c("x", "y")
  ), retain_paths = TRUE)
  expect_identical(safe$method, "demo")
  expect_identical(safe$n_cells, 3L)
  expect_identical(safe$backend_version, "1.2.3")
  expect_true(safe$output_retained)
  expect_false(any(c("output_h5ad", "result_files", "versions", "arbitrary") %in% names(safe)))
})

test_that("CellPhoneDB artifact manifests accept only standard result tables", {
  output_dir <- tempfile("cellphonedb-results-")
  dir.create(output_dir)
  withr::defer(unlink(output_dir, recursive = TRUE, force = TRUE))
  resource <- file.path(output_dir, "resource-usage.txt")
  writeLines("private diagnostics", resource)
  expect_error(
    Shennong:::.sn_collect_python_artifacts(
      list(method = "cellphonedb", result_files = resource),
      output_dir
    ),
    "non-standard result artifact"
  )
  result <- file.path(output_dir, "significant_means.txt")
  writeLines("id\tvalue", result)
  artifacts <- Shennong:::.sn_collect_python_artifacts(
    list(method = "cellphonedb", result_files = result),
    output_dir
  )
  expect_identical(unname(artifacts), normalizePath(result))
})

test_that("temporary Python cleanup is verified before reporting success", {
  parent <- tempfile("python-cleanup-check-parent-")
  run_dir <- Shennong:::.sn_create_owned_run_dir(parent)
  withr::defer(unlink(parent, recursive = TRUE, force = TRUE))
  expect_error(
    Shennong:::.sn_remove_python_run_directory(
      run_dir,
      unlink_fn = function(...) 0L
    ),
    "not reported as cleaned"
  )
  expect_true(dir.exists(run_dir))

  unmarked <- file.path(parent, "unmarked")
  dir.create(unmarked)
  writeLines("keep", file.path(unmarked, "user-sentinel.txt"))
  expect_error(
    Shennong:::.sn_remove_python_run_directory(unmarked),
    "without a Shennong ownership marker"
  )
  expect_true(file.exists(file.path(unmarked, "user-sentinel.txt")))
  expect_error(
    Shennong:::.sn_sanitize_failed_python_run(unmarked, "demo", "execute"),
    "without a Shennong ownership marker"
  )
  expect_true(file.exists(file.path(unmarked, "user-sentinel.txt")))
})

test_that("infercnvpy controls fail before dependent stages can be launched", {
  object <- .make_python_contract_object()
  invoke <- function(...) {
    defaults <- list(
      object = object,
      reference_key = NULL,
      reference_cat = NULL,
      run_pca = TRUE,
      run_neighbors = TRUE,
      run_leiden = TRUE,
      run_umap = FALSE,
      score = TRUE,
      cnv_score_group_by = NULL
    )
    do.call(
      Shennong:::.sn_validate_infercnvpy_controls,
      utils::modifyList(defaults, list(...))
    )
  }
  expect_error(invoke(run_pca = FALSE), "run_neighbors.*requires.*run_pca")
  expect_error(invoke(run_neighbors = FALSE, run_leiden = TRUE), "run_leiden.*requires.*run_neighbors")
  expect_error(
    invoke(run_neighbors = FALSE, run_leiden = FALSE, run_umap = TRUE, score = FALSE),
    "run_umap.*requires.*run_neighbors"
  )
  expect_error(
    invoke(run_leiden = FALSE, score = TRUE),
    "score.*requires.*run_leiden"
  )
})

test_that("scCAD output contract rejects misaligned identities and scores", {
  valid <- list(
    cell_ids = as.list(c("c1", "c2")),
    rare_sets = list(list("c1")),
    scores = list(0.8),
    sub_clusters = as.list(c("rare", "major")),
    degs_list = list(list("g1"))
  )
  expect_no_error(Shennong:::.sn_validate_sccad_result(valid, c("c1", "c2"), c("g1", "g2")))
  invalid_cells <- valid
  invalid_cells$cell_ids <- rev(invalid_cells$cell_ids)
  expect_error(
    Shennong:::.sn_validate_sccad_result(invalid_cells, c("c1", "c2"), c("g1", "g2")),
    "exactly match"
  )
  invalid_scores <- valid
  invalid_scores$scores <- list()
  expect_error(
    Shennong:::.sn_validate_sccad_result(invalid_scores, c("c1", "c2"), c("g1", "g2")),
    "one finite score"
  )
})
