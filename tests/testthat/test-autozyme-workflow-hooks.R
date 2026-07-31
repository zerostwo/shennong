.sn_autozyme_workflow_body <- function(name) {
  fun <- get(name, envir = asNamespace("Shennong"), inherits = FALSE)
  paste(deparse(body(fun), width.cutoff = 500L), collapse = "\n")
}

.sn_autozyme_fixed_positions <- function(text, pattern) {
  positions <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (identical(positions, -1L)) integer(0) else positions
}

.sn_expect_autozyme_wrapper_before <- function(name,
                                               wrapper,
                                               targets,
                                               minimum_wrappers = 1L,
                                               patch = NULL) {
  text <- .sn_autozyme_workflow_body(name)
  wrappers <- .sn_autozyme_fixed_positions(
    text,
    wrapper
  )

  expect_true(length(wrappers) >= minimum_wrappers, info = name)
  if (!is.null(patch)) {
    expect_true(
      grepl(paste0('"', patch, '"'), text, fixed = TRUE),
      info = paste(name, patch)
    )
  }

  for (target in targets) {
    target_positions <- .sn_autozyme_fixed_positions(text, target)
    expect_true(length(target_positions) > 0L, info = paste(name, target))
    expect_true(
      all(vapply(
        target_positions,
        function(position) any(wrappers < position),
        logical(1)
      )),
      info = paste(name, "must enter the AutoZyme wrapper before", target)
    )
  }
}

.sn_autozyme_target_is_wrapped <- function(name, wrapper, target) {
  fun <- get(name, envir = asNamespace("Shennong"), inherits = FALSE)

  inspect_call <- function(expr, inside_wrapper = FALSE) {
    if (!is.call(expr)) return(FALSE)

    call_name <- paste(deparse(expr[[1L]], width.cutoff = 500L), collapse = "")
    inside_wrapper <- inside_wrapper || identical(call_name, wrapper)
    if (inside_wrapper && identical(call_name, target)) return(TRUE)

    arguments <- as.list(expr)[-1L]
    length(arguments) > 0L && any(vapply(
      arguments,
      inspect_call,
      logical(1),
      inside_wrapper = inside_wrapper
    ))
  }

  inspect_call(body(fun))
}

test_that("enrichment, WGCNA, and tradeSeq use scoped default patches", {
  enrichment <- .sn_autozyme_workflow_body("sn_enrich")
  enrichment_wrapper <- .sn_autozyme_fixed_positions(
    enrichment,
    ".sn_with_default_autozyme"
  )
  expect_length(enrichment_wrapper, 2L)
  expect_length(.sn_autozyme_fixed_positions(enrichment, '"clusterprofiler"'), 2L)
  expect_length(.sn_autozyme_fixed_positions(enrichment, '"fgsea"'), 2L)
  expect_length(.sn_autozyme_fixed_positions(enrichment, "with_enrichment_autozyme"), 10L)

  for (target in c(
    "clusterProfiler::gseGO", "clusterProfiler::enrichGO",
    "clusterProfiler::gseKEGG", "clusterProfiler::enrichKEGG",
    "clusterProfiler::GSEA", "clusterProfiler::enricher",
    "clusterProfiler::compareCluster"
  )) {
    positions <- .sn_autozyme_fixed_positions(enrichment, target)
    expect_true(length(positions) > 0L, info = target)
    expect_true(all(vapply(
      positions,
      function(position) {
        local_wrapper <- .sn_autozyme_fixed_positions(
          substr(enrichment, 1L, position),
          "with_enrichment_autozyme"
        )
        length(local_wrapper) > 1L
      },
      logical(1)
    )), info = target)
  }
  expect_true(
    max(enrichment_wrapper) <
      .sn_autozyme_fixed_positions(enrichment, "sn_store_enrichment")[[1L]]
  )

  .sn_expect_autozyme_wrapper_before(
    "sn_run_wgcna",
    ".sn_with_default_autozyme",
    c("WGCNA::goodSamplesGenes", "WGCNA::blockwiseModules", ".sn_bulk_result"),
    minimum_wrappers = 1L,
    patch = "wgcna"
  )
  .sn_expect_autozyme_wrapper_before(
    ".sn_fit_trajectory_dynamics",
    ".sn_with_default_autozyme",
    "tradeSeq::fitGAM",
    patch = "tradeseq"
  )
  .sn_expect_autozyme_wrapper_before(
    "sn_run_trajectory",
    ".sn_with_default_autozyme",
    c(".sn_fit_trajectory_dynamics", ".sn_analysis_provenance"),
    patch = "tradeseq"
  )
  .sn_expect_autozyme_wrapper_before(
    "sn_run_trajectory",
    ".sn_with_autozyme_provenance_context",
    c(".sn_trajectory_slingshot", ".sn_analysis_provenance"),
    patch = "slingshot"
  )
})

test_that("Seurat preprocessing and DE enable acceleration only on Seurat paths", {
  .sn_expect_autozyme_wrapper_before(
    "sn_normalize_data",
    ".sn_with_default_seurat_autozyme",
    c("Seurat::NormalizeData", "Seurat::SCTransform"),
    minimum_wrappers = 2L
  )
  .sn_expect_autozyme_wrapper_before(
    ".sn_run_seurat_de",
    ".sn_with_default_seurat_autozyme",
    c("Seurat::FindAllMarkers", "Seurat::FindMarkers"),
    minimum_wrappers = 2L
  )
  .sn_expect_autozyme_wrapper_before(
    "sn_find_de",
    ".sn_with_autozyme_provenance_context",
    c(".sn_run_seurat_de", ".sn_contextual_analysis_provenance"),
    patch = "seurat"
  )
  .sn_expect_autozyme_wrapper_before(
    ".sn_run_nichenetr",
    ".sn_with_default_seurat_autozyme",
    "Seurat::FindMarkers"
  )
})

test_that("Seurat clustering, priority, and plotting hooks precede patched targets", {
  .sn_expect_autozyme_wrapper_before(
    ".sn_run_seurat_layer_integration",
    ".sn_with_default_seurat_autozyme",
    c("Seurat::CCAIntegration", "Seurat::IntegrateLayers"),
    minimum_wrappers = 2L
  )
  .sn_expect_autozyme_wrapper_before(
    ".sn_make_temporary_grouping",
    ".sn_with_default_seurat_autozyme",
    c("Seurat::ScaleData", "Seurat::RunPCA"),
    minimum_wrappers = 2L
  )
  .sn_expect_autozyme_wrapper_before(
    ".sn_select_variable_features",
    ".sn_with_default_seurat_autozyme",
    "Seurat::FindVariableFeatures",
    minimum_wrappers = 2L
  )
  .sn_expect_autozyme_wrapper_before(
    ".sn_run_cluster_impl",
    ".sn_with_default_seurat_autozyme",
    c(
      "Seurat::SCTransform", "Seurat::NormalizeData",
      "Seurat::ScaleData", "Seurat::RunPCA"
    ),
    minimum_wrappers = 8L
  )
  .sn_expect_autozyme_wrapper_before(
    ".sn_get_seurat_logcounts_sce",
    ".sn_with_default_seurat_autozyme",
    "Seurat::NormalizeData"
  )
  .sn_expect_autozyme_wrapper_before(
    ".sn_priority_scissor_impl",
    ".sn_with_default_seurat_autozyme",
    c("Seurat::ScaleData", "Seurat::RunPCA"),
    minimum_wrappers = 2L
  )
  .sn_expect_autozyme_wrapper_before(
    ".sn_scale_heatmap_features",
    ".sn_with_default_seurat_autozyme",
    "Seurat::ScaleData"
  )

  for (name in c(
    ".sn_make_temporary_grouping",
    ".sn_run_cluster_impl",
    ".sn_priority_rareq",
    ".sn_priority_scissor_impl"
  )) {
    expect_false(
      .sn_autozyme_target_is_wrapped(
        name,
        ".sn_with_default_seurat_autozyme",
        "Seurat::FindNeighbors"
      ),
      info = paste(name, "must not wrap disabled FindNeighbors patch target")
    )
  }
})

test_that("the Seurat wrapper suspends AutoZyme for BPCells-backed objects", {
  uses_bpcells <- TRUE
  disabled_calls <- 0L
  scoped_calls <- 0L
  executions <- 0L

  testthat::local_mocked_bindings(
    .sn_seurat_uses_bpcells = function(object, assay = NULL) uses_bpcells,
    .sn_with_autozyme_disabled = function(expr) {
      disabled_calls <<- disabled_calls + 1L
      force(expr)
    },
    .sn_with_default_autozyme = function(expr, patches) {
      scoped_calls <<- scoped_calls + 1L
      expect_identical(patches, "seurat")
      force(expr)
    },
    .package = "Shennong"
  )

  value <- Shennong:::.sn_with_default_seurat_autozyme(
    {
      executions <- executions + 1L
      41L
    },
    object = structure(list(), class = "mock_seurat"),
    assay = "RNA"
  )
  expect_identical(value, 41L)
  expect_identical(executions, 1L)
  expect_identical(disabled_calls, 1L)
  expect_identical(scoped_calls, 0L)

  uses_bpcells <- FALSE
  value <- Shennong:::.sn_with_default_seurat_autozyme(
    {
      executions <- executions + 1L
      42L
    },
    object = structure(list(), class = "mock_seurat"),
    assay = "RNA"
  )
  expect_identical(value, 42L)
  expect_identical(executions, 2L)
  expect_identical(disabled_calls, 1L)
  expect_identical(scoped_calls, 1L)
})

test_that("BPCells-backed Seurat assays are detected before acceleration", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    matrix(
      as.integer((seq_len(48L) %% 5L) > 0L),
      nrow = 8L,
      ncol = 6L
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))

  path <- tempfile("shennong-autozyme-bpcells-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  BPCells::write_matrix_dir(
    BPCells::convert_matrix_type(counts, "uint32_t"),
    dir = path
  )

  bpcells_object <- SeuratObject::CreateSeuratObject(
    counts = BPCells::open_matrix_dir(path)
  )
  memory_object <- SeuratObject::CreateSeuratObject(counts = counts)

  expect_true(
    Shennong:::.sn_seurat_uses_bpcells(bpcells_object, assay = "RNA")
  )
  expect_false(
    Shennong:::.sn_seurat_uses_bpcells(memory_object, assay = "RNA")
  )
})
