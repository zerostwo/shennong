library(testthat)

# AUDIT-06: make inst/methods/*.yml load-bearing.
#
# Direction A: every registry entry must be accepted by its owning function's
# dispatch surface (match.arg choice set) unless explicitly waived.
# Direction B: every dispatch value must be present in the registry unless
# explicitly waived with a reason.
#
# The hardcoded surfaces that predate the registry (integration,
# normalization, doublets, pseudobulk DE) are pinned by an exact snapshot so
# any change fails closed until the snapshot is consciously updated.

.sn_parity_choice_set <- function(fn_name) {
  fun <- tryCatch(getExportedValue("Shennong", fn_name), error = function(e) NULL)
  if (is.null(fun)) {
    fun <- getFromNamespace(fn_name, "Shennong")
  }
  unname(eval(formals(fun)$method))
}

.sn_registry_entries <- function(file) {
  path <- system.file(file.path("methods", file), package = "Shennong")
  if (!nzchar(path)) {
    path <- file.path("..", "..", "inst", "methods", file)
  }
  entries <- jsonlite::fromJSON(path, simplifyDataFrame = FALSE)
  vapply(entries, `[[`, character(1), "name")
}

test_that("registry entries are accepted by their owning functions", {
  owners <- list(
    list(file = "annotation.yml", task = "annotation", fn = "sn_run_annotation"),
    list(file = "bulk.yml", task = "bulk_de", fn = ".sn_find_bulk_de"),
    list(file = "bulk.yml", task = "bulk_deconvolution", fn = "sn_run_bulk_deconvolution"),
    list(file = "cnv.yml", task = "cnv", fn = "sn_run_cnv"),
    list(file = "communication.yml", task = "communication", fn = "sn_run_cell_communication"),
    list(file = "differential_abundance.yml", task = "differential_abundance", fn = "sn_test_abundance"),
    list(file = "grn.yml", task = "grn", fn = "sn_run_grn"),
    list(file = "metabolism.yml", task = "metabolism", fn = "sn_run_metabolism"),
    list(file = "program_discovery.yml", task = "program_discovery", fn = "sn_discover_programs"),
    list(file = "program_scoring.yml", task = "program_scoring", fn = "sn_score_programs"),
    list(file = "spatial.yml", task = "spatial_svg", fn = "sn_find_spatial_features"),
    list(file = "spatial.yml", task = "spatial_domain", fn = "sn_find_spatial_domains"),
    list(file = "spatial.yml", task = "spatial_neighborhood", fn = "sn_run_spatial_neighborhood"),
    list(file = "state_priority.yml", task = "state_priority", fn = "sn_prioritize_states"),
    list(file = "trajectory.yml", task = "trajectory", fn = "sn_run_trajectory"),
    list(file = "trajectory.yml", task = "velocity", fn = "sn_run_velocity"),
    list(file = "trajectory.yml", task = "fate", fn = "sn_run_fate")
  )
  # task::name values accepted outside a `method=` selector or pending a
  # conformance contract; each needs a concrete reason.
  waivers <- c(
    "tradeseq",        # dynamic-genes stage toggled by test_dynamic, not method=
    "cibersortx",      # container CLI tool, not a managed runtime backend yet
    "cell2location",   # reached through deprecated spatial deconvolution shim
    "tangram"          # reached through deprecated spatial mapping shim
  )
  files <- unique(vapply(owners, `[[`, character(1), "file"))
  for (owner in owners) {
    names <- .sn_registry_entries(owner$file)
    tasks <- vapply(
      jsonlite::fromJSON(system.file(file.path("methods", owner$file), package = "Shennong"), simplifyDataFrame = FALSE),
      `[[`, character(1), "task"
    )
    entry_names <- names[tasks == owner$task]
    expect_true(length(entry_names) > 0L, info = paste("task has entries:", owner$task))
    if (identical(owner$fn, ".sn_find_bulk_de")) {
      choices <- unname(eval(formals(getFromNamespace(".sn_find_bulk_de", "Shennong"))$method))
    } else {
      choices <- .sn_parity_choice_set(owner$fn)
    }
    expect_false(is.null(choices), info = paste(owner$fn, "has no method= selector"))
    outside <- setdiff(entry_names, choices)
    outside <- setdiff(outside, waivers)
    expect_true(length(outside) == 0L, info = paste0(
      owner$task, ": registry methods missing from ", owner$fn, " dispatch: ",
      paste(outside, collapse = ", ")
    ))
  }
})

test_that("dispatch values are covered by the registry", {
  coverage <- list(
    list(task = "annotation", file = "annotation.yml", fn = "sn_run_annotation"),
    list(task = "bulk_de", file = "bulk.yml", fn = ".sn_find_bulk_de"),
    list(task = "communication", file = "communication.yml", fn = "sn_run_cell_communication"),
    list(task = "cnv", file = "cnv.yml", fn = "sn_run_cnv"),
    list(task = "differential_abundance", file = "differential_abundance.yml", fn = "sn_test_abundance"),
    list(task = "grn", file = "grn.yml", fn = "sn_run_grn"),
    list(task = "metabolism", file = "metabolism.yml", fn = "sn_run_metabolism"),
    list(task = "program_discovery", file = "program_discovery.yml", fn = "sn_discover_programs"),
    list(task = "program_scoring", file = "program_scoring.yml", fn = "sn_score_programs"),
    list(task = "spatial_svg", file = "spatial.yml", fn = "sn_find_spatial_features"),
    list(task = "spatial_domain", file = "spatial.yml", fn = "sn_find_spatial_domains"),
    list(task = "spatial_neighborhood", file = "spatial.yml", fn = "sn_run_spatial_neighborhood"),
    list(task = "state_priority", file = "state_priority.yml", fn = "sn_prioritize_states"),
    list(task = "trajectory", file = "trajectory.yml", fn = "sn_run_trajectory"),
    list(task = "velocity", file = "trajectory.yml", fn = "sn_run_velocity"),
    list(task = "fate", file = "trajectory.yml", fn = "sn_run_fate")
  )
  # Surfaces deliberately absent from inst/methods/ today. Adding entries here
  # requires an admitted conformance contract per AGENTS.md, so these remain
  # explicit, reviewed gaps instead of silent drift.
  unregistered <- list(
    list(fn = NULL, choices = c("dorothea", "progeny"),
         reason = "regulatory activity has no registry file yet"),
    list(fn = "sn_test_programs", choices = NULL,
         reason = "program comparison tests have no registry file yet")
  )
  for (item in coverage) {
    if (identical(item$fn, ".sn_find_bulk_de")) {
      choices <- unname(eval(formals(getFromNamespace(".sn_find_bulk_de", "Shennong"))$method))
    } else {
      choices <- .sn_parity_choice_set(item$fn)
    }
    names <- .sn_registry_entries(item$file)
    tasks <- vapply(
      jsonlite::fromJSON(system.file(file.path("methods", item$file), package = "Shennong"), simplifyDataFrame = FALSE),
      `[[`, character(1), "task"
    )
    entry_names <- names[tasks == item$task]
    missing <- setdiff(choices, entry_names)
    # 'sct' style lowercase aliases resolve to canonical registry names.
    missing <- setdiff(tolower(missing), tolower(entry_names))
    expect_true(length(missing) == 0L, info = paste0(
      item$task, ": dispatch values missing from registry: ",
      paste(missing, collapse = ", ")
    ))
  }
})

test_that("hardcoded method surfaces stay pinned (fail-closed harvest)", {
  cluster_formals <- formals(getExportedValue("Shennong", "sn_run_cluster"))
  integration_choices <- eval(cluster_formals$integration_method)
  normalization_choices <- eval(cluster_formals$normalization_method)
  expect_identical(
    integration_choices,
    c("harmony", "unintegrated", "coralysis", "seurat_cca", "seurat_rpca",
      "scvi", "scanvi", "scpoli", "bbknn", "totalvi", "mmochi")
  )
  expect_identical(normalization_choices, c("seurat", "scran", "sctransform"))

  doublet_choices <- eval(formals(getExportedValue("Shennong", "sn_find_doublets"))$method)
  expect_identical(doublet_choices, c("scdblfinder", "scrublet"))

  pseudobulk_choices <- eval(formals(getFromNamespace(".sn_run_pseudobulk_de", "Shennong"))$method)
  expect_identical(pseudobulk_choices, c("DESeq2", "edgeR", "limma"))
})
