.sn_run_scdesign3_backend <- function(...) {
  scDesign3::scdesign3(...)
}

.sn_seurat_to_sce_for_scdesign3 <- function(object,
                                            assay = "RNA",
                                            layer = "counts") {
  counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  SingleCellExperiment::SingleCellExperiment(
    list(counts = .sn_as_sparse_matrix(counts)),
    colData = object[[]]
  )
}

.sn_validate_scdesign3_columns <- function(sce,
                                           celltype = NULL,
                                           pseudotime = NULL,
                                           spatial = NULL,
                                           other_covariates = NULL) {
  available <- colnames(SummarizedExperiment::colData(sce))
  requested <- unique(c(celltype, pseudotime, spatial, other_covariates))
  requested <- requested[!is.na(requested) & nzchar(requested)]
  missing <- setdiff(requested, available)
  if (length(missing) > 0L) {
    stop(glue("Missing scDesign3 covariate column(s): {paste(missing, collapse = ', ')}."), call. = FALSE)
  }
  invisible(TRUE)
}

.sn_default_scdesign3_formula <- function(celltype = NULL,
                                          pseudotime = NULL,
                                          spatial = NULL,
                                          other_covariates = NULL,
                                          k = 10) {
  if (!is.null(pseudotime) && length(pseudotime) > 0L) {
    return(glue("s({pseudotime[[1]]}, bs = 'cr', k = {k})"))
  }
  if (!is.null(spatial) && length(spatial) == 2L) {
    return(glue("s({spatial[[1]]}, {spatial[[2]]}, bs = 'gp', k = {k})"))
  }
  terms <- unique(c(celltype, other_covariates))
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) == 0L) {
    return("1")
  }
  paste(terms, collapse = " + ")
}

.sn_default_scdesign3_corr_formula <- function(celltype = NULL,
                                               pseudotime = NULL,
                                               spatial = NULL,
                                               other_covariates = NULL) {
  terms <- unique(c(pseudotime, spatial, celltype, other_covariates))
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) == 0L) {
    return("1")
  }
  paste(terms, collapse = " + ")
}

.sn_extract_scdesign3_counts <- function(result, expected_genes = NULL) {
  if (!is.list(result) || !"new_count" %in% names(result)) {
    stop("scDesign3 did not return a `new_count` result.", call. = FALSE)
  }
  counts <- result$new_count
  if (is.list(counts)) {
    counts <- counts[[1]]
  }
  counts <- .sn_as_sparse_matrix(counts)
  if (is.null(rownames(counts)) && !is.null(expected_genes) && length(expected_genes) == nrow(counts)) {
    rownames(counts) <- expected_genes
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("sim_cell_", seq_len(ncol(counts)))
  }
  counts
}

#' Simulate single-cell counts with scDesign3
#'
#' \code{sn_simulate()} provides a method-based simulation entry point.
#' Currently \code{method = "scdesign3"} delegates to
#' \code{.sn_simulate_scdesign3()}.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param method Simulation backend. Currently supports \code{"scdesign3"}.
#' @param ... Additional arguments passed to the selected backend.
#'
#' @return Simulated data in the backend's requested format.
#'
#' @examples
#' \dontrun{
#' sim <- sn_simulate(seurat_obj, method = "scdesign3", celltype = "cell_type")
#' }
#' @export
sn_simulate <- function(object,
                        method = c("scdesign3"),
                        ...) {
  method <- match.arg(method)
  switch(
    method,
    scdesign3 = .sn_simulate_scdesign3(object = object, ...)
  )
}

.sn_simulate_scdesign3 <- function(object,
                                  celltype = NULL,
                                  pseudotime = NULL,
                                  spatial = NULL,
                                  other_covariates = NULL,
                                  ncell = NULL,
                                  mu_formula = NULL,
                                  sigma_formula = "1",
                                  corr_formula = NULL,
                                  family_use = "nb",
                                  n_cores = 2,
                                  assay = "RNA",
                                  layer = "counts",
                                  assay_use = "counts",
                                  return = c("seurat", "sce", "counts", "result"),
                                  project = "scdesign3",
                                  combine_original = FALSE,
                                  seed = 717,
                                  ...) {
  check_installed("scDesign3", reason = "to simulate single-cell data with scDesign3.")
  check_installed(c("SingleCellExperiment", "SummarizedExperiment"))

  return <- match.arg(return)
  if (inherits(object, "Seurat")) {
    sce <- .sn_seurat_to_sce_for_scdesign3(object = object, assay = assay, layer = layer)
  } else if (inherits(object, "SingleCellExperiment")) {
    sce <- object
  } else {
    stop("`object` must be a Seurat or SingleCellExperiment object.", call. = FALSE)
  }

  coldata <- SummarizedExperiment::colData(sce)
  if (is.null(celltype) && "seurat_clusters" %in% colnames(coldata)) {
    celltype <- "seurat_clusters"
  }
  .sn_validate_scdesign3_columns(
    sce = sce,
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates
  )

  ncell <- ncell %||% ncol(sce)
  formula_k <- max(3L, min(10L, floor(ncell / 5)))
  mu_formula <- mu_formula %||% .sn_default_scdesign3_formula(
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    k = formula_k
  )
  corr_formula <- corr_formula %||% .sn_default_scdesign3_corr_formula(
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }
  result <- .sn_run_scdesign3_backend(
    sce = sce,
    assay_use = assay_use,
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    ncell = ncell,
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    family_use = family_use,
    n_cores = n_cores,
    corr_formula = corr_formula,
    ...
  )
  result$shennong <- list(
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    corr_formula = corr_formula,
    family_use = family_use,
    ncell = ncell
  )

  if (return == "result") {
    return(result)
  }

  sim_counts <- .sn_extract_scdesign3_counts(
    result = result,
    expected_genes = rownames(sce)
  )
  sim_meta <- result$new_covariate
  if (is.null(sim_meta)) {
    sim_meta <- data.frame(row.names = colnames(sim_counts))
  } else {
    sim_meta <- as.data.frame(sim_meta)
    if (nrow(sim_meta) == ncol(sim_counts)) {
      rownames(sim_meta) <- colnames(sim_counts)
    }
  }

  if (return == "counts") {
    return(sim_counts)
  }
  if (return == "sce") {
    return(SingleCellExperiment::SingleCellExperiment(
      list(counts = sim_counts),
      colData = sim_meta
    ))
  }

  check_installed("Seurat")
  sim_object <- Seurat::CreateSeuratObject(
    counts = sim_counts,
    meta.data = sim_meta,
    project = project
  )
  sim_object@misc$scdesign3 <- result$shennong
  if (isTRUE(combine_original) && inherits(object, "Seurat")) {
    original <- object
    original$simulation_source <- "original"
    sim_object$simulation_source <- "simulated"
    combined <- .sn_with_default_acceleration(
      merge(
        x = original,
        y = sim_object,
        add.cell.ids = c("original", "simulated"),
        project = project
      ),
      patches = "seurat_merge",
      strict = TRUE,
      operation = "merge"
    )
    combined@misc$scdesign3 <- result$shennong
    return(combined)
  }
  sim_object
}

#' Deprecated direct scDesign3 wrapper
#'
#' `sn_simulate_scdesign3()` is a deprecated compatibility entry point. Use
#' [sn_simulate(method = "scdesign3")] instead.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param ... Additional arguments passed to the scDesign3 backend.
#'
#' @return Simulated data in the backend's requested format.
#'
#' @export
sn_simulate_scdesign3 <- function(object, ...) {
  .Deprecated("sn_simulate(method = \"scdesign3\")", package = "Shennong")
  .sn_simulate_scdesign3(object = object, ...)
}
