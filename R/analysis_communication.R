.sn_validate_seurat_object <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }
  invisible(TRUE)
}

.sn_expression_matrix <- function(object,
                                  assay = NULL,
                                  layer = "data",
                                  features = NULL,
                                  cells = NULL) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  mat <- SeuratObject::LayerData(object = object, assay = assay, layer = layer)
  if (!is.null(features)) {
    features <- intersect(features, rownames(mat))
    mat <- mat[features, , drop = FALSE]
  }
  if (!is.null(cells)) {
    cells <- intersect(cells, colnames(mat))
    mat <- mat[, cells, drop = FALSE]
  }
  mat
}

.sn_resolve_species <- function(object, species = NULL) {
  species <- species %||% tryCatch(sn_get_species(object), error = function(...) NULL)
  species <- tolower(species %||% "human")
  if (!species %in% c("human", "mouse")) {
    stop("`species` must be either 'human' or 'mouse' for this backend.", call. = FALSE)
  }
  species
}

.sn_group_average_matrix <- function(mat, groups) {
  groups <- as.factor(groups)
  keep <- !is.na(groups)
  mat <- mat[, keep, drop = FALSE]
  groups <- droplevels(groups[keep])
  group_levels <- levels(groups)
  out <- vapply(group_levels, function(level) {
    cells <- which(groups == level)
    Matrix::rowMeans(mat[, cells, drop = FALSE])
  }, numeric(nrow(mat)))
  rownames(out) <- rownames(mat)
  colnames(out) <- group_levels
  out
}

.sn_expressed_genes <- function(mat, cells, min_pct = 0.1) {
  cells <- intersect(cells, colnames(mat))
  if (length(cells) == 0L) {
    return(character())
  }
  pct <- Matrix::rowMeans(mat[, cells, drop = FALSE] > 0)
  names(pct)[pct >= min_pct]
}

.sn_lr_columns <- function(lr_network) {
  ligand_col <- .sn_communication_column(lr_network, c("ligand", "from", "source"))
  receptor_col <- .sn_communication_column(lr_network, c("receptor", "to", "target"))
  if (is.null(ligand_col) || is.null(receptor_col)) {
    stop("`lr_network` must contain ligand/receptor columns such as `ligand`/`receptor` or `from`/`to`.", call. = FALSE)
  }
  list(ligand = ligand_col, receptor = receptor_col)
}

.sn_run_cellchat <- function(object,
                             group_by,
                             assay = NULL,
                             layer = "data",
                             species = NULL,
                             cellchat_db = NULL,
                             min_cells = 10,
                             population_size = FALSE,
                             raw_use = TRUE,
                             ...) {
  check_installed("CellChat", reason = "to run CellChat cell-cell communication inference.")
  species <- .sn_resolve_species(object, species)
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  mat <- .sn_expression_matrix(object = object, assay = assay, layer = layer)
  meta <- object[[]][colnames(mat), , drop = FALSE]
  meta$shennong_group <- as.character(meta[[group_by]])
  if (!"samples" %in% names(meta)) meta$samples <- "sample1"
  keep <- !is.na(meta$shennong_group) & nzchar(meta$shennong_group)
  mat <- mat[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  db <- cellchat_db
  if (is_null(db)) {
    database_name <- if (identical(species, "human")) "CellChatDB.human" else "CellChatDB.mouse"
    database_environment <- new.env(parent = emptyenv())
    utils::data(list = database_name, package = "CellChat", envir = database_environment)
    db <- database_environment[[database_name]]
    if (is_null(db)) stop("Could not load the CellChat database `", database_name, "`.", call. = FALSE)
  }
  backend_warnings <- character()
  cellchat <- withCallingHandlers(
    {
      current <- CellChat::createCellChat(
        object = mat,
        meta = meta,
        group.by = "shennong_group",
        assay = assay
      )
      current@DB <- db
      current <- CellChat::subsetData(current)
      current <- CellChat::identifyOverExpressedGenes(current)
      current <- CellChat::identifyOverExpressedInteractions(current)
      current <- CellChat::computeCommunProb(
        current,
        raw.use = raw_use,
        population.size = population_size,
        ...
      )
      current <- CellChat::filterCommunication(current, min.cells = min_cells)
      current <- CellChat::computeCommunProbPathway(current)
      CellChat::aggregateNet(current)
    },
    warning = function(warning) {
      backend_warnings <<- c(backend_warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )

  table <- tibble::as_tibble(CellChat::subsetCommunication(cellchat))
  list(
    table = table,
    backend_result = cellchat,
    warnings = unique(backend_warnings),
    artifacts = list(cellchat = cellchat)
  )
}

.sn_run_nichenetr <- function(object,
                              group_by,
                              assay = NULL,
                              layer = "data",
                              sender,
                              receiver,
                              geneset = NULL,
                              background_genes = NULL,
                              condition_col = NULL,
                              condition_oi = NULL,
                              condition_reference = NULL,
                              ligand_target_matrix = NULL,
                              lr_network = NULL,
                              expressed_pct = 0.1,
                              top_n = 50,
                              ...) {
  check_installed("nichenetr", reason = "to run NicheNet ligand-activity inference.")
  if (is.null(ligand_target_matrix)) {
    stop("`ligand_target_matrix` is required for `method = \"nichenetr\"`.", call. = FALSE)
  }
  if (is.null(lr_network)) {
    stop("`lr_network` is required for `method = \"nichenetr\"`.", call. = FALSE)
  }
  if (is.null(sender) || is.null(receiver)) {
    stop("`sender` and `receiver` are required for `method = \"nichenetr\"`.", call. = FALSE)
  }

  mat <- .sn_expression_matrix(object = object, assay = assay, layer = layer)
  meta <- object[[]][colnames(mat), , drop = FALSE]
  sender_cells <- rownames(meta)[as.character(meta[[group_by]]) %in% sender]
  receiver_cells <- rownames(meta)[as.character(meta[[group_by]]) %in% receiver]

  if (is.null(geneset)) {
    if (is.null(condition_col) || is.null(condition_oi) || is.null(condition_reference)) {
      stop("Supply `geneset`, or supply `condition_col`, `condition_oi`, and `condition_reference` to derive receiver DE genes.", call. = FALSE)
    }
    receiver_object <- subset(object, cells = receiver_cells)
    Seurat::Idents(receiver_object) <- receiver_object[[condition_col, drop = TRUE]]
    de <- Seurat::FindMarkers(
      receiver_object,
      ident.1 = condition_oi,
      ident.2 = condition_reference,
      assay = assay,
      slot = layer,
      ...
    )
    de$gene <- rownames(de)
    p_col <- intersect(c("p_val_adj", "p_val"), colnames(de))[[1]] %||% NULL
    fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(de))[[1]] %||% NULL
    if (is.null(p_col) || is.null(fc_col)) {
      stop("Could not identify p-value and log-fold-change columns in receiver DE results.", call. = FALSE)
    }
    geneset <- de$gene[de[[p_col]] <= 0.05 & de[[fc_col]] > 0]
  }
  geneset <- intersect(as.character(geneset), rownames(ligand_target_matrix))
  background_genes <- background_genes %||% .sn_expressed_genes(mat, receiver_cells, min_pct = expressed_pct)
  background_genes <- intersect(as.character(background_genes), rownames(ligand_target_matrix))

  expressed_sender <- .sn_expressed_genes(mat, sender_cells, min_pct = expressed_pct)
  expressed_receiver <- .sn_expressed_genes(mat, receiver_cells, min_pct = expressed_pct)
  lr_cols <- .sn_lr_columns(lr_network)
  lr_tbl <- tibble::as_tibble(lr_network)
  potential_ligands <- lr_tbl |>
    dplyr::filter(.data[[lr_cols$ligand]] %in% expressed_sender, .data[[lr_cols$receptor]] %in% expressed_receiver) |>
    dplyr::pull(dplyr::all_of(lr_cols$ligand)) |>
    unique()
  potential_ligands <- intersect(potential_ligands, colnames(ligand_target_matrix))
  if (length(potential_ligands) == 0L) {
    stop("No potential ligands remain after sender/receiver expression filtering.", call. = FALSE)
  }

  activities <- nichenetr::predict_ligand_activities(
    geneset = geneset,
    background_expressed_genes = background_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands,
    ...
  )
  activities <- tibble::as_tibble(activities) |>
    dplyr::arrange(dplyr::desc(.data$pearson)) |>
    dplyr::mutate(source = paste(sender, collapse = ","), target = paste(receiver, collapse = ",")) |>
    dplyr::relocate(dplyr::all_of(c("source", "target")))
  ligand_column <- .sn_communication_column(activities, c("test_ligand", "ligand"))
  ligand_receptors <- lr_tbl |>
    dplyr::filter(.data[[lr_cols$ligand]] %in% activities[[ligand_column]]) |>
    dplyr::transmute(
      test_ligand = as.character(.data[[lr_cols$ligand]]),
      receptor = as.character(.data[[lr_cols$receptor]])
    ) |>
    dplyr::distinct()
  if (!identical(ligand_column, "test_ligand")) names(ligand_receptors)[[1]] <- ligand_column
  activities <- dplyr::left_join(activities, ligand_receptors, by = ligand_column)
  if (!is.null(top_n)) {
    activities <- utils::head(activities, top_n)
  }
  ligand_targets <- .sn_nichenet_target_links(
    ligand_target_matrix,
    ligands = activities[[ligand_column]],
    geneset = geneset,
    top_n = max(20L, top_n %||% 50L)
  )
  list(
    table = activities,
    backend_result = activities,
    ligand_targets = ligand_targets,
    artifacts = list(
      geneset = geneset,
      background_genes = background_genes,
      potential_ligands = potential_ligands
    )
  )
}

.sn_run_liana <- function(object,
                          group_by,
                          assay = NULL,
                          layer = "data",
                          resource = NULL,
                          ...) {
  check_installed("liana", reason = "to run LIANA cell-cell communication inference.")
  liana_wrap <- get("liana_wrap", envir = asNamespace("liana"), inherits = FALSE)
  if (identical(tolower(resource %||% ""), "consensus")) resource <- "Consensus"
  input_name <- if ("sce" %in% names(formals(liana_wrap))) "sce" else "seurat_object"
  input_object <- object
  liana_assay <- assay
  if (identical(input_name, "sce")) {
    check_installed("SingleCellExperiment", reason = "to adapt Seurat v5 objects for current LIANA releases.")
    assay <- assay %||% Seurat::DefaultAssay(object)
    input_object <- suppressWarnings(Seurat::as.SingleCellExperiment(object, assay = assay))
    available_assays <- SummarizedExperiment::assayNames(input_object)
    liana_assay <- if (identical(layer, "counts") && "counts" %in% available_assays) {
      "counts"
    } else if ("logcounts" %in% available_assays) {
      "logcounts"
    } else {
      available_assays[[1]]
    }
  }
  input <- list(input_object)
  names(input) <- input_name
  args <- utils::modifyList(
    c(input, list(idents_col = group_by, assay = liana_assay, resource = resource, verbose = FALSE)),
    list(...),
    keep.null = TRUE
  )
  backend_warnings <- character()
  result <- withCallingHandlers(
    do.call(liana_wrap, args),
    warning = function(warning) {
      backend_warnings <<- c(backend_warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  table <- if (is.data.frame(result)) tibble::as_tibble(result) else tibble::as_tibble(result[[1]])
  list(
    table = table,
    backend_result = result,
    warnings = unique(backend_warnings),
    artifacts = list(liana = result)
  )
}

.sn_multinichenet_label_map <- function(values) {
  original <- unique(as.character(values))
  encoded <- make.unique(make.names(original))
  list(
    encoded = stats::setNames(encoded, original),
    decoded = stats::setNames(original, encoded)
  )
}

.sn_run_multinichenet <- function(object,
                                  group_by,
                                  sample_by,
                                  condition_by,
                                  assay,
                                  sender,
                                  receiver,
                                  condition_oi,
                                  condition_reference,
                                  contrast,
                                  ligand_target_matrix,
                                  lr_network,
                                  min_cells = 10,
                                  top_n = 250,
                                  ...) {
  check_installed("multinichenetr", reason = "to run the MultiNicheNet communication backend.")
  check_installed("SingleCellExperiment", reason = "to prepare MultiNicheNet input.")
  if (is_null(sample_by) || is_null(condition_by)) {
    stop("MultiNicheNet requires both `sample_by` and `condition_by`.", call. = FALSE)
  }
  if (is_null(ligand_target_matrix) || is_null(lr_network)) {
    stop("MultiNicheNet requires `ligand_target_matrix` and `lr_network`.", call. = FALSE)
  }
  metadata <- object[[]]
  .sn_validate_constant_within_sample(metadata, sample_col = sample_by, group_col = condition_by)
  conditions <- unique(as.character(metadata[[condition_by]]))
  contrast <- contrast %||% if (!is_null(condition_oi) && !is_null(condition_reference)) {
    c(condition_oi, condition_reference)
  } else if (length(conditions) == 2L) {
    conditions
  } else {
    NULL
  }
  if (is_null(contrast) || length(contrast) != 2L || any(!contrast %in% conditions)) {
    stop("MultiNicheNet requires a two-condition `contrast` present in `condition_by`.", call. = FALSE)
  }

  assay <- assay %||% Seurat::DefaultAssay(object)
  sce <- suppressWarnings(Seurat::as.SingleCellExperiment(object, assay = assay))
  group_map <- .sn_multinichenet_label_map(metadata[[group_by]])
  sample_map <- .sn_multinichenet_label_map(metadata[[sample_by]])
  condition_map <- .sn_multinichenet_label_map(metadata[[condition_by]])
  SummarizedExperiment::colData(sce)$.sn_celltype <- unname(group_map$encoded[as.character(metadata[[group_by]])])
  SummarizedExperiment::colData(sce)$.sn_sample <- unname(sample_map$encoded[as.character(metadata[[sample_by]])])
  SummarizedExperiment::colData(sce)$.sn_condition <- unname(condition_map$encoded[as.character(metadata[[condition_by]])])
  encoded_contrast <- unname(condition_map$encoded[contrast])
  contrast_name <- paste(encoded_contrast, collapse = "-")
  args <- utils::modifyList(list(
    sce = sce,
    celltype_id = ".sn_celltype",
    sample_id = ".sn_sample",
    group_id = ".sn_condition",
    batches = NA_character_,
    covariates = NA_character_,
    lr_network = tibble::as_tibble(lr_network),
    ligand_target_matrix = ligand_target_matrix,
    contrasts_oi = paste0("'", contrast_name, "'"),
    contrast_tbl = tibble::tibble(contrast = contrast_name, group = encoded_contrast[[1]]),
    senders_oi = if (is_null(sender)) NULL else unname(group_map$encoded[as.character(sender)]),
    receivers_oi = if (is_null(receiver)) NULL else unname(group_map$encoded[as.character(receiver)]),
    min_cells = min_cells,
    top_n_target = top_n,
    verbose = FALSE
  ), list(...), keep.null = TRUE)
  backend_warnings <- character()
  result <- withCallingHandlers(
    do.call(multinichenetr::multi_nichenet_analysis, args),
    warning = function(warning) {
      backend_warnings <<- c(backend_warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  table <- result$combined_prioritization_tables$group_prioritization_tbl %||%
    result$prioritization_tables$group_prioritization_tbl
  if (!is.data.frame(table)) stop("MultiNicheNet did not return a prioritization table.", call. = FALSE)
  table <- tibble::as_tibble(table)
  for (column in intersect(c("sender", "receiver"), names(table))) {
    table[[column]] <- unname(group_map$decoded[as.character(table[[column]])])
  }
  if ("group" %in% names(table)) table$group <- unname(condition_map$decoded[as.character(table$group)])
  ligand_targets <- result$lr_target_prior_cor %||% tibble::tibble()
  if (is.data.frame(ligand_targets)) ligand_targets <- tibble::as_tibble(ligand_targets)
  list(
    table = table,
    backend_result = result,
    ligand_targets = ligand_targets,
    warnings = unique(backend_warnings),
    artifacts = list(multinichenet = result, contrast = contrast)
  )
}

.sn_run_communication_backend <- function(object,
                                          method,
                                          group_by,
                                          assay,
                                          layer,
                                          species,
                                          sender,
                                          receiver,
                                          geneset,
                                          background_genes,
                                          condition_by,
                                          condition_oi,
                                          condition_reference,
                                          sample_by,
                                          contrast,
                                          ligand_target_matrix,
                                          lr_network,
                                          expressed_pct,
                                          top_n,
                                          cellchat_db,
                                          min_cells,
                                          population_size,
                                          raw_use,
                                          resource,
                                          controls) {
  switch(
    method,
    cellchat = do.call(.sn_run_cellchat, utils::modifyList(list(
      object = object, group_by = group_by, assay = assay, layer = layer,
      species = species, cellchat_db = cellchat_db, min_cells = min_cells,
      population_size = population_size, raw_use = raw_use
    ), controls, keep.null = TRUE)),
    nichenet = do.call(.sn_run_nichenetr, utils::modifyList(list(
      object = object, group_by = group_by, assay = assay, layer = layer,
      sender = sender, receiver = receiver, geneset = geneset,
      background_genes = background_genes, condition_col = condition_by,
      condition_oi = condition_oi, condition_reference = condition_reference,
      ligand_target_matrix = ligand_target_matrix, lr_network = lr_network,
      expressed_pct = expressed_pct, top_n = top_n
    ), controls, keep.null = TRUE)),
    liana = do.call(.sn_run_liana, utils::modifyList(list(
      object = object, group_by = group_by, assay = assay, layer = layer,
      resource = resource
    ), controls, keep.null = TRUE)),
    cellphonedb = {
      manifest <- do.call(sn_run_cellphonedb, utils::modifyList(list(
        object = object,
        assay = assay,
        layer = "counts",
        group_by = group_by,
        result_name = paste0("communication_", format(Sys.time(), "%Y%m%d%H%M%S")),
        return_object = FALSE,
        method_control = controls$method_control %||% list()
      ), controls[setdiff(names(controls), "method_control")], keep.null = TRUE))
      parsed <- .sn_read_cellphonedb_output(manifest$output_dir)
      list(table = parsed$table, backend_result = parsed$raw, artifacts = list(manifest = manifest))
    },
    multinichenet = do.call(.sn_run_multinichenet, utils::modifyList(list(
      object = object, group_by = group_by, sample_by = sample_by,
      condition_by = condition_by, assay = assay, sender = sender,
      receiver = receiver, condition_oi = condition_oi,
      condition_reference = condition_reference, contrast = contrast,
      ligand_target_matrix = ligand_target_matrix, lr_network = lr_network,
      min_cells = min_cells, top_n = top_n
    ), controls, keep.null = TRUE))
  )
}

#' Run cell-cell communication inference
#'
#' \code{sn_run_cell_communication()} wraps established communication
#' backends and stores a comparable ligand-receptor schema. Multiple backends
#' can be run together to calculate method concordance and a consensus rank.
#'
#' @param object A Seurat object.
#' @param method One or more of \code{"liana"}, \code{"cellchat"},
#'   \code{"cellphonedb"}, \code{"nichenet"}, or \code{"multinichenet"}.
#'   The legacy alias \code{"nichenetr"} is accepted.
#' @param group_by Metadata column defining sender/receiver cell groups.
#' @param assay,layer Assay and layer used to retrieve expression.
#' @param species Species passed to species-aware resources. Defaults to
#'   \code{sn_get_species(object)} when available.
#' @param sender,receiver Sender and receiver group labels required by
#'   \code{method = "nichenetr"}.
#' @param geneset Receiver gene set of interest for NicheNet. If omitted,
#'   receiver differential expression is computed from \code{condition_by},
#'   \code{condition_oi}, and \code{condition_reference}.
#' @param background_genes Background expressed genes for NicheNet.
#' @param condition_by Metadata column containing receiver conditions.
#' @param condition_oi,condition_reference Receiver condition contrast used to
#'   derive \code{geneset} for NicheNet.
#' @param ligand_target_matrix,lr_network NicheNet prior matrices/networks.
#'   These must be supplied for \code{method = "nichenet"}.
#' @param expressed_pct Minimum fraction of cells expressing a gene before it
#'   is considered expressed in NicheNet filtering.
#' @param top_n Number of top NicheNet ligands to keep.
#' @param cellchat_db Optional CellChat database object.
#' @param min_cells Minimum cells per group passed to CellChat filtering.
#' @param population_size,raw_use CellChat communication-probability controls.
#' @param resource Optional LIANA resource.
#' @param sample_by Metadata column defining biological samples. When supplied,
#'   ligand and receptor expression is aggregated within each sample before
#'   condition comparison.
#' @param consensus If \code{TRUE} and multiple methods are requested, use the
#'   cross-method consensus ranking as the primary result table.
#' @param contrast Optional length-two condition contrast, with case first and
#'   reference second.
#' @param backend_control Named list of method-specific argument lists.
#' @param store_name Name used under \code{object@misc$cell_communication_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object;
#'   otherwise return the stored-result list.
#' @param ... Backend-specific arguments.
#'
#' @return A Seurat object or stored-result list.
#' @export
sn_run_cell_communication <- function(object,
                                      method = c("liana", "cellchat", "cellphonedb", "nichenet", "multinichenet"),
                                      group_by,
                                      assay = NULL,
                                      layer = "data",
                                      species = NULL,
                                      sender = NULL,
                                      receiver = NULL,
                                      geneset = NULL,
                                      background_genes = NULL,
                                      condition_by = NULL,
                                      condition_oi = NULL,
                                      condition_reference = NULL,
                                      ligand_target_matrix = NULL,
                                      lr_network = NULL,
                                      expressed_pct = 0.1,
                                      top_n = 50,
                                      cellchat_db = NULL,
                                      min_cells = 10,
                                      population_size = FALSE,
                                      raw_use = TRUE,
                                      resource = NULL,
                                      sample_by = NULL,
                                      consensus = TRUE,
                                      contrast = NULL,
                                      backend_control = list(),
                                      store_name = "default",
                                      return_object = TRUE,
                                      ...) {
  .sn_validate_seurat_object(object)
  if (missing(method)) method <- "liana"
  requested_method <- tolower(as.character(method))
  method <- tolower(as.character(method))
  method[method == "nichenetr"] <- "nichenet"
  supported <- c("liana", "cellchat", "cellphonedb", "nichenet", "multinichenet")
  invalid <- setdiff(method, supported)
  if (length(invalid) > 0L) stop("Unsupported communication method(s): ", paste(invalid, collapse = ", "), ".", call. = FALSE)
  method <- unique(method)
  if (is.null(group_by) || !group_by %in% colnames(object[[]])) {
    stop("`group_by` must identify a metadata column in `object`.", call. = FALSE)
  }
  if (!is_null(sample_by) && !sample_by %in% colnames(object[[]])) stop("`sample_by` was not found in object metadata.", call. = FALSE)
  if (!is_null(condition_by) && !condition_by %in% colnames(object[[]])) stop("`condition_by` was not found in object metadata.", call. = FALSE)
  dots <- list(...)
  if (length(method) > 1L && length(dots) > 0L) {
    stop("For multiple communication methods, place method-specific arguments under `backend_control`.", call. = FALSE)
  }
  results <- lapply(method, function(current) {
    controls <- backend_control[[current]] %||% if (length(method) == 1L) dots else list()
    .sn_run_communication_backend(
      object, current, group_by, assay, layer, species, sender, receiver,
      geneset, background_genes, condition_by, condition_oi,
      condition_reference, sample_by, contrast, ligand_target_matrix, lr_network, expressed_pct,
      top_n, cellchat_db, min_cells, population_size, raw_use, resource,
      controls
    )
  })
  names(results) <- method
  standardized <- dplyr::bind_rows(lapply(method, function(current) {
    .sn_standardize_communication(
      results[[current]]$table,
      method = current,
      sender = paste(sender %||% NA_character_, collapse = ","),
      receiver = paste(receiver %||% NA_character_, collapse = ",")
    )
  }))
  consensus_table <- .sn_communication_consensus(standardized)
  primary <- if (isTRUE(consensus) && length(method) > 1L) consensus_table else standardized
  concordance <- .sn_communication_concordance(standardized)
  assay_resolved <- assay %||% Seurat::DefaultAssay(object)
  sample_evidence <- .sn_communication_sample_evidence(
    object, primary, group_by, sample_by, condition_by, assay_resolved, layer
  )
  comparison_contrast <- contrast %||% if (!is_null(condition_oi) && !is_null(condition_reference)) c(condition_oi, condition_reference) else NULL
  comparison <- .sn_compare_communication_samples(sample_evidence, contrast = comparison_contrast)
  ligand_targets <- dplyr::bind_rows(lapply(results, function(result) result$ligand_targets %||% tibble::tibble()))
  backend_warnings <- unique(unlist(lapply(results, function(result) result$warnings %||% character()), use.names = FALSE))
  artifacts <- lapply(results, `[[`, "artifacts")
  stored_method <- if (length(method) > 1L && isTRUE(consensus)) {
    "consensus"
  } else if (length(requested_method) == 1L && identical(requested_method, "nichenetr")) {
    "nichenetr"
  } else {
    paste(method, collapse = "+")
  }
  stored <- sn_store_cell_communication(
    object = object,
    result = primary,
    store_name = store_name,
    method = stored_method,
    backend = paste(method, collapse = "+"),
    group_by = group_by,
    sender = sender,
    receiver = receiver,
    species = species,
    artifacts = artifacts,
    raw_result = standardized,
    consensus_result = consensus_table,
    sample_evidence = sample_evidence,
    comparison = comparison,
    concordance = concordance,
    ligand_targets = ligand_targets,
    warnings = backend_warnings,
    sample_by = sample_by,
    condition_by = condition_by,
    return_object = return_object
  )
  stored
}

#' Store a cell-cell communication result on a Seurat object
#'
#' @param object A Seurat object.
#' @param result Communication result table.
#' @param store_name Name used under
#'   \code{object@misc$cell_communication_results}.
#' @param method User-facing communication method label.
#' @param backend Canonical backend identifier. This differs from \code{method}
#'   only when preserving a legacy alias such as \code{"nichenetr"}.
#' @param group_by Metadata column used for groups.
#' @param sender,receiver Optional sender/receiver labels.
#' @param species Optional species label.
#' @param artifacts Optional backend-specific artifacts.
#' @param raw_result,consensus_result,sample_evidence,comparison,concordance,ligand_targets
#'   Optional standardized secondary tables stored in the unified result.
#' @param warnings Character vector of backend warnings retained for audit.
#' @param sample_by,condition_by Optional sample and condition metadata columns.
#' @param return_object If \code{TRUE}, return the updated object.
#'
#' @return A Seurat object or stored-result list.
#' @export
sn_store_cell_communication <- function(object,
                                        result,
                                        store_name = "default",
                                        method = "cellchat",
                                        backend = method,
                                        group_by = NULL,
                                        sender = NULL,
                                        receiver = NULL,
                                        species = NULL,
                                        artifacts = NULL,
                                        raw_result = NULL,
                                        consensus_result = NULL,
                                        sample_evidence = NULL,
                                        comparison = NULL,
                                        concordance = NULL,
                                        ligand_targets = NULL,
                                        warnings = character(),
                                        sample_by = NULL,
                                        condition_by = NULL,
                                        return_object = TRUE) {
  .sn_validate_seurat_object(object)
  result <- if (all(c("source", "target", "ligand", "receptor", "score", "method") %in% names(result))) {
    tibble::as_tibble(result)
  } else {
    .sn_standardize_communication(result, method = method, sender = sender, receiver = receiver)
  }
  method <- paste(as.character(method), collapse = "+")
  backend <- paste(as.character(backend), collapse = "+")
  sample_count <- if (is.data.frame(sample_evidence) && "sample" %in% names(sample_evidence)) {
    length(unique(sample_evidence$sample))
  } else {
    0L
  }
  stored_result <- list(
    schema_version = "1.0",
    analysis_type = "cell_communication",
    name = store_name,
    backend = backend,
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = result,
    input = list(group_by = group_by, sample_by = sample_by, condition_by = condition_by, species = species),
    parameters = list(sender = sender, receiver = receiver),
    tables = list(
      primary = result,
      backend_raw = raw_result %||% result,
      consensus = consensus_result %||% tibble::tibble(),
      sample_evidence = sample_evidence %||% tibble::tibble(),
      condition_comparison = comparison %||% tibble::tibble(),
      method_concordance = concordance %||% tibble::tibble(),
      ligand_targets = ligand_targets %||% tibble::tibble()
    ),
    embeddings = list(),
    graphs = list(),
    models = list(artifact_names = names(artifacts %||% list())),
    diagnostics = list(
      interactions = nrow(result),
      methods = unique(result$method),
      samples = sample_count
    ),
    warnings = as.character(warnings),
    provenance = .sn_analysis_provenance(),
    analysis = "cell_communication",
    method = method,
    group_by = group_by,
    sender = sender,
    receiver = receiver,
    species = species,
    artifacts = artifacts
  )
  object <- .sn_store_misc_result(
    object = object,
    collection = "cell_communication_results",
    store_name = store_name,
    result = stored_result
  )
  if (isTRUE(return_object)) {
    return(.sn_log_seurat_command(object = object, name = "sn_store_cell_communication"))
  }
  stored_result
}

#' Retrieve a stored cell-cell communication result
#'
#' @param object A Seurat object.
#' @param communication_name Name of the stored result.
#' @param sources,targets Optional source/target labels to keep.
#' @param with_metadata If \code{TRUE}, return the full stored-result list.
#'
#' @return A tibble or stored-result list.
#' @export
sn_get_cell_communication_result <- function(object,
                                             communication_name = "default",
                                             sources = NULL,
                                             targets = NULL,
                                             with_metadata = FALSE) {
  .sn_validate_seurat_object(object)
  stored <- .sn_get_misc_result(
    object = object,
    collection = "cell_communication_results",
    store_name = communication_name
  )
  if (isTRUE(with_metadata)) {
    return(stored)
  }
  table <- tibble::as_tibble(stored$table)
  source_col <- .sn_communication_column(table, c("source", "sender", "source_cell", "cell_type1"))
  target_col <- .sn_communication_column(table, c("target", "receiver", "target_cell", "cell_type2"))
  if (!is.null(sources) && !is.null(source_col)) {
    table <- dplyr::filter(table, .data[[source_col]] %in% sources)
  }
  if (!is.null(targets) && !is.null(target_col)) {
    table <- dplyr::filter(table, .data[[target_col]] %in% targets)
  }
  table
}
