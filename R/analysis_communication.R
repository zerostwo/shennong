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
  ligand_col <- intersect(c("ligand", "from", "source"), colnames(lr_network))[[1]] %||% NULL
  receptor_col <- intersect(c("receptor", "to", "target"), colnames(lr_network))[[1]] %||% NULL
  if (is.null(ligand_col) || is.null(receptor_col)) {
    stop("`lr_network` must contain ligand/receptor columns such as `ligand`/`receptor` or `from`/`to`.", call. = FALSE)
  }
  list(ligand = ligand_col, receptor = receptor_col)
}

.sn_resource_species_name <- function(species) {
  if (identical(species, "human")) "human" else "mouse"
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
  keep <- !is.na(meta$shennong_group) & nzchar(meta$shennong_group)
  mat <- mat[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  cellchat <- CellChat::createCellChat(
    object = mat,
    meta = meta,
    group.by = "shennong_group",
    assay = assay
  )
  db <- cellchat_db %||% get(
    if (identical(species, "human")) "CellChatDB.human" else "CellChatDB.mouse",
    envir = asNamespace("CellChat")
  )
  cellchat@DB <- db
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(
    cellchat,
    raw.use = raw_use,
    population.size = population_size,
    ...
  )
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  table <- tibble::as_tibble(CellChat::subsetCommunication(cellchat))
  list(
    table = table,
    backend_result = cellchat,
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
  if (!is.null(top_n)) {
    activities <- utils::head(activities, top_n)
  }
  list(
    table = activities,
    backend_result = activities,
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
  result <- do.call(
    liana_wrap,
    c(
      list(
        seurat_object = object,
        idents_col = group_by,
        assay = assay,
        layer = layer,
        resource = resource
      ),
      list(...)
    )
  )
  table <- if (is.data.frame(result)) tibble::as_tibble(result) else tibble::as_tibble(result[[1]])
  list(table = table, backend_result = result, artifacts = list(liana = result))
}

#' Run cell-cell communication inference
#'
#' \code{sn_run_cell_communication()} wraps established communication
#' backends: CellChat for global interaction networks, NicheNet for
#' sender-to-receiver ligand activity, and LIANA for consensus ligand-receptor
#' scoring when the optional package is installed.
#'
#' @param object A Seurat object.
#' @param method One of \code{"cellchat"}, \code{"nichenetr"}, or
#'   \code{"liana"}.
#' @param group_by Metadata column defining sender/receiver cell groups.
#' @param assay,layer Assay and layer used to retrieve expression.
#' @param species Species passed to species-aware resources. Defaults to
#'   \code{sn_get_species(object)} when available.
#' @param sender,receiver Sender and receiver group labels required by
#'   \code{method = "nichenetr"}.
#' @param geneset Receiver gene set of interest for NicheNet. If omitted,
#'   receiver differential expression is computed from \code{condition_col},
#'   \code{condition_oi}, and \code{condition_reference}.
#' @param background_genes Background expressed genes for NicheNet.
#' @param condition_col,condition_oi,condition_reference Receiver condition
#'   contrast used to derive \code{geneset} for NicheNet.
#' @param ligand_target_matrix,lr_network NicheNet prior matrices/networks.
#'   These must be supplied for \code{method = "nichenetr"}.
#' @param expressed_pct Minimum fraction of cells expressing a gene before it
#'   is considered expressed in NicheNet filtering.
#' @param top_n Number of top NicheNet ligands to keep.
#' @param cellchat_db Optional CellChat database object.
#' @param min_cells Minimum cells per group passed to CellChat filtering.
#' @param population_size,raw_use CellChat communication-probability controls.
#' @param resource Optional LIANA resource.
#' @param store_name Name used under \code{object@misc$cell_communication_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object;
#'   otherwise return the stored-result list.
#' @param ... Backend-specific arguments.
#'
#' @return A Seurat object or stored-result list.
#' @export
sn_run_cell_communication <- function(object,
                                      method = c("cellchat", "nichenetr", "liana"),
                                      group_by,
                                      assay = NULL,
                                      layer = "data",
                                      species = NULL,
                                      sender = NULL,
                                      receiver = NULL,
                                      geneset = NULL,
                                      background_genes = NULL,
                                      condition_col = NULL,
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
                                      store_name = "default",
                                      return_object = TRUE,
                                      ...) {
  .sn_validate_seurat_object(object)
  method <- match.arg(method)
  if (is.null(group_by) || !group_by %in% colnames(object[[]])) {
    stop("`group_by` must identify a metadata column in `object`.", call. = FALSE)
  }

  result <- switch(
    method,
    cellchat = .sn_run_cellchat(
      object = object,
      group_by = group_by,
      assay = assay,
      layer = layer,
      species = species,
      cellchat_db = cellchat_db,
      min_cells = min_cells,
      population_size = population_size,
      raw_use = raw_use,
      ...
    ),
    nichenetr = .sn_run_nichenetr(
      object = object,
      group_by = group_by,
      assay = assay,
      layer = layer,
      sender = sender,
      receiver = receiver,
      geneset = geneset,
      background_genes = background_genes,
      condition_col = condition_col,
      condition_oi = condition_oi,
      condition_reference = condition_reference,
      ligand_target_matrix = ligand_target_matrix,
      lr_network = lr_network,
      expressed_pct = expressed_pct,
      top_n = top_n,
      ...
    ),
    liana = .sn_run_liana(
      object = object,
      group_by = group_by,
      assay = assay,
      layer = layer,
      resource = resource,
      ...
    )
  )

  stored <- sn_store_cell_communication(
    object = object,
    result = result$table,
    store_name = store_name,
    method = method,
    group_by = group_by,
    sender = sender,
    receiver = receiver,
    species = species,
    artifacts = result$artifacts,
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
#' @param method Communication backend.
#' @param group_by Metadata column used for groups.
#' @param sender,receiver Optional sender/receiver labels.
#' @param species Optional species label.
#' @param artifacts Optional backend-specific artifacts.
#' @param return_object If \code{TRUE}, return the updated object.
#'
#' @return A Seurat object or stored-result list.
#' @export
sn_store_cell_communication <- function(object,
                                        result,
                                        store_name = "default",
                                        method = "cellchat",
                                        group_by = NULL,
                                        sender = NULL,
                                        receiver = NULL,
                                        species = NULL,
                                        artifacts = NULL,
                                        return_object = TRUE) {
  .sn_validate_seurat_object(object)
  stored_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = tibble::as_tibble(result),
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
  source_col <- intersect(c("source", "sender", "source_cell", "cell_type1"), colnames(table))[[1]] %||% NULL
  target_col <- intersect(c("target", "receiver", "target_cell", "cell_type2"), colnames(table))[[1]] %||% NULL
  if (!is.null(sources) && !is.null(source_col)) {
    table <- dplyr::filter(table, .data[[source_col]] %in% sources)
  }
  if (!is.null(targets) && !is.null(target_col)) {
    table <- dplyr::filter(table, .data[[target_col]] %in% targets)
  }
  table
}

.sn_dorothea_network <- function(species, confidence_levels = c("A", "B", "C")) {
  check_installed("dorothea", reason = "to use the bundled DoRothEA regulons.")
  data_name <- if (identical(species, "human")) "dorothea_hs" else "dorothea_mm"
  env <- new.env(parent = emptyenv())
  utils::data(list = data_name, package = "dorothea", envir = env)
  network <- get(data_name, envir = env)
  if ("confidence" %in% colnames(network)) {
    network <- network[network$confidence %in% confidence_levels, , drop = FALSE]
  }
  network
}

.sn_progeny_network <- function(species, top = 500) {
  check_installed("progeny", reason = "to use the PROGENy pathway model.")
  model_fun <- get("getModel", envir = asNamespace("progeny"), inherits = FALSE)
  model_fun(organism = if (identical(species, "human")) "Human" else "Mouse", top = top)
}

.sn_normalize_regulatory_network <- function(network, method) {
  network <- tibble::as_tibble(network)
  if (identical(method, "dorothea")) {
    source_col <- intersect(c("tf", "source"), colnames(network))[[1]] %||% NULL
    target_col <- intersect(c("target", "gene"), colnames(network))[[1]] %||% NULL
    mor_col <- intersect(c("mor", "mor_tf", "weight"), colnames(network))[[1]] %||% NULL
  } else {
    source_col <- intersect(c("pathway", "source"), colnames(network))[[1]] %||% NULL
    target_col <- intersect(c("target", "gene"), colnames(network))[[1]] %||% NULL
    mor_col <- intersect(c("weight", "mor"), colnames(network))[[1]] %||% NULL
  }
  if (is.null(source_col) || is.null(target_col)) {
    stop("`network` must contain source and target columns.", call. = FALSE)
  }
  if (is.null(mor_col)) {
    network$mor <- 1
    mor_col <- "mor"
  }
  out <- data.frame(
    source = as.character(network[[source_col]]),
    target = as.character(network[[target_col]]),
    mor = as.numeric(network[[mor_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$source) & !is.na(out$target) & !is.na(out$mor), , drop = FALSE]
  unique(out)
}

#' Infer transcription-factor or pathway activity
#'
#' \code{sn_run_regulatory_activity()} runs fast footprint-style activity
#' inference using DoRothEA regulons or PROGENy pathway models through
#' \pkg{decoupleR}. DoRothEA returns transcription-factor activities; PROGENy
#' returns pathway activities.
#'
#' @param object A Seurat object.
#' @param method One of \code{"dorothea"} or \code{"progeny"}.
#' @param assay,layer Assay and layer used to retrieve expression.
#' @param group_by Optional metadata column. When supplied, expression is
#'   averaged by group before activity inference.
#' @param species Species used to choose built-in resources.
#' @param network Optional user-supplied regulatory network. If omitted,
#'   DoRothEA uses the \pkg{dorothea} package data and PROGENy uses
#'   \pkg{progeny} or \pkg{decoupleR} resources when installed.
#' @param confidence_levels DoRothEA confidence levels to keep.
#' @param progeny_top Number of target genes per PROGENy pathway.
#' @param minsize Minimum target-set size passed to \code{decoupleR::run_ulm()}.
#' @param store_name Name used under
#'   \code{object@misc$regulatory_activity_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#' @param ... Additional arguments passed to \code{decoupleR::run_ulm()}.
#'
#' @return A Seurat object or stored-result list.
#' @export
sn_run_regulatory_activity <- function(object,
                                       method = c("dorothea", "progeny"),
                                       assay = NULL,
                                       layer = "data",
                                       group_by = NULL,
                                       species = NULL,
                                       network = NULL,
                                       confidence_levels = c("A", "B", "C"),
                                       progeny_top = 500,
                                       minsize = 5L,
                                       store_name = "default",
                                       return_object = TRUE,
                                       ...) {
  .sn_validate_seurat_object(object)
  check_installed("decoupleR", reason = "to run fast DoRothEA/PROGENy activity inference.")
  method <- match.arg(method)
  species <- .sn_resolve_species(object, species)
  mat <- .sn_expression_matrix(object = object, assay = assay, layer = layer)
  if (!is.null(group_by)) {
    if (!group_by %in% colnames(object[[]])) {
      stop("`group_by` must identify a metadata column in `object`.", call. = FALSE)
    }
    mat <- .sn_group_average_matrix(mat, object[[group_by, drop = TRUE]][colnames(mat)])
  }

  network <- network %||% if (identical(method, "dorothea")) {
    .sn_dorothea_network(species, confidence_levels = confidence_levels)
  } else {
    .sn_progeny_network(species, top = progeny_top)
  }
  network <- .sn_normalize_regulatory_network(network, method = method)

  table <- decoupleR::run_ulm(
    mat = mat,
    network = network,
    .source = source,
    .target = target,
    .mor = mor,
    minsize = minsize,
    ...
  ) |>
    tibble::as_tibble()
  table$analysis_type <- if (identical(method, "dorothea")) "transcription_factor" else "pathway"

  sn_store_regulatory_activity(
    object = object,
    result = table,
    store_name = store_name,
    method = method,
    group_by = group_by,
    species = species,
    network = network,
    return_object = return_object
  )
}

#' Store regulatory activity results on a Seurat object
#'
#' @param object A Seurat object.
#' @param result Regulatory activity table.
#' @param store_name Name used under
#'   \code{object@misc$regulatory_activity_results}.
#' @param method Inference method.
#' @param group_by Optional grouping column.
#' @param species Optional species label.
#' @param network Optional regulatory network used for inference.
#' @param return_object If \code{TRUE}, return the updated object.
#'
#' @return A Seurat object or stored-result list.
#' @export
sn_store_regulatory_activity <- function(object,
                                         result,
                                         store_name = "default",
                                         method = "dorothea",
                                         group_by = NULL,
                                         species = NULL,
                                         network = NULL,
                                         return_object = TRUE) {
  .sn_validate_seurat_object(object)
  stored_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = tibble::as_tibble(result),
    analysis = "regulatory_activity",
    method = method,
    group_by = group_by,
    species = species,
    network = network
  )
  object <- .sn_store_misc_result(
    object = object,
    collection = "regulatory_activity_results",
    store_name = store_name,
    result = stored_result
  )
  if (isTRUE(return_object)) {
    return(.sn_log_seurat_command(object = object, name = "sn_store_regulatory_activity"))
  }
  stored_result
}

#' Retrieve stored regulatory activity results
#'
#' @param object A Seurat object.
#' @param activity_name Name of the stored result.
#' @param sources Optional TF or pathway names to keep.
#' @param conditions Optional cell or group names to keep.
#' @param with_metadata If \code{TRUE}, return the full stored-result list.
#'
#' @return A tibble or stored-result list.
#' @export
sn_get_regulatory_activity_result <- function(object,
                                              activity_name = "default",
                                              sources = NULL,
                                              conditions = NULL,
                                              with_metadata = FALSE) {
  .sn_validate_seurat_object(object)
  stored <- .sn_get_misc_result(
    object = object,
    collection = "regulatory_activity_results",
    store_name = activity_name
  )
  if (isTRUE(with_metadata)) {
    return(stored)
  }
  table <- tibble::as_tibble(stored$table)
  if (!is.null(sources) && "source" %in% colnames(table)) {
    table <- dplyr::filter(table, .data$source %in% sources)
  }
  if (!is.null(conditions) && "condition" %in% colnames(table)) {
    table <- dplyr::filter(table, .data$condition %in% conditions)
  }
  table
}
