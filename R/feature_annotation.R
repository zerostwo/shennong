.sn_feature_class_specs <- function() {
  list(
    transcription_factor = list(
      label = "Transcription factor",
      collection = "C5",
      subcollection = "GO:MF",
      terms = c("GOMF_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY")
    ),
    surface_membrane = list(
      label = "Cell-surface or plasma-membrane",
      collection = "C5",
      subcollection = "GO:CC",
      terms = c(
        "GOCC_CELL_SURFACE",
        "GOCC_EXTERNAL_SIDE_OF_PLASMA_MEMBRANE",
        "GOCC_EXTERNAL_SIDE_OF_APICAL_PLASMA_MEMBRANE",
        "GOCC_EXTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE",
        "GOCC_APICAL_PLASMA_MEMBRANE",
        "GOCC_APICOLATERAL_PLASMA_MEMBRANE",
        "GOCC_BASOLATERAL_PLASMA_MEMBRANE",
        "GOCC_LATERAL_PLASMA_MEMBRANE",
        "GOCC_PLASMA_MEMBRANE_RAFT"
      )
    ),
    cytokine = list(
      label = "Cytokine",
      collection = "C5",
      subcollection = "GO:MF",
      terms = c("GOMF_CYTOKINE_ACTIVITY")
    ),
    chemokine = list(
      label = "Chemokine",
      collection = "C5",
      subcollection = "GO:MF",
      terms = c("GOMF_CHEMOKINE_ACTIVITY")
    )
  )
}

.sn_resolve_feature_classes <- function(feature_classes = NULL) {
  specs <- .sn_feature_class_specs()
  feature_classes <- feature_classes %||% names(specs)
  feature_classes <- unique(as.character(feature_classes))
  unknown <- setdiff(feature_classes, names(specs))
  if (length(unknown) > 0) {
    stop(
      glue(
        "Unknown feature class: {paste(unknown, collapse = ', ')}. ",
        "Use one of: {paste(names(specs), collapse = ', ')}."
      ),
      call. = FALSE
    )
  }
  feature_classes
}

.sn_msigdbr_species <- function(species) {
  species <- rlang::arg_match(species, c("human", "mouse"))
  if (identical(species, "human")) "Homo sapiens" else "Mus musculus"
}

.sn_feature_class_resource_msigdbr <- function(species, feature_classes) {
  check_installed("msigdbr", reason = "to annotate DE genes by feature class with MSigDB GO resources.")
  specs <- .sn_feature_class_specs()[feature_classes]
  species_name <- .sn_msigdbr_species(species)
  keys <- unique(vapply(
    specs,
    function(spec) paste(spec$collection, spec$subcollection, sep = "\r"),
    character(1)
  ))
  msig_tables <- stats::setNames(lapply(keys, function(key) {
    parts <- strsplit(key, "\r", fixed = TRUE)[[1]]
    msigdbr::msigdbr(
      species = species_name,
      collection = parts[[1]],
      subcollection = parts[[2]]
    )
  }), keys)

  rows <- lapply(names(specs), function(feature_class) {
    spec <- specs[[feature_class]]
    key <- paste(spec$collection, spec$subcollection, sep = "\r")
    table <- msig_tables[[key]]
    table <- table[table$gs_name %in% spec$terms, , drop = FALSE]
    if (nrow(table) == 0L) {
      return(tibble::tibble())
    }
    tibble::tibble(
      species = species,
      gene = as.character(table$gene_symbol),
      feature_class = feature_class,
      feature_class_label = spec$label,
      feature_class_term = as.character(table$gs_name),
      feature_class_source = paste0("MSigDB ", table$db_version)
    ) |>
      dplyr::distinct()
  })

  dplyr::bind_rows(rows)
}

.sn_normalize_feature_class_resource <- function(resource,
                                                 species,
                                                 feature_classes) {
  if (is.null(resource)) {
    stop("`custom_resource` must be supplied when `resource = 'custom'`.", call. = FALSE)
  }
  if (!is.data.frame(resource)) {
    stop("`custom_resource` must be a data frame.", call. = FALSE)
  }
  if (!all(c("gene", "feature_class") %in% colnames(resource))) {
    stop("`custom_resource` must contain `gene` and `feature_class` columns.", call. = FALSE)
  }

  specs <- .sn_feature_class_specs()
  table <- tibble::as_tibble(resource)
  if ("species" %in% colnames(table)) {
    table <- dplyr::filter(table, is.na(.data$species) | .data$species == !!species)
  }
  table <- dplyr::filter(table, .data$feature_class %in% feature_classes)
  if (!"feature_class_label" %in% colnames(table)) {
    labels <- vapply(specs, function(x) x$label, character(1))
    table$feature_class_label <- unname(labels[as.character(table$feature_class)])
  }
  if (!"feature_class_term" %in% colnames(table)) {
    table$feature_class_term <- NA_character_
  }
  if (!"feature_class_source" %in% colnames(table)) {
    table$feature_class_source <- "custom"
  }
  table |>
    dplyr::mutate(
      species = species,
      gene = as.character(.data$gene),
      feature_class = as.character(.data$feature_class),
      feature_class_label = as.character(.data$feature_class_label),
      feature_class_term = as.character(.data$feature_class_term),
      feature_class_source = as.character(.data$feature_class_source)
    ) |>
    dplyr::select(
      "species",
      "gene",
      "feature_class",
      "feature_class_label",
      "feature_class_term",
      "feature_class_source"
    ) |>
    dplyr::distinct()
}

.sn_feature_class_resource <- function(species,
                                       feature_classes = NULL,
                                       resource = c("auto", "msigdbr", "custom"),
                                       custom_resource = NULL) {
  species <- rlang::arg_match(species, c("human", "mouse"))
  feature_classes <- .sn_resolve_feature_classes(feature_classes)
  resource <- match.arg(resource)
  if (identical(resource, "auto")) {
    resource <- if (!is.null(custom_resource)) "custom" else "msigdbr"
  }

  table <- if (identical(resource, "custom")) {
    .sn_normalize_feature_class_resource(
      resource = custom_resource,
      species = species,
      feature_classes = feature_classes
    )
  } else {
    .sn_feature_class_resource_msigdbr(species = species, feature_classes = feature_classes)
  }

  table |>
    dplyr::filter(.data$feature_class %in% feature_classes, !is.na(.data$gene), nzchar(.data$gene)) |>
    dplyr::mutate(.sn_gene_key = toupper(.data$gene)) |>
    dplyr::distinct()
}

.sn_collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & nzchar(x)]))
  if (length(x) == 0L) NA_character_ else paste(x, collapse = ";")
}

.sn_de_table_feature_annotation <- function(table,
                                            species,
                                            gene_col,
                                            feature_classes = NULL,
                                            resource = c("auto", "msigdbr", "custom"),
                                            custom_resource = NULL) {
  if (!is.data.frame(table)) {
    stop("`x` must be a data frame or a Seurat object with a stored DE result.", call. = FALSE)
  }
  if (!gene_col %in% colnames(table)) {
    stop(glue("Column '{gene_col}' was not found in the DE table."), call. = FALSE)
  }

  feature_classes <- .sn_resolve_feature_classes(feature_classes)
  class_resource <- .sn_feature_class_resource(
    species = species,
    feature_classes = feature_classes,
    resource = resource,
    custom_resource = custom_resource
  )
  input <- tibble::as_tibble(table)
  input$.sn_row_id <- seq_len(nrow(input))
  input$.sn_gene_key <- toupper(as.character(input[[gene_col]]))

  joined <- merge(
    input[, c(".sn_row_id", ".sn_gene_key"), drop = FALSE],
    class_resource,
    by = ".sn_gene_key",
    all.x = TRUE,
    sort = FALSE
  ) |>
    tibble::as_tibble()
  summary <- joined |>
    dplyr::group_by(.data$.sn_row_id) |>
    dplyr::summarise(
      feature_classes = .sn_collapse_unique(.data$feature_class),
      feature_class_labels = .sn_collapse_unique(.data$feature_class_label),
      feature_class_terms = .sn_collapse_unique(.data$feature_class_term),
      feature_class_sources = .sn_collapse_unique(.data$feature_class_source),
      .groups = "drop"
    )

  annotated <- dplyr::left_join(
    input,
    summary,
    by = ".sn_row_id"
  )
  matched_ids <- split(joined$.sn_row_id, joined$feature_class)
  for (feature_class in feature_classes) {
    col <- paste0("is_", feature_class)
    annotated[[col]] <- annotated$.sn_row_id %in% (matched_ids[[feature_class]] %||% integer(0))
  }

  annotated |>
    dplyr::select(-".sn_row_id", -".sn_gene_key")
}

#' Annotate DE or marker genes by feature class
#'
#' `sn_annotate_de_features()` flags genes from marker or differential
#' expression tables that encode transcription factors, cell-surface or
#' plasma-membrane proteins, cytokines, or chemokines. It can annotate a direct
#' data frame or a stored Shennong DE result under
#' `object@misc$de_results[[de_name]]`.
#'
#' By default the feature classes are derived from MSigDB GO sets through the
#' optional \pkg{msigdbr} package. For stricter project-specific definitions,
#' supply `resource = "custom"` with a table containing `gene` and
#' `feature_class` columns.
#'
#' @param x A DE/marker data frame, or a Seurat object containing a stored DE
#'   result.
#' @param de_name Stored DE result name when `x` is a Seurat object.
#' @param species One of `"human"` or `"mouse"`. When `x` is a Seurat object
#'   and `species` is `NULL`, Shennong tries `sn_get_species(x)`.
#' @param gene_col Column containing gene symbols in the DE table.
#' @param feature_classes Feature classes to annotate. Defaults to
#'   `"transcription_factor"`, `"surface_membrane"`, `"cytokine"`, and
#'   `"chemokine"`.
#' @param resource Annotation resource. `"msigdbr"` uses bundled MSigDB gene
#'   sets through \pkg{msigdbr}; `"custom"` uses `custom_resource`; `"auto"`
#'   chooses `"custom"` when `custom_resource` is supplied and otherwise
#'   `"msigdbr"`.
#' @param custom_resource Optional data frame with at least `gene` and
#'   `feature_class` columns. Optional columns are `species`,
#'   `feature_class_label`, `feature_class_term`, and
#'   `feature_class_source`.
#' @param store_name Optional stored DE result name for the annotated table
#'   when `x` is a Seurat object and `return_object = TRUE`. Defaults to
#'   `paste0(de_name, "_feature_classes")`.
#' @param return_object If `TRUE`, return the updated Seurat object with the
#'   annotated table stored under `object@misc$de_results[[store_name]]`.
#'   Otherwise return the annotated table.
#'
#' @return A tibble with feature-class columns, or an updated Seurat object.
#'
#' @examples
#' marker_tbl <- tibble::tibble(
#'   gene = c("TBX21", "CXCL10", "IL7R", "ACTB"),
#'   cluster = c("Tcell", "Myeloid", "Tcell", "Bcell"),
#'   avg_log2FC = c(2.1, 1.8, 1.2, 0.4)
#' )
#' custom_classes <- tibble::tibble(
#'   gene = c("TBX21", "CXCL10", "IL7R"),
#'   feature_class = c("transcription_factor", "chemokine", "surface_membrane")
#' )
#' sn_annotate_de_features(
#'   marker_tbl,
#'   species = "human",
#'   resource = "custom",
#'   custom_resource = custom_classes
#' )
#'
#' \dontrun{
#' obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'   store_name = "celltype_markers", return_object = TRUE
#' )
#' obj <- sn_annotate_de_features(obj, de_name = "celltype_markers")
#' sn_get_de_result(obj, de_name = "celltype_markers_feature_classes")
#' }
#' @export
sn_annotate_de_features <- function(x,
                                    de_name = "default",
                                    species = NULL,
                                    gene_col = "gene",
                                    feature_classes = NULL,
                                    resource = c("auto", "msigdbr", "custom"),
                                    custom_resource = NULL,
                                    store_name = NULL,
                                    return_object = inherits(x, "Seurat")) {
  resource <- match.arg(resource)
  if (inherits(x, "Seurat")) {
    species <- species %||% tryCatch(sn_get_species(x), error = function(...) NULL)
    species <- species %||% "human"
    species <- rlang::arg_match(species, c("human", "mouse"))
    stored <- .sn_get_misc_result(object = x, collection = "de_results", store_name = de_name)
    annotated <- .sn_de_table_feature_annotation(
      table = stored$table,
      species = species,
      gene_col = gene_col,
      feature_classes = feature_classes,
      resource = resource,
      custom_resource = custom_resource
    )
    if (!isTRUE(return_object)) {
      return(annotated)
    }
    store_name <- store_name %||% paste0(de_name, "_feature_classes")
    updated <- stored
    updated$table <- annotated
    updated$schema_version <- "1.1.0"
    updated$feature_annotation <- list(
      source_de_name = de_name,
      species = species,
      gene_col = gene_col,
      feature_classes = .sn_resolve_feature_classes(feature_classes),
      resource = if (identical(resource, "auto") && !is.null(custom_resource)) "custom" else resource
    )
    x <- .sn_store_misc_result(
      object = x,
      collection = "de_results",
      store_name = store_name,
      result = updated
    )
    return(.sn_log_seurat_command(object = x, name = "sn_annotate_de_features"))
  }

  species <- species %||% "human"
  species <- rlang::arg_match(species, c("human", "mouse"))
  .sn_de_table_feature_annotation(
    table = x,
    species = species,
    gene_col = gene_col,
    feature_classes = feature_classes,
    resource = resource,
    custom_resource = custom_resource
  )
}
