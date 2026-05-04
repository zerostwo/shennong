.sn_store_misc_result <- function(object, collection, store_name, result) {
  misc_data <- methods::slot(object, "misc")
  misc_data[[collection]] <- misc_data[[collection]] %||% list()
  misc_data[[collection]][[store_name]] <- result
  methods::slot(object, "misc") <- misc_data
  object
}

.sn_get_misc_result <- function(object, collection, store_name) {
  misc_data <- methods::slot(object, "misc")
  collection_data <- misc_data[[collection]] %||% list()
  if (!store_name %in% names(collection_data)) {
    stop(glue("No stored result named '{store_name}' was found in `object@misc${collection}`."))
  }
  collection_data[[store_name]]
}

.sn_resolve_misc_result_name <- function(object,
                                         collection,
                                         store_name = NULL,
                                         preferred_analysis = NULL,
                                         arg_name = "store_name") {
  misc_data <- methods::slot(object, "misc")
  collection_data <- misc_data[[collection]] %||% list()

  if (!is_null(store_name) && nzchar(store_name)) {
    if (!store_name %in% names(collection_data)) {
      stop(glue("No stored result named '{store_name}' was found in `object@misc${collection}`."))
    }
    return(store_name)
  }

  if (length(collection_data) == 0L) {
    stop(glue("No stored results were found in `object@misc${collection}`; please supply `{arg_name}`."), call. = FALSE)
  }

  available_names <- names(collection_data)
  latest_name <- function(candidates) {
    if (length(candidates) == 0L) {
      return(NULL)
    }
    created_at <- vapply(
      candidates,
      function(candidate) collection_data[[candidate]]$created_at %||% "",
      character(1)
    )
    candidates[[order(created_at, decreasing = TRUE, na.last = TRUE)[[1]]]]
  }

  preferred_names <- if (is_null(preferred_analysis)) {
    character()
  } else {
    available_names[vapply(
      available_names,
      function(candidate) identical(collection_data[[candidate]]$analysis %||% NULL, preferred_analysis),
      logical(1)
    )]
  }

  resolved_name <- if ("default" %in% available_names) {
    "default"
  } else if (length(available_names) == 1L) {
    available_names[[1]]
  } else {
    latest_name(preferred_names) %||% latest_name(available_names)
  }

  .sn_log_info("`{arg_name}` was not supplied; using stored result '{resolved_name}' from `object@misc${collection}`.")
  resolved_name
}

.sn_compact_collection_summary <- function(object, collection, type_label) {
  misc_data <- methods::slot(object, "misc")
  collection_data <- misc_data[[collection]] %||% list()
  if (length(collection_data) == 0) {
    return(tibble::tibble())
  }

  entries <- names(collection_data)
  tibble::tibble(
    collection = collection,
    type = type_label,
    name = entries,
    analysis = vapply(entries, function(entry) collection_data[[entry]]$analysis %||% NA_character_, character(1)),
    method = vapply(entries, function(entry) collection_data[[entry]]$method %||% NA_character_, character(1)),
    created_at = vapply(entries, function(entry) collection_data[[entry]]$created_at %||% NA_character_, character(1)),
    n_rows = vapply(entries, function(entry) nrow(collection_data[[entry]]$table %||% tibble::tibble()), integer(1)),
    source = vapply(entries, function(entry) collection_data[[entry]]$source_de_name %||% NA_character_, character(1))
  )
}

.sn_subset_ranked_table <- function(table,
                                    rank_col = NULL,
                                    group_col = NULL,
                                    top_n = NULL,
                                    direction = c("all", "up", "down"),
                                    groups = NULL) {
  direction <- match.arg(direction)
  table <- tibble::as_tibble(table)

  if (!is_null(groups) && !is_null(group_col) && group_col %in% colnames(table)) {
    table <- dplyr::filter(table, .data[[group_col]] %in% groups)
  }

  if (is_null(top_n) || is_null(rank_col) || !rank_col %in% colnames(table)) {
    return(table)
  }

  ranking <- table[[rank_col]]
  if (direction == "up") {
    table <- table[ranking > 0, , drop = FALSE]
    ordering <- ranking[ranking > 0]
  } else if (direction == "down") {
    table <- table[ranking < 0, , drop = FALSE]
    ordering <- abs(ranking[ranking < 0])
  } else {
    ordering <- abs(ranking)
  }

  table$..ranking_value <- ordering

  out <- if (!is_null(group_col) && group_col %in% colnames(table)) {
    table |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
      dplyr::slice_max(order_by = .data$..ranking_value, n = top_n, with_ties = FALSE) |>
      dplyr::ungroup()
  } else {
    table |>
      dplyr::slice_max(order_by = .data$..ranking_value, n = top_n, with_ties = FALSE)
  }

  dplyr::select(out, -dplyr::any_of("..ranking_value"))
}

#' List stored Shennong analysis and interpretation results on a Seurat object
#'
#' @param object A \code{Seurat} object.
#'
#' @return A tibble describing stored DE, enrichment, and interpretation
#'   results.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(10 * 12, lambda = 1), nrow = 10, ncol = 12)
#'   rownames(counts) <- c(
#'     "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
#'     "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
#'   )
#'   colnames(counts) <- paste0("cell", 1:12)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(
#'     obj,
#'     analysis = "markers",
#'     group_by = NULL,
#'     layer = "data",
#'     min_pct = 0,
#'     logfc_threshold = 0,
#'     return_object = TRUE,
#'     verbose = FALSE
#'   )
#'   sn_list_results(obj)
#' }
#' @export
sn_list_results <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  dplyr::bind_rows(
    .sn_compact_collection_summary(object, "de_results", "de"),
    .sn_compact_collection_summary(object, "cell_communication_results", "cell_communication"),
    .sn_compact_collection_summary(object, "deconvolution_results", "deconvolution"),
    .sn_compact_collection_summary(object, "milo_results", "milo"),
    .sn_compact_collection_summary(object, "regulatory_activity_results", "regulatory_activity"),
    .sn_compact_collection_summary(object, "enrichment_results", "enrichment"),
    .sn_compact_collection_summary(object, "interpretation_results", "interpretation")
  ) |>
    dplyr::arrange(.data$collection, .data$name)
}

#' Retrieve a stored DE result from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param de_name Name of the stored DE result.
#' @param top_n Optional number of rows to keep. When supplied together with a
#'   ranking column, results are reduced to the top rows overall or per group.
#' @param direction One of \code{"all"}, \code{"up"}, or \code{"down"}.
#' @param groups Optional subset of group labels to keep.
#' @param with_metadata If \code{TRUE}, return the full stored result list
#'   instead of just the result table.
#'
#' @return A tibble or stored-result list.
#'
#' @examples
#' \dontrun{
#' markers <- sn_get_de_result(seurat_obj, de_name = "cluster_markers", top_n = 5)
#' }
#' @export
sn_get_de_result <- function(object,
                             de_name = "default",
                             top_n = NULL,
                             direction = c("all", "up", "down"),
                             groups = NULL,
                             with_metadata = FALSE) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  stored <- .sn_get_misc_result(object = object, collection = "de_results", store_name = de_name)
  if (isTRUE(with_metadata)) {
    return(stored)
  }

  .sn_subset_ranked_table(
    table = stored$table,
    rank_col = stored$rank_col,
    group_col = stored$group_col,
    top_n = top_n,
    direction = direction,
    groups = groups
  )
}

#' Retrieve a stored enrichment result from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param enrichment_name Name of the stored enrichment result.
#' @param top_n Optional number of top terms to keep.
#' @param groups Optional subset of cluster/group labels when the stored table
#'   includes a \code{Cluster} column.
#' @param with_metadata If \code{TRUE}, return the full stored result list
#'   instead of just the term table.
#'
#' @return A tibble or stored-result list.
#'
#' @examples
#' \dontrun{
#' terms <- sn_get_enrichment_result(seurat_obj, enrichment_name = "cluster_gsea", top_n = 10)
#' }
#' @export
sn_get_enrichment_result <- function(object,
                                     enrichment_name = "default",
                                     top_n = NULL,
                                     groups = NULL,
                                     with_metadata = FALSE) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  stored <- .sn_get_misc_result(object = object, collection = "enrichment_results", store_name = enrichment_name)
  if (isTRUE(with_metadata)) {
    return(stored)
  }

  table <- tibble::as_tibble(stored$table)
  group_col <- c("Cluster", "cluster", ".sign")[
    c("Cluster", "cluster", ".sign") %in% colnames(table)
  ][1] %||% NULL
  rank_col <- c("NES", "Count", "GeneRatio", "p.adjust", "pvalue")[
    c("NES", "Count", "GeneRatio", "p.adjust", "pvalue") %in% colnames(table)
  ][1] %||% NULL

  if (!is_null(groups) && !is_null(group_col)) {
    table <- dplyr::filter(table, .data[[group_col]] %in% groups)
  }

  if (is_null(top_n) || is_null(rank_col)) {
    return(table)
  }

  if (rank_col %in% c("p.adjust", "pvalue")) {
    table <- table[order(table[[rank_col]], decreasing = FALSE), , drop = FALSE]
  } else {
    table <- table[order(abs(table[[rank_col]]), decreasing = TRUE), , drop = FALSE]
  }

  utils::head(table, top_n)
}

#' Retrieve a stored interpretation result from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param interpretation_name Name of the stored interpretation result.
#'
#' @return The stored interpretation-result list.
#'
#' @examples
#' \dontrun{
#' interpretation <- sn_get_interpretation_result(seurat_obj, "annotation_note")
#' }
#' @export
sn_get_interpretation_result <- function(object, interpretation_name = "default") {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  .sn_get_misc_result(
    object = object,
    collection = "interpretation_results",
    store_name = interpretation_name
  )
}

.sn_as_enrichment_table <- function(result) {
  if (is.data.frame(result)) {
    return(tibble::as_tibble(result))
  }

  tibble::as_tibble(as.data.frame(result))
}

.sn_compact_value <- function(x, max_items = 8) {
  if (length(x) == 0) {
    return("")
  }

  x <- unique(as.character(stats::na.omit(x)))
  if (length(x) == 0) {
    return("")
  }

  x <- utils::head(x, max_items)
  paste(x, collapse = ", ")
}

.sn_render_prompt_value <- function(x, max_rows = 8) {
  if (is.data.frame(x)) {
    rendered <- if (tibble::is_tibble(x)) {
      if (is.null(max_rows) || is.infinite(max_rows)) {
        utils::capture.output(print(x, n = Inf, width = Inf))
      } else {
        utils::capture.output(print(x, n = max_rows, width = Inf))
      }
    } else if (is.null(max_rows) || is.infinite(max_rows)) {
      utils::capture.output(print(x))
    } else {
      utils::capture.output(print(utils::head(x, max_rows)))
    }
    return(paste(rendered, collapse = "\n"))
  }

  if (is.list(x)) {
    rendered <- vapply(names(x), function(current_name) {
      paste0(current_name, ":\n", .sn_render_prompt_value(x[[current_name]], max_rows = max_rows))
    }, character(1))
    return(paste(rendered, collapse = "\n\n"))
  }

  paste(as.character(x), collapse = ", ")
}

.sn_extract_json_text <- function(text) {
  if (!is.character(text) || length(text) != 1 || is.na(text)) {
    return(NULL)
  }

  fenced <- stringr::str_match(text, "```json\\s*(\\{.*\\}|\\[.*\\])\\s*```")[, 2]
  if (!is.na(fenced) && nzchar(fenced)) {
    return(fenced)
  }

  json_like <- stringr::str_match(text, "(\\{[\\s\\S]*\\}|\\[[\\s\\S]*\\])")[, 2]
  if (!is.na(json_like) && nzchar(json_like)) {
    return(json_like)
  }

  NULL
}

.sn_http_content_type <- function(headers_text) {
  if (!is.character(headers_text) || length(headers_text) != 1L || is.na(headers_text)) {
    return(NA_character_)
  }
  match <- stringr::str_match(
    headers_text,
    "(?im)^content-type:\\s*([^;\\r\\n]+)"
  )[, 2]
  match %||% NA_character_
}

.sn_http_body_preview <- function(body_text, max_chars = 200) {
  if (!is.character(body_text) || length(body_text) != 1L || is.na(body_text) || !nzchar(body_text)) {
    return("")
  }
  preview <- stringr::str_squish(body_text)
  if (nchar(preview) > max_chars) {
    preview <- paste0(substr(preview, 1, max_chars), "...")
  }
  preview
}

.sn_is_probably_html <- function(body_text, content_type = NA_character_) {
  if (is.character(content_type) && length(content_type) == 1L && grepl("html", content_type, ignore.case = TRUE)) {
    return(TRUE)
  }
  if (!is.character(body_text) || length(body_text) != 1L || is.na(body_text)) {
    return(FALSE)
  }
  grepl("^\\s*<!DOCTYPE html|^\\s*<html\\b", body_text, ignore.case = TRUE)
}

.sn_parse_openai_http_response <- function(response,
                                           endpoint,
                                           wire_api,
                                           resolved_model) {
  status_code <- response$status_code %||% 200L
  headers_text <- tryCatch(rawToChar(response$headers), error = function(...) "")
  body_text <- tryCatch(rawToChar(response$content), error = function(...) "")
  content_type <- .sn_http_content_type(headers_text)
  preview <- .sn_http_body_preview(body_text)

  if (status_code >= 500L) {
    return(list(
      ok = FALSE,
      retryable = TRUE,
      message = glue(
        "OpenAI-style endpoint '{endpoint}' returned HTTP {status_code}. ",
        "Body preview: {preview}"
      )
    ))
  }

  if (.sn_is_probably_html(body_text = body_text, content_type = content_type)) {
    return(list(
      ok = FALSE,
      retryable = TRUE,
      message = glue(
        "OpenAI-style endpoint '{endpoint}' returned HTML instead of JSON",
        if (!is.na(content_type) && nzchar(content_type)) glue(" (content-type: {content_type})") else "",
        ". This usually means the base URL is wrong or the upstream gateway returned an error page. ",
        "Body preview: {preview}"
      )
    ))
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(body_text, simplifyVector = FALSE),
    error = identity
  )
  if (inherits(parsed, "error")) {
    return(list(
      ok = FALSE,
      retryable = TRUE,
      message = glue(
        "Failed to parse JSON from '{endpoint}'",
        if (!is.na(content_type) && nzchar(content_type)) glue(" (content-type: {content_type})") else "",
        ": {conditionMessage(parsed)}. Body preview: {preview}"
      )
    ))
  }

  text <- .sn_text_scalar(.sn_extract_openai_response_text(parsed = parsed, wire_api = wire_api))
  if (is.null(text)) {
    return(list(
      ok = FALSE,
      retryable = FALSE,
      message = glue("OpenAI-style response from '{endpoint}' did not contain extractable text.")
    ))
  }

  list(
    ok = TRUE,
    retryable = FALSE,
    parsed = parsed,
    text = text,
    model = parsed$model %||% resolved_model
  )
}

.sn_parse_annotation_response <- function(response) {
  if (is.list(response) && !is.null(response$structured)) {
    parsed <- response$structured
    annotations <- parsed$cluster_annotations %||% parsed$annotations %||% NULL
    if (is.null(annotations)) {
      return(NULL)
    }

    annotation_tbl <- tibble::as_tibble(annotations)
    if (!"cluster" %in% colnames(annotation_tbl)) {
      return(NULL)
    }

    for (col_name in intersect(c("risk_flags", "alternatives", "supporting_markers", "supporting_functions", "recommended_checks"), colnames(annotation_tbl))) {
      annotation_tbl[[col_name]] <- vapply(annotation_tbl[[col_name]], function(x) {
        if (is.null(x) || (length(x) == 1L && is.na(x))) {
          return(NA_character_)
        }
        if (is.list(x)) {
          x <- unlist(x, recursive = TRUE, use.names = FALSE)
        }
        x <- as.character(x)
        x <- x[!is.na(x) & nzchar(x)]
        if (length(x) == 0L) {
          return(NA_character_)
        }
        paste(x, collapse = "; ")
      }, character(1))
    }

    for (col_name in intersect(c("cluster", "primary_label", "broad_label", "confidence", "status", "note"), colnames(annotation_tbl))) {
      annotation_tbl[[col_name]] <- as.character(annotation_tbl[[col_name]])
    }

    return(list(
      table = annotation_tbl,
      narrative = parsed$narrative_summary %||% parsed$summary %||% NULL,
      raw = parsed
    ))
  }

  text <- response$text %||% NULL
  json_text <- .sn_extract_json_text(text)
  if (is.null(json_text)) {
    return(NULL)
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(json_text, simplifyDataFrame = TRUE),
    error = function(...) NULL
  )
  if (is.null(parsed)) {
    return(NULL)
  }

  annotations <- parsed$cluster_annotations %||% parsed$annotations %||% NULL
  if (is.null(annotations)) {
    return(NULL)
  }

  annotation_tbl <- tibble::as_tibble(annotations)
  if (!"cluster" %in% colnames(annotation_tbl)) {
    return(NULL)
  }

  if ("risk_flags" %in% colnames(annotation_tbl)) {
    annotation_tbl$risk_flags <- vapply(annotation_tbl$risk_flags, function(x) {
      if (is.list(x)) {
        paste(unlist(x), collapse = "; ")
      } else {
        paste(as.character(x), collapse = "; ")
      }
    }, character(1))
  }

  for (col_name in intersect(c("primary_label", "broad_label", "status"), colnames(annotation_tbl))) {
    annotation_tbl[[col_name]] <- as.character(annotation_tbl[[col_name]])
  }
  for (col_name in intersect(c("confidence", "alternatives", "supporting_markers", "supporting_functions", "note", "recommended_checks"), colnames(annotation_tbl))) {
    annotation_tbl[[col_name]] <- vapply(annotation_tbl[[col_name]], function(x) {
      if (is.list(x)) {
        paste(unlist(x), collapse = "; ")
      } else {
        paste(as.character(x), collapse = "; ")
      }
    }, character(1))
  }

  list(
    table = annotation_tbl,
    narrative = parsed$narrative_summary %||% parsed$summary %||% NULL,
    raw = parsed
  )
}

.sn_normalize_annotation_label <- function(label, style = c("title", "snake", "asis")) {
  style <- match.arg(style)
  label <- as.character(label %||% "")
  label <- stringr::str_replace_all(label, "[_\\-]+", " ")
  label <- stringr::str_squish(label)
  if (!nzchar(label)) {
    return(NA_character_)
  }

  if (identical(style, "asis")) {
    return(label)
  }
  if (identical(style, "snake")) {
    return(gsub("\\s+", "_", tolower(label)))
  }

  stringr::str_to_title(tolower(label))
}

.sn_normalize_annotation_table <- function(annotation_tbl, label_style = c("title", "snake", "asis")) {
  label_style <- match.arg(label_style)
  annotation_tbl <- tibble::as_tibble(annotation_tbl)
  for (col_name in intersect(c("primary_label", "broad_label", "alternatives"), colnames(annotation_tbl))) {
    values <- strsplit(as.character(annotation_tbl[[col_name]] %||% ""), ";", fixed = TRUE)
    annotation_tbl[[col_name]] <- vapply(values, function(x) {
      x <- trimws(x)
      x <- x[nzchar(x)]
      if (length(x) == 0) {
        return(NA_character_)
      }
      paste(vapply(x, .sn_normalize_annotation_label, character(1), style = label_style), collapse = "; ")
    }, character(1))
  }

  if ("status" %in% colnames(annotation_tbl)) {
    status_values <- tolower(annotation_tbl$status)
    status_values <- gsub("[^a-z0-9]+", "_", status_values)
    annotation_tbl$status <- status_values
  }

  annotation_tbl
}

.sn_apply_annotation_metadata <- function(object,
                                          annotation_tbl,
                                          cluster_col,
                                          metadata_prefix = "sn_annotation",
                                          metadata_fields = c("primary_label", "broad_label", "confidence", "status", "risk_flags")) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }
  if (!cluster_col %in% colnames(object[[]])) {
    stop(glue("Column '{cluster_col}' was not found in metadata."))
  }

  annotation_tbl <- tibble::as_tibble(annotation_tbl)
  if (!"cluster" %in% colnames(annotation_tbl)) {
    stop("`annotation_tbl` must contain a `cluster` column.")
  }

  meta <- object[[]] |>
    tibble::rownames_to_column("barcode")
  annotation_tbl$cluster <- as.character(annotation_tbl$cluster)
  meta[[cluster_col]] <- as.character(meta[[cluster_col]])
  merged <- dplyr::left_join(meta, annotation_tbl, by = stats::setNames("cluster", cluster_col))
  rownames(merged) <- merged$barcode

  field_map <- c(
    primary_label = paste0(metadata_prefix, "_label"),
    broad_label = paste0(metadata_prefix, "_broad_label"),
    confidence = paste0(metadata_prefix, "_confidence"),
    status = paste0(metadata_prefix, "_status"),
    alternatives = paste0(metadata_prefix, "_alternatives"),
    risk_flags = paste0(metadata_prefix, "_risk_flags"),
    supporting_markers = paste0(metadata_prefix, "_supporting_markers"),
    supporting_functions = paste0(metadata_prefix, "_supporting_functions"),
    note = paste0(metadata_prefix, "_note"),
    recommended_checks = paste0(metadata_prefix, "_recommended_checks")
  )
  metadata_fields <- unique(as.character(metadata_fields %||% character()))
  metadata_fields <- metadata_fields[nzchar(metadata_fields)]
  if (length(metadata_fields) > 0L) {
    unknown_fields <- setdiff(metadata_fields, names(field_map))
    if (length(unknown_fields) > 0L) {
      stop(glue("Unsupported `metadata_fields`: {paste(unknown_fields, collapse = ', ')}."), call. = FALSE)
    }
  }

  metadata_to_add <- data.frame(row.names = colnames(object))
  for (source_col in metadata_fields) {
    if (source_col %in% colnames(merged)) {
      metadata_to_add[[field_map[[source_col]]]] <- merged[colnames(object), source_col, drop = TRUE]
    }
  }
  metadata_to_add[[paste0(metadata_prefix, "_cluster")]] <- merged[colnames(object), "cluster", drop = TRUE]
  SeuratObject::AddMetaData(object, metadata = metadata_to_add)
}

.sn_compose_annotation_background <- function(background = NULL,
                                              label_candidates = NULL) {
  parts <- character()
  if (!is_null(background) && nzchar(background)) {
    parts <- c(parts, as.character(background))
  }
  label_candidates <- unique(as.character(label_candidates %||% character()))
  label_candidates <- label_candidates[nzchar(label_candidates)]
  if (length(label_candidates) > 0L) {
    parts <- c(
      parts,
      paste0(
        "Annotation priors / candidate labels: ",
        paste(label_candidates, collapse = ", "),
        ". Prefer these labels or closely related broad lineages/states when they are supported by the evidence. ",
        "If direct marker support is weak, do not force unrelated T-cell or NK-cell labels only because they are common blood populations."
      )
    )
  }
  prior_notes <- .sn_annotation_domain_priors(
    background = background,
    label_candidates = label_candidates
  )
  if (length(prior_notes) > 0L) {
    parts <- c(parts, prior_notes)
  }
  if (length(parts) == 0L) {
    return(NULL)
  }
  paste(parts, collapse = "\n")
}

.sn_annotation_domain_priors <- function(background = NULL,
                                         label_candidates = NULL) {
  context_text <- paste(
    as.character(background %||% ""),
    paste(as.character(label_candidates %||% character()), collapse = " "),
    collapse = " "
  )
  context_text <- tolower(context_text)
  notes <- character()

  if (grepl("blood|pbmc|peripheral blood", context_text) &&
      grepl("(^|[^a-z])ilc([0-9]|[^a-z]|$)|innate lymphoid", context_text)) {
    notes <- c(
      notes,
      paste(
        "Blood ILC prior:",
        "KIT-positive IL7R-positive KLRB1-positive helper-like clusters with weak or incomplete RORC-IL23R-NCR2 support",
        "are more plausibly KIT+ ILC / ILCP-like than mature ILC3.",
        "Reserve mature ILC3 labels for stronger type-3 programs."
      )
    )
  }
  if (grepl("tonsil|mucosa|intestinal|ileal|iel|nkp44", context_text)) {
    notes <- c(
      notes,
      paste(
        "Mucosal ILC3 prior:",
        "RORC-IL23R-AHR with KIT and especially NCR2/NKp44 support can justify ILC3-like labels in tonsil or intestinal samples,",
        "provided the cluster is not dominated by strong cytotoxic NK programs."
      )
    )
  }

  notes <- c(
    notes,
    paste(
      "Mixed-lineage prior:",
      "when strong TCR genes and strong NK cytotoxic genes coexist, consider T/NK mixed, contamination, or transitional states before forcing a pure T-cell or NK label."
    ),
    paste(
      "Contamination prior:",
      "when hemoglobin or erythroid genes dominate, label erythroid contamination explicitly rather than generic contamination."
    )
  )

  unique(notes[nzchar(notes)])
}

.sn_annotation_label_family <- function(label) {
  label <- tolower(as.character(label %||% ""))
  if (!nzchar(label)) {
    return(NA_character_)
  }
  if (grepl("erythroid|hemoglobin|epithelial|contamination", label)) {
    return("contamination")
  }
  if (grepl("innate lymphoid|(^|[^a-z])ilc([0-9]|[^a-z]|$)|kit\\+ ilc|ilc precursor", label)) {
    return("ilc")
  }
  if (grepl("mast|basophil", label)) {
    return("mast_basophil")
  }
  if (grepl("cytotoxic t/nk", label)) {
    return("cytotoxic_mixed")
  }
  if (grepl("(^|[^a-z])(t cell|t-cell|cd4|cd8)([^a-z]|$)", label)) {
    return("t")
  }
  if (grepl("(^|[^a-z])(nk|natural killer|cytotoxic lymphocyte)([^a-z]|$)", label)) {
    return("nk")
  }
  if (grepl("b cell|b-cell|plasma", label)) {
    return("b")
  }
  if (grepl("dendritic|apc|myeloid|monocyte", label)) {
    return("apc")
  }
  NA_character_
}

.sn_annotation_hint_family <- function(hint) {
  .sn_annotation_label_family(hint)
}

.sn_annotation_broad_label_from_family <- function(family) {
  switch(
    family %||% "",
    ilc = "Innate Lymphoid",
    t = "Lymphoid",
    nk = "Lymphoid",
    cytotoxic_mixed = "Lymphoid",
    b = "Lymphoid",
    apc = "Myeloid/APC",
    mast_basophil = "Mast/Basophil",
    contamination = "Contamination",
    NA_character_
  )
}

.sn_annotation_has_ilc_prior <- function(label_candidates = NULL,
                                         background = NULL) {
  candidate_text <- paste(as.character(label_candidates %||% character()), collapse = " ")
  background_text <- paste(as.character(background %||% character()), collapse = " ")
  grepl("(^|[^a-z])ilc([0-9]|[^a-z]|$)|innate lymphoid", paste(candidate_text, background_text), ignore.case = TRUE)
}

.sn_sort_annotation_subset <- function(tbl, clusters) {
  tbl <- tibble::as_tibble(tbl)
  if (nrow(tbl) == 0L || !"cluster" %in% colnames(tbl)) {
    return(tbl)
  }
  clusters <- unique(as.character(clusters))
  tbl$cluster <- as.character(tbl$cluster)
  tbl <- tbl[tbl$cluster %in% clusters, , drop = FALSE]
  tbl$..cluster_order <- match(tbl$cluster, clusters)
  tbl <- tbl[order(tbl$..cluster_order, tbl$cluster), , drop = FALSE]
  dplyr::select(tbl, -dplyr::any_of("..cluster_order"))
}

.sn_subset_annotation_evidence <- function(evidence,
                                           clusters) {
  clusters <- unique(as.character(clusters %||% character()))
  if (length(clusters) == 0L) {
    return(evidence)
  }

  out <- evidence
  for (name in c(
    "cluster_summary",
    "top_marker_table",
    "enrichment_summary",
    "qc_summary",
    "lineage_hints",
    "canonical_marker_snapshot",
    "geometry_summary"
  )) {
    value <- out[[name]]
    if (is.data.frame(value) && "cluster" %in% colnames(value)) {
      out[[name]] <- .sn_sort_annotation_subset(value, clusters = clusters)
    }
  }
  out$focus_clusters <- clusters
  out
}

.sn_prepare_annotation_compact_evidence <- function(evidence,
                                                    clusters = NULL,
                                                    include_marker_table = FALSE,
                                                    include_enrichment_table = FALSE,
                                                    include_qc = TRUE,
                                                    include_canonical_snapshot = FALSE,
                                                    include_geometry = TRUE) {
  compact <- .sn_subset_annotation_evidence(evidence, clusters = clusters)
  keep_cols <- c(
    "cluster", "n_cells", "fraction",
    grep("_distribution$", colnames(compact$cluster_summary %||% tibble::tibble()), value = TRUE),
    grep("(_predicted_labels|_majority_voting)$", colnames(compact$cluster_summary %||% tibble::tibble()), value = TRUE),
    "top_markers", "top_functions",
    "heuristic_hint", "heuristic_rationale", "heuristic_top_signatures",
    if (isTRUE(include_geometry)) "nearest_clusters",
    if (isTRUE(include_qc)) c(
      "median_nFeature_RNA", "median_nCount_RNA",
      grep("^median_percent\\.", colnames(compact$cluster_summary %||% tibble::tibble()), value = TRUE),
      "max_failed_qc_fraction", "doublet_fraction", "max_zero_count_fraction"
    )
  )
  keep_cols <- unique(stats::na.omit(unlist(keep_cols)))
  keep_cols <- intersect(keep_cols, colnames(compact$cluster_summary %||% tibble::tibble()))
  if (length(keep_cols) > 0L) {
    compact$cluster_summary <- dplyr::select(compact$cluster_summary, dplyr::all_of(keep_cols))
  }

  if (!isTRUE(include_marker_table)) {
    compact$top_marker_table <- tibble::tibble()
  }
  if (!isTRUE(include_enrichment_table)) {
    compact$enrichment_summary <- tibble::tibble()
  }
  if (!isTRUE(include_qc)) {
    compact$qc_summary <- tibble::tibble()
  }
  if (!isTRUE(include_canonical_snapshot)) {
    compact$canonical_marker_snapshot <- tibble::tibble()
  }
  if (!isTRUE(include_geometry)) {
    compact$geometry_summary <- tibble::tibble()
    compact$geometry_reduction <- NULL
  }
  compact
}

.sn_prepare_annotation_focus_clusters <- function(evidence,
                                                  broad_table = NULL,
                                                  label_candidates = NULL,
                                                  background = NULL,
                                                  max_clusters = 8L) {
  cluster_summary <- tibble::as_tibble(evidence$cluster_summary %||% tibble::tibble())
  if (nrow(cluster_summary) == 0L || !"cluster" %in% colnames(cluster_summary)) {
    return(character())
  }
  cluster_summary$cluster <- as.character(cluster_summary$cluster)

  broad_table <- tibble::as_tibble(broad_table %||% tibble::tibble())
  if (nrow(broad_table) > 0L && "cluster" %in% colnames(broad_table)) {
    broad_table$cluster <- as.character(broad_table$cluster)
    broad_keep <- intersect(c("cluster", "primary_label", "broad_label", "confidence", "status"), colnames(broad_table))
    broad_join <- broad_table[, broad_keep, drop = FALSE]
    rename_map <- c(
      cluster = "cluster",
      primary_label = "broad_pass_primary_label",
      broad_label = "broad_pass_broad_label",
      confidence = "broad_pass_confidence",
      status = "broad_pass_status"
    )
    colnames(broad_join) <- unname(rename_map[colnames(broad_join)])
    cluster_summary <- dplyr::left_join(
      cluster_summary,
      broad_join,
      by = "cluster"
    )
  }

  has_ilc_prior <- .sn_annotation_has_ilc_prior(
    label_candidates = label_candidates,
    background = background
  )
  predicted_family <- vapply(
    cluster_summary$broad_pass_primary_label %||% rep(NA_character_, nrow(cluster_summary)),
    .sn_annotation_label_family,
    character(1)
  )
  heuristic_family <- vapply(
    cluster_summary$heuristic_hint %||% rep(NA_character_, nrow(cluster_summary)),
    .sn_annotation_hint_family,
    character(1)
  )
  low_confidence <- tolower(cluster_summary$broad_pass_confidence %||% rep("", nrow(cluster_summary))) %in% c("", "low", "medium")
  non_confident <- tolower(cluster_summary$broad_pass_status %||% rep("", nrow(cluster_summary))) != "confident"
  ambiguous_family <- predicted_family %in% c("ilc", "t", "nk", "cytotoxic_mixed", "mast_basophil") |
    heuristic_family %in% c("ilc", "t", "nk", "cytotoxic_mixed", "mast_basophil")
  transitional_hint <- grepl(
    "mixed|transition|ilc3|erythroid",
    cluster_summary$heuristic_hint %||% "",
    ignore.case = TRUE
  )
  heuristic_conflict <- !is.na(predicted_family) & !is.na(heuristic_family) &
    predicted_family != heuristic_family &
    predicted_family %in% c("ilc", "t", "nk", "cytotoxic_mixed") &
    heuristic_family %in% c("ilc", "t", "nk", "cytotoxic_mixed")

  focus_flag <- low_confidence | non_confident | ambiguous_family | heuristic_conflict | transitional_hint
  if (isTRUE(has_ilc_prior)) {
    focus_flag <- focus_flag | grepl(
      "ILC|innate lymphoid|KIT\\+",
      cluster_summary$heuristic_hint %||% "",
      ignore.case = TRUE
    )
  }

  cluster_summary$..focus_priority <- 0
  cluster_summary$..focus_priority <- cluster_summary$..focus_priority + ifelse(heuristic_conflict, 4, 0)
  cluster_summary$..focus_priority <- cluster_summary$..focus_priority + ifelse(transitional_hint, 3, 0)
  cluster_summary$..focus_priority <- cluster_summary$..focus_priority + ifelse(low_confidence, 3, 0)
  cluster_summary$..focus_priority <- cluster_summary$..focus_priority + ifelse(non_confident, 2, 0)
  cluster_summary$..focus_priority <- cluster_summary$..focus_priority + ifelse(
    heuristic_family %in% c("ilc", "t", "nk", "cytotoxic_mixed", "mast_basophil"),
    1,
    0
  )
  if (isTRUE(has_ilc_prior)) {
    cluster_summary$..focus_priority <- cluster_summary$..focus_priority + ifelse(
      grepl("ILC|innate lymphoid|KIT\\+", cluster_summary$heuristic_hint %||% "", ignore.case = TRUE),
      2,
      0
    )
  }

  ranked_focus <- cluster_summary[focus_flag, , drop = FALSE]
  if (nrow(ranked_focus) == 0L) {
    return(character())
  }
  ranked_focus$..n_cells <- suppressWarnings(as.numeric(ranked_focus$n_cells %||% NA_real_))
  ranked_focus <- ranked_focus[order(-ranked_focus$..focus_priority, -ranked_focus$..n_cells, ranked_focus$cluster), , drop = FALSE]
  focus <- unique(as.character(stats::na.omit(ranked_focus$cluster)))
  if (length(focus) > max_clusters) {
    focus <- utils::head(focus, max_clusters)
  }
  focus
}

.sn_split_annotation_focus_batches <- function(clusters,
                                               batch_size = 4L) {
  clusters <- unique(as.character(clusters %||% character()))
  if (length(clusters) == 0L) {
    return(list())
  }
  split(clusters, ceiling(seq_along(clusters) / max(1L, as.integer(batch_size))))
}

.sn_make_annotation_evidence_tool <- function(evidence) {
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    return(NULL)
  }

  ellmer::tool(
    function(cluster_ids) {
      cluster_ids <- unique(as.character(cluster_ids %||% character()))
      if (length(cluster_ids) == 0L) {
        return(list(message = "No cluster IDs requested."))
      }
      .sn_prepare_annotation_compact_evidence(
        evidence = evidence,
        clusters = cluster_ids,
        include_marker_table = TRUE,
        include_enrichment_table = TRUE,
        include_qc = TRUE,
        include_canonical_snapshot = TRUE,
        include_geometry = !is.null(evidence$geometry_reduction)
      )
    },
    name = "lookup_cluster_annotation_evidence",
    description = paste(
      "Retrieve detailed Shennong evidence for one or more clusters, including",
      "specific markers, canonical marker snapshots, functional terms, QC, and geometry summaries."
    ),
    arguments = list(
      cluster_ids = ellmer::type_array(
        ellmer::type_string("Cluster identifier."),
        description = "One or more cluster IDs to inspect."
      )
    )
  )
}

.sn_build_annotation_stage_background <- function(background,
                                                  label_candidates = NULL,
                                                  stage = c("single_pass", "broad_pass", "focused_analysis", "focused_refinement"),
                                                  broad_pass_table = NULL,
                                                  focus_clusters = NULL) {
  stage <- match.arg(stage)
  stage_note <- switch(
    stage,
    single_pass = NULL,
    broad_pass = paste(
      "Annotation stage: broad pass.",
      "First assign the broadest defensible lineage or state for every cluster.",
      "Prefer lineage/state labels over over-specific subtypes when direct canonical support is limited."
    ),
    focused_analysis = paste(
      "Annotation stage: focused analysis.",
      "Compare the listed focus clusters and summarize the most discriminative lineage programs, conflicts, and caveats.",
      "Do not emit final JSON labels in this stage; produce only a concise comparison note."
    ),
    focused_refinement = paste(
      "Annotation stage: focused refinement.",
      "Reassess only the listed focus clusters by comparing them against each other.",
      "Use canonical marker differences to refine ILC-related, T-cell, NK-like, or mast/basophil-like identities.",
      "If subtype evidence is still weak, keep a conservative ILC-like or lineage-level label."
    )
  )

  parts <- c(background, stage_note)
  if (identical(stage, "focused_refinement") && length(focus_clusters %||% character()) > 0L) {
    parts <- c(parts, paste0("Focus clusters: ", paste(unique(as.character(focus_clusters)), collapse = ", "), "."))
  }
  if (identical(stage, "focused_refinement") && is.data.frame(broad_pass_table) && nrow(broad_pass_table) > 0L) {
    compact_broad <- broad_pass_table
    keep_cols <- intersect(c("cluster", "primary_label", "broad_label", "confidence", "status"), colnames(compact_broad))
    compact_broad <- compact_broad[, keep_cols, drop = FALSE]
    parts <- c(
      parts,
      paste(
        "Broad-pass annotations to refine:",
        .sn_render_prompt_value(compact_broad, max_rows = 200L)
      )
    )
  }

  .sn_compose_annotation_background(
    background = paste(parts[nzchar(parts %||% "")], collapse = "\n"),
    label_candidates = label_candidates
  )
}

.sn_reconcile_annotation_table <- function(annotation_tbl,
                                           evidence,
                                           label_candidates = NULL,
                                           background = NULL) {
  annotation_tbl <- tibble::as_tibble(annotation_tbl)
  if (nrow(annotation_tbl) == 0L || !"cluster" %in% colnames(annotation_tbl)) {
    return(annotation_tbl)
  }

  hints <- tibble::as_tibble(evidence$lineage_hints %||% tibble::tibble())
  if (nrow(hints) == 0L || !"cluster" %in% colnames(hints)) {
    return(annotation_tbl)
  }
  hints$cluster <- as.character(hints$cluster)
  annotation_tbl$cluster <- as.character(annotation_tbl$cluster)
  annotation_tbl <- dplyr::left_join(annotation_tbl, hints, by = "cluster")

  has_ilc_prior <- .sn_annotation_has_ilc_prior(
    label_candidates = label_candidates,
    background = background
  )

  for (i in seq_len(nrow(annotation_tbl))) {
    predicted_family <- .sn_annotation_label_family(annotation_tbl$primary_label[[i]])
    heuristic_hint <- annotation_tbl$heuristic_hint[[i]] %||% NA_character_
    heuristic_family <- .sn_annotation_hint_family(heuristic_hint)
    confidence <- tolower(annotation_tbl$confidence[[i]] %||% "")
    status <- tolower(annotation_tbl$status[[i]] %||% "")

    if (is.na(heuristic_family) || is.na(predicted_family)) {
      next
    }
    relaxed_conflict <- predicted_family != heuristic_family &&
      confidence %in% c("", "low", "medium") &&
      status %in% c("", "ambiguous", "possible_transition", "possible_contamination", "possible_low_quality")

    if (isTRUE(has_ilc_prior) &&
        heuristic_family == "ilc" &&
        predicted_family %in% c("t", "nk", "cytotoxic_mixed", "mast_basophil") &&
        relaxed_conflict) {
      annotation_tbl$status[[i]] <- "ambiguous"
      alternatives <- trimws(unlist(strsplit(annotation_tbl$alternatives[[i]] %||% "", ";", fixed = TRUE)))
      alternatives <- unique(c(alternatives, heuristic_hint))
      alternatives <- alternatives[nzchar(alternatives)]
      if (length(alternatives) > 0L) {
        annotation_tbl$alternatives[[i]] <- paste(alternatives, collapse = "; ")
      }
      note_text <- trimws(annotation_tbl$note[[i]] %||% "")
      reconcile_note <- paste(
        "Canonical lineage guardrail suggests an ILC-like alternative that should be reviewed against the model label.",
        "The automatic result was not force-overwritten."
      )
      annotation_tbl$note[[i]] <- paste(c(note_text, reconcile_note)[nzchar(c(note_text, reconcile_note))], collapse = " ")
    }
  }

  dplyr::select(annotation_tbl, -dplyr::any_of(c("heuristic_hint", "heuristic_rationale", "heuristic_top_signatures")))
}

.sn_default_llm_api_key <- function() {
  key <- Sys.getenv("OPENAI_API_KEY", unset = "")
  key
}

.sn_default_llm_base_url <- function() {
  base_url <- Sys.getenv("OPENAI_BASE_URL", unset = "")
  if (!nzchar(base_url)) {
    base_url <- "https://api.openai.com/v1"
  }
  base_url
}

.sn_default_llm_model <- function() {
  model <- Sys.getenv("OPENAI_MODEL", unset = "")
  if (!nzchar(model)) {
    model <- "gpt-4.1"
  }
  model
}

.sn_default_reasoning_effort <- function() {
  effort <- Sys.getenv("OPENAI_REASONING_EFFORT", unset = "")
  if (!nzchar(effort)) {
    return(NULL)
  }
  effort
}

.sn_guess_provider_type <- function(base_url) {
  if (grepl("^https://api\\.openai\\.com(/|$)", base_url)) {
    return("openai")
  }
  "openai_compatible"
}

.sn_is_retryable_llm_error <- function(message) {
  if (!is.character(message) || length(message) != 1L || is.na(message)) {
    return(FALSE)
  }
  grepl(
    "lexical error|invalid char in json text|<!DOCTYPE html>|<html|502|503|504|gateway|temporar|timeout|timed out|connection reset",
    message,
    ignore.case = TRUE
  )
}

.sn_harden_llm_error <- function(message) {
  if (!is.character(message) || length(message) != 1L || is.na(message) || !nzchar(message)) {
    return("LLM request failed.")
  }
  if (.sn_is_retryable_llm_error(message)) {
    return(paste(
      "LLM request failed because the upstream provider returned an HTML/error page or another non-JSON proxy response.",
      "This is common with third-party proxy endpoints such as Sub2API when authentication, gateway routing, or temporary upstream health is unstable.",
      "Check the endpoint and credentials, then retry."
    ))
  }
  message
}

.sn_get_default_ellmer_provider <- function() {
  sn_make_ellmer_provider(
    api_key = .sn_default_llm_api_key(),
    base_url = .sn_default_llm_base_url(),
    model = .sn_default_llm_model(),
    provider_type = .sn_guess_provider_type(.sn_default_llm_base_url()),
    reasoning_effort = .sn_default_reasoning_effort()
  )
}

.sn_annotation_json_schema_text <- function() {
  paste(
    "Return valid JSON with keys `cluster_annotations` and `narrative_summary`.",
    "`cluster_annotations` must be an array of objects with:",
    "`cluster`, `primary_label`, `broad_label`, `confidence`, and `status`.",
    "Optional fields may include `alternatives`, `supporting_markers`, `supporting_functions`, `risk_flags`, `note`, and `recommended_checks`.",
    "`confidence` should use one of: high, medium, low.",
    "`status` should use one of: confident, ambiguous, possible_contamination, possible_low_quality, possible_transition.",
    "`risk_flags` should be an array that may include: contamination, low_quality, transitional_state, doublet_like, proliferating_state."
  )
}

.sn_annotation_structured_type <- function() {
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    return(NULL)
  }

  ellmer::type_object(
    "Structured cluster annotation result.",
    cluster_annotations = ellmer::type_array(
      ellmer::type_object(
        "One annotation record per cluster.",
        cluster = ellmer::type_string("Cluster identifier."),
        primary_label = ellmer::type_string("Most plausible cell type or state label."),
        broad_label = ellmer::type_string("Broader lineage or family label."),
        confidence = ellmer::type_enum(
          c("high", "medium", "low"),
          "Calibrated confidence label."
        ),
        status = ellmer::type_enum(
          c("confident", "ambiguous", "possible_contamination", "possible_low_quality", "possible_transition"),
          "Annotation status."
        ),
        alternatives = ellmer::type_array(
          ellmer::type_string("Alternative plausible label."),
          description = "Alternative labels when ambiguity remains.",
          required = FALSE
        ),
        supporting_markers = ellmer::type_array(
          ellmer::type_string("Marker gene supporting the annotation."),
          description = "Directly supportive marker genes.",
          required = FALSE
        ),
        supporting_functions = ellmer::type_array(
          ellmer::type_string("Functional term supporting the annotation."),
          description = "Directly supportive pathways or biological processes.",
          required = FALSE
        ),
        risk_flags = ellmer::type_array(
          ellmer::type_enum(
            c("contamination", "low_quality", "transitional_state", "doublet_like", "proliferating_state")
          ),
          description = "Risk flags or caveats.",
          required = FALSE
        ),
        note = ellmer::type_string(
          "Brief cluster-level rationale.",
          required = FALSE
        ),
        recommended_checks = ellmer::type_array(
          ellmer::type_string("Suggested manual follow-up check."),
          description = "Suggested validation checks.",
          required = FALSE
        )
      ),
      description = "Array of cluster annotation records."
    ),
    narrative_summary = ellmer::type_string(
      "Optional concise narrative summary across all clusters.",
      required = FALSE
    )
  )
}

.sn_build_messages <- function(system_prompt, user_prompt) {
  list(
    list(role = "system", content = system_prompt),
    list(role = "user", content = user_prompt)
  )
}

.sn_render_table_markdown <- function(x, max_rows = 10) {
  if (!is.data.frame(x) || nrow(x) == 0) {
    return(NULL)
  }

  if (requireNamespace("knitr", quietly = TRUE)) {
    return(paste(utils::capture.output(knitr::kable(utils::head(x, max_rows), format = "pipe")), collapse = "\n"))
  }

  paste(utils::capture.output(print(utils::head(x, max_rows))), collapse = "\n")
}

.sn_render_evidence_markdown <- function(x, max_rows = 8) {
  if (is.data.frame(x)) {
    table_text <- .sn_render_table_markdown(x, max_rows = max_rows)
    return(table_text %||% "_No rows available._")
  }

  if (is.list(x)) {
    sections <- unlist(lapply(names(x), function(current_name) {
      current_value <- x[[current_name]]
      c(
        paste0("### ", current_name),
        .sn_render_evidence_markdown(current_value, max_rows = max_rows)
      )
    }), use.names = FALSE)
    return(paste(sections, collapse = "\n\n"))
  }

  value <- .sn_render_prompt_value(x, max_rows = max_rows)
  if (!nzchar(value)) {
    return("_No value available._")
  }
  paste0("```text\n", value, "\n```")
}

.sn_compose_background_parts <- function(...) {
  parts <- unlist(list(...), recursive = FALSE, use.names = FALSE)
  parts <- as.character(parts)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (length(parts) == 0L) {
    return(NULL)
  }
  paste(parts, collapse = "\n\n")
}

.sn_interpretation_task_instructions <- function(task, evidence) {
  template_name <- switch(
    task,
    annotation = "task_annotation.txt",
    de = if (identical(evidence$summary$analysis, "markers")) "task_de_markers.txt" else "task_de_contrast.txt",
    enrichment = "task_enrichment.txt",
    results = "task_results.txt",
    figure_legend = "task_figure_legend.txt",
    presentation_summary = "task_presentation_summary.txt"
  )

  paste(.sn_render_template(file.path("interpretation", template_name)), collapse = " ")
}

.sn_human_readable_prompt <- function(task,
                                      evidence,
                                      background = NULL,
                                      instruction = NULL) {
  sections <- c(
    paste0("# Shennong Interpretation Request: ", task),
    if (!is_null(background) && nzchar(background)) paste0("## Background\n", background),
    paste0("## Goal\n", instruction),
    "## Evidence Summary"
  )

  if (is.list(evidence)) {
    evidence_sections <- unlist(lapply(names(evidence), function(current_name) {
      current_value <- evidence[[current_name]]
      if (is.data.frame(current_value)) {
        table_text <- .sn_render_table_markdown(current_value)
        return(c(
          paste0("### ", current_name),
          if (!is_null(table_text)) table_text else "_No rows available._"
        ))
      }

      if (is.list(current_value) && !is.data.frame(current_value)) {
        return(c(
          paste0("### ", current_name),
          .sn_render_prompt_value(current_value)
        ))
      }

      c(
        paste0("### ", current_name),
        .sn_render_prompt_value(current_value)
      )
    }))
    sections <- c(sections, evidence_sections)
  }

  list(
    output_format = "human",
    task = task,
    text = paste(sections, collapse = "\n\n"),
    evidence = evidence
  )
}

.sn_interpret_elapsed_text <- function(started_at) {
  if (is.null(started_at)) {
    return(NA_character_)
  }
  elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
  paste0(formatC(elapsed, format = "f", digits = 1), "s")
}

.sn_interpret_progress_start <- function(task,
                                         enabled = interactive(),
                                         total_steps = 5L) {
  label <- switch(
    task,
    interpret_annotation = "sn_interpret_annotation",
    interpret_de = "sn_interpret_de",
    interpret_enrichment = "sn_interpret_enrichment",
    write_results = "sn_write_results",
    write_figure_legend = "sn_write_figure_legend",
    write_presentation_summary = "sn_write_presentation_summary",
    task
  )

  state <- list(
    enabled = isTRUE(enabled),
    task = task,
    label = label,
    started_at = Sys.time(),
    step = 0L,
    total_steps = as.integer(total_steps),
    cli_id = NULL
  )

  if (isTRUE(state$enabled) && requireNamespace("cli", quietly = TRUE)) {
    state$cli_id <- tryCatch(
      cli::cli_progress_bar(
        name = label,
        total = state$total_steps,
        clear = FALSE
      ),
      error = function(...) NULL
    )
  }

  .sn_log_info("[{label}] Starting interpretation workflow.")
  state
}

.sn_interpret_progress_step <- function(state, status) {
  if (is.null(state)) {
    return(state)
  }
  state$step <- min(state$total_steps, state$step + 1L)
  elapsed <- .sn_interpret_elapsed_text(state$started_at)
  if (!is.null(state$cli_id) && requireNamespace("cli", quietly = TRUE)) {
    update_result <- tryCatch(
      cli::cli_progress_update(
        id = state$cli_id,
        set = state$step,
        status = status
      ),
      error = identity
    )
    if (inherits(update_result, "error")) {
      state$cli_id <- NULL
    }
  }
  .sn_log_info("[{state$label}] Step {state$step}/{state$total_steps}: {status} (elapsed {elapsed}).")
  state
}

.sn_interpret_progress_done <- function(state, status = "Completed") {
  if (is.null(state)) {
    return(invisible(NULL))
  }
  elapsed <- .sn_interpret_elapsed_text(state$started_at)
  if (!is.null(state$cli_id) && requireNamespace("cli", quietly = TRUE)) {
    tryCatch(
      cli::cli_progress_done(id = state$cli_id),
      error = function(...) invisible(NULL)
    )
  }
  .sn_log_info("[{state$label}] {status} (total elapsed {elapsed}).")
  invisible(NULL)
}

.sn_finish_interpretation <- function(object,
                                      task,
                                      evidence,
                                      prompt,
                                      provider = NULL,
                                      model = NULL,
                                      store_name = "default",
                                      cluster_col = NULL,
                                      metadata_prefix = "sn_annotation",
                                      metadata_fields = c("primary_label", "broad_label", "confidence", "status", "risk_flags"),
                                      label_style = c("title", "snake", "asis"),
                                      apply_metadata = FALSE,
                                      return_prompt = FALSE,
                                      return_object = TRUE,
                                      progress_state = NULL,
                                      ...) {
  label_style <- match.arg(label_style)
  if (return_prompt || identical(prompt$output_format, "human")) {
    .sn_interpret_progress_done(progress_state, status = "Prompt prepared")
    return(prompt)
  }

  if (is_null(provider)) {
    provider <- tryCatch(
      .sn_get_default_ellmer_provider(),
      error = function(...) NULL
    )
  }
  if (is_null(provider)) {
    stop(
      "`provider` must be supplied unless `return_prompt = TRUE`, or ellmer-compatible credentials must be available via environment variables.",
      call. = FALSE
    )
  }

  progress_state <- .sn_interpret_progress_step(progress_state, "Waiting for LLM response")
  response <- sn_run_llm(
    messages = prompt$messages,
    provider = provider,
    model = model,
    structured_type = if (identical(task, "interpret_annotation")) .sn_annotation_structured_type() else NULL,
    ...
  )

  parsed_annotation <- NULL
  if (identical(task, "interpret_annotation")) {
    progress_state <- .sn_interpret_progress_step(progress_state, "Parsing structured annotation")
    parsed_annotation <- .sn_parse_annotation_response(response)
    if (!is.null(parsed_annotation)) {
      parsed_annotation$table <- .sn_normalize_annotation_table(
        parsed_annotation$table,
        label_style = label_style
      )
      if (isTRUE(apply_metadata) && !is.null(cluster_col)) {
        object <- .sn_apply_annotation_metadata(
          object = object,
          annotation_tbl = parsed_annotation$table,
          cluster_col = cluster_col,
          metadata_prefix = metadata_prefix,
          metadata_fields = metadata_fields
        )
      }
    }
  }

  progress_state <- .sn_interpret_progress_step(progress_state, "Storing interpretation result")
  interpretation_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    task = task,
    evidence = evidence,
    prompt = prompt,
    response = response,
    model_info = list(model = model),
    annotation_table = parsed_annotation$table %||% NULL,
    narrative_summary = parsed_annotation$narrative %||% NULL,
    metadata_prefix = if (isTRUE(apply_metadata) && identical(task, "interpret_annotation")) metadata_prefix else NULL
  )

  object <- .sn_store_misc_result(
    object = object,
    collection = "interpretation_results",
    store_name = store_name,
    result = interpretation_result
  )

  if (return_object) {
    .sn_interpret_progress_done(progress_state)
    return(.sn_log_seurat_command(object = object, name = paste0("sn_", task)))
  }

  .sn_interpret_progress_done(progress_state)
  response
}

.sn_merge_annotation_tables <- function(base_tbl, update_tbl) {
  base_tbl <- tibble::as_tibble(base_tbl)
  update_tbl <- tibble::as_tibble(update_tbl)
  if (nrow(base_tbl) == 0L) {
    return(update_tbl)
  }
  if (nrow(update_tbl) == 0L || !"cluster" %in% colnames(update_tbl)) {
    return(base_tbl)
  }

  base_tbl$cluster <- as.character(base_tbl$cluster)
  update_tbl$cluster <- as.character(update_tbl$cluster)
  extra_cols <- setdiff(colnames(update_tbl), colnames(base_tbl))
  for (col_name in extra_cols) {
    base_tbl[[col_name]] <- NA
  }
  base_tbl <- base_tbl[, union(colnames(base_tbl), colnames(update_tbl)), drop = FALSE]
  keep_base <- !base_tbl$cluster %in% update_tbl$cluster
  merged <- dplyr::bind_rows(base_tbl[keep_base, , drop = FALSE], update_tbl)
  .sn_sort_annotation_subset(merged, clusters = unique(c(base_tbl$cluster, update_tbl$cluster)))
}

.sn_finish_annotation_agentic <- function(object,
                                          evidence,
                                          broad_prompt,
                                          provider = NULL,
                                          model = NULL,
                                          store_name = "default",
                                          cluster_col = "seurat_clusters",
                                          metadata_prefix = "sn_annotation",
                                          metadata_fields = c("primary_label", "broad_label", "confidence", "status", "risk_flags"),
                                          label_style = c("title", "snake", "asis"),
                                          apply_metadata = FALSE,
                                          label_candidates = NULL,
                                          background = NULL,
                                          return_prompt = FALSE,
                                          return_object = TRUE,
                                          progress_state = NULL,
                                          ...) {
  label_style <- match.arg(label_style)
  if (return_prompt || identical(broad_prompt$output_format, "human")) {
    .sn_interpret_progress_done(progress_state, status = "Prompt prepared")
    return(list(
      output_format = broad_prompt$output_format,
      task = "annotation",
      annotation_mode = "agentic",
      broad_prompt = broad_prompt,
      focused_prompt = NULL,
      note = "Focused refinement prompt is generated after the broad-pass structured response is available.",
      evidence = evidence
    ))
  }

  if (is_null(provider)) {
    provider <- tryCatch(
      .sn_get_default_ellmer_provider(),
      error = function(...) NULL
    )
  }
  if (is_null(provider)) {
    stop(
      "`provider` must be supplied unless `return_prompt = TRUE`, or ellmer-compatible credentials must be available via environment variables.",
      call. = FALSE
    )
  }

  workflow <- list(
    annotation_mode = "agentic",
    broad_pass = list(
      prompt = broad_prompt,
      evidence = broad_prompt$evidence %||% NULL
    )
  )

  progress_state <- .sn_interpret_progress_step(progress_state, "Waiting for broad-pass annotation")
  broad_response <- sn_run_llm(
    messages = broad_prompt$messages,
    provider = provider,
    model = model,
    structured_type = .sn_annotation_structured_type(),
    ...
  )
  workflow$broad_pass$response <- broad_response

  progress_state <- .sn_interpret_progress_step(progress_state, "Parsing broad-pass annotation")
  broad_parsed <- .sn_parse_annotation_response(broad_response)
  if (is.null(broad_parsed) || is.null(broad_parsed$table) || nrow(broad_parsed$table) == 0L) {
    stop(
      "Agentic annotation requires a structured JSON response in the broad-pass stage, but the provider output could not be parsed.",
      call. = FALSE
    )
  }
  broad_parsed$table <- .sn_normalize_annotation_table(
    broad_parsed$table,
    label_style = label_style
  )
  workflow$broad_pass$annotation_table <- broad_parsed$table
  workflow$broad_pass$narrative_summary <- broad_parsed$narrative %||% NULL

  focus_clusters <- .sn_prepare_annotation_focus_clusters(
    evidence = evidence,
    broad_table = broad_parsed$table,
    label_candidates = label_candidates,
    background = background
  )
  workflow$focus_clusters <- focus_clusters

  final_table <- broad_parsed$table
  final_narrative <- broad_parsed$narrative %||% NULL
  focused_prompt <- NULL
  analysis_note <- NULL

  if (length(focus_clusters) > 0L) {
    progress_state <- .sn_interpret_progress_step(progress_state, "Running focused comparison analysis")
    analysis_tool <- .sn_make_annotation_evidence_tool(evidence)
    analysis_prompt <- sn_build_prompt(
      evidence = .sn_prepare_annotation_compact_evidence(
        evidence = evidence,
        clusters = focus_clusters,
        include_marker_table = FALSE,
        include_enrichment_table = FALSE,
        include_qc = TRUE,
        include_canonical_snapshot = TRUE,
        include_geometry = !is.null(evidence$geometry_reduction)
      ),
      task = "annotation",
      audience = "scientist",
      style = "focused cluster comparison note",
      background = paste(
        .sn_build_annotation_stage_background(
          background = background,
          label_candidates = label_candidates,
          stage = "focused_analysis",
          broad_pass_table = broad_parsed$table,
          focus_clusters = focus_clusters
        ),
        "Use the evidence lookup tool when you need richer per-cluster details before forming the refinement summary.",
        "Do not return JSON here; return only a concise comparison note that highlights the most discriminative lineage programs and likely failure modes.",
        sep = "\n"
      ),
      output_format = "llm",
      include_json_schema = FALSE
    )
    workflow$analysis_pass <- list(
      prompt = analysis_prompt,
      focus_clusters = focus_clusters
    )
    analysis_response <- tryCatch(
      sn_run_llm(
        messages = analysis_prompt$messages,
        provider = provider,
        model = model,
        tools = Filter(Negate(is.null), list(analysis_tool)),
        ...
      ),
      error = identity
    )
    workflow$analysis_pass$response <- analysis_response
    if (!inherits(analysis_response, "error")) {
      analysis_note <- analysis_response$text %||% NULL
      workflow$analysis_pass$text <- analysis_note
    } else {
      workflow$analysis_pass$error <- conditionMessage(analysis_response)
    }

    focus_batches <- .sn_split_annotation_focus_batches(focus_clusters, batch_size = 4L)
    workflow$focused_passes <- vector("list", length(focus_batches))

    for (batch_idx in seq_along(focus_batches)) {
      current_batch <- focus_batches[[batch_idx]]
      progress_state <- .sn_interpret_progress_step(
        progress_state,
        paste0("Building focused refinement prompt (batch ", batch_idx, "/", length(focus_batches), ")")
      )
      focused_evidence <- .sn_prepare_annotation_compact_evidence(
        evidence = evidence,
        clusters = current_batch,
        include_marker_table = TRUE,
        include_enrichment_table = TRUE,
        include_qc = TRUE,
        include_canonical_snapshot = TRUE,
        include_geometry = !is.null(evidence$geometry_reduction)
      )
      focused_prompt <- sn_build_prompt(
        evidence = focused_evidence,
        task = "annotation",
        audience = "scientist",
        style = "cell type annotation refinement",
        background = .sn_compose_background_parts(
          .sn_build_annotation_stage_background(
            background = background,
            label_candidates = label_candidates,
            stage = "focused_refinement",
            broad_pass_table = broad_parsed$table,
            focus_clusters = current_batch
          ),
          if (!is.null(analysis_note) && nzchar(analysis_note)) {
            paste0("Tool-assisted comparison note:\n", analysis_note)
          } else {
            NULL
          }
        ),
        output_format = "llm",
        include_json_schema = TRUE
      )
      workflow$focused_passes[[batch_idx]] <- list(
        batch = current_batch,
        prompt = focused_prompt,
        evidence = focused_evidence
      )

      progress_state <- .sn_interpret_progress_step(
        progress_state,
        paste0("Waiting for focused refinement annotation (batch ", batch_idx, "/", length(focus_batches), ")")
      )
      focused_response <- tryCatch(
        sn_run_llm(
          messages = focused_prompt$messages,
          provider = provider,
          model = model,
          structured_type = .sn_annotation_structured_type(),
          ...
        ),
        error = identity
      )
      workflow$focused_passes[[batch_idx]]$response <- focused_response
      if (inherits(focused_response, "error")) {
        workflow$focused_passes[[batch_idx]]$error <- conditionMessage(focused_response)
        next
      }

      focused_parsed <- .sn_parse_annotation_response(focused_response)
      if (!is.null(focused_parsed) && !is.null(focused_parsed$table) && nrow(focused_parsed$table) > 0L) {
        focused_parsed$table <- .sn_normalize_annotation_table(
          focused_parsed$table,
          label_style = label_style
        )
        workflow$focused_passes[[batch_idx]]$annotation_table <- focused_parsed$table
        workflow$focused_passes[[batch_idx]]$narrative_summary <- focused_parsed$narrative %||% NULL
        final_table <- .sn_merge_annotation_tables(
          base_tbl = final_table,
          update_tbl = focused_parsed$table
        )
        final_narrative <- focused_parsed$narrative %||% final_narrative
      }
    }

    if (length(workflow$focused_passes) == 1L) {
      workflow$focused_pass <- workflow$focused_passes[[1]]
    }
  }

  final_table <- .sn_reconcile_annotation_table(
    annotation_tbl = final_table,
    evidence = evidence,
    label_candidates = label_candidates,
    background = background
  )
  final_table <- .sn_normalize_annotation_table(final_table, label_style = label_style)

  if (isTRUE(apply_metadata) && !is.null(cluster_col)) {
    object <- .sn_apply_annotation_metadata(
      object = object,
      annotation_tbl = final_table,
      cluster_col = cluster_col,
      metadata_prefix = metadata_prefix,
      metadata_fields = metadata_fields
    )
  }

  progress_state <- .sn_interpret_progress_step(progress_state, "Storing interpretation result")
  interpretation_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    task = "interpret_annotation",
    annotation_mode = "agentic",
    evidence = evidence,
    prompt = broad_prompt,
    response = broad_response,
    model_info = list(model = model),
    annotation_table = final_table,
    narrative_summary = final_narrative,
    metadata_prefix = if (isTRUE(apply_metadata)) metadata_prefix else NULL,
    workflow = workflow
  )

  object <- .sn_store_misc_result(
    object = object,
    collection = "interpretation_results",
    store_name = store_name,
    result = interpretation_result
  )

  if (return_object) {
    .sn_interpret_progress_done(progress_state)
    return(.sn_log_seurat_command(object = object, name = "sn_interpret_annotation"))
  }

  .sn_interpret_progress_done(progress_state)
  list(
    text = final_narrative %||% broad_response$text %||% "",
    model = model,
    raw = workflow,
    annotation_table = final_table,
    broad_prompt = broad_prompt,
    focused_prompt = focused_prompt
  )
}

.sn_prepare_cluster_summary <- function(object, cluster_col = "seurat_clusters") {
  metadata <- object[[]]
  if (!cluster_col %in% colnames(metadata)) {
    stop(glue("Column '{cluster_col}' was not found in metadata."))
  }

  summary <- dplyr::count(metadata, .data[[cluster_col]], name = "n_cells") |>
    dplyr::rename(cluster = dplyr::all_of(cluster_col)) |>
    dplyr::mutate(fraction = .data$n_cells / sum(.data$n_cells))

  sample_columns <- intersect(c("sample", "study", "tissue"), colnames(metadata))
  if (length(sample_columns) > 0) {
    sample_summary <- lapply(sample_columns, function(current_col) {
      dplyr::count(metadata, .data[[cluster_col]], .data[[current_col]], name = "n_cells") |>
        dplyr::rename(cluster = dplyr::all_of(cluster_col), value = dplyr::all_of(current_col)) |>
        dplyr::group_by(.data$cluster) |>
        dplyr::summarise(
          !!paste0(current_col, "_distribution") := .sn_compact_value(
            paste0(.data$value, " (", .data$n_cells, ")")
          ),
          .groups = "drop"
        )
    })
    summary <- Reduce(function(x, y) dplyr::left_join(x, y, by = "cluster"), sample_summary, init = summary)
  }

  tibble::as_tibble(summary)
}

.sn_marker_logfc_col <- function(marker_table, rank_col = NULL) {
  candidates <- c(rank_col, "avg_log2FC", "avg_logFC", "log2FoldChange", "logFC")
  candidates <- unique(stats::na.omit(candidates))
  candidates[candidates %in% colnames(marker_table)][1] %||% NULL
}

.sn_specificity_frequency <- function(table,
                                      group_col,
                                      feature_col) {
  if (!group_col %in% colnames(table) || !feature_col %in% colnames(table) || nrow(table) == 0L) {
    return(stats::setNames(numeric(), character()))
  }
  freq <- table |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c(group_col, feature_col)))) |>
    dplyr::count(dplyr::across(dplyr::all_of(feature_col)), name = ".specificity_freq")
  stats::setNames(freq$.specificity_freq, freq[[feature_col]])
}

.sn_prepare_marker_candidates <- function(de_result,
                                          positive_only = TRUE) {
  marker_table <- de_result$table
  group_col <- de_result$group_col
  rank_col <- de_result$rank_col
  p_col <- de_result$p_col

  if (nrow(marker_table) == 0) {
    return(tibble::tibble())
  }
  if (is_null(group_col) || !group_col %in% colnames(marker_table)) {
    return(tibble::tibble())
  }
  if (is_null(rank_col) || !rank_col %in% colnames(marker_table)) {
    return(tibble::tibble())
  }

  working <- tibble::as_tibble(marker_table)
  if (!is_null(p_col) && p_col %in% colnames(working)) {
    working <- working[working[[p_col]] <= de_result$p_val_cutoff, , drop = FALSE]
  }
  if (nrow(working) == 0) {
    working <- tibble::as_tibble(marker_table)
  }

  logfc_col <- .sn_marker_logfc_col(working, rank_col = rank_col)
  if (isTRUE(positive_only) && !is.null(logfc_col)) {
    positive <- !is.na(working[[logfc_col]]) & working[[logfc_col]] > 0
    if (any(positive)) {
      working <- working[positive, , drop = FALSE]
    }
  }
  if (nrow(working) == 0) {
    return(tibble::tibble())
  }

  ranking <- abs(working[[rank_col]])
  if (identical(de_result$analysis, "markers") && !is.null(logfc_col)) {
    ranking <- pmax(working[[logfc_col]], 0)
  }
  working$..ranking_value <- ranking
  working$..specificity_freq <- .sn_specificity_frequency(
    table = working,
    group_col = group_col,
    feature_col = "gene"
  )[as.character(working$gene)] %||% rep(NA_real_, nrow(working))
  working$..specificity_freq[is.na(working$..specificity_freq)] <- Inf
  working
}

.sn_prepare_marker_summary <- function(de_result,
                                       n_markers = 10,
                                       selection = c("specific", "top")) {
  selection <- match.arg(selection)
  working <- .sn_prepare_marker_candidates(
    de_result = de_result,
    positive_only = identical(de_result$analysis, "markers")
  )
  group_col <- de_result$group_col

  if (nrow(working) == 0L) {
    return(tibble::tibble(cluster = character(), top_markers = character()))
  }

  top_markers <- if (identical(selection, "specific")) {
    working |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
      dplyr::arrange(.data$..specificity_freq, dplyr::desc(.data$..ranking_value), .by_group = TRUE) |>
      dplyr::slice_head(n = n_markers) |>
      dplyr::summarise(
        top_markers = .sn_compact_value(.data$gene, max_items = n_markers),
        .groups = "drop"
      )
  } else {
    working |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
      dplyr::slice_max(order_by = .data$..ranking_value, n = n_markers, with_ties = FALSE) |>
      dplyr::summarise(
        top_markers = .sn_compact_value(.data$gene, max_items = n_markers),
        .groups = "drop"
      )
  } |>
    dplyr::rename(cluster = dplyr::all_of(group_col))

  tibble::as_tibble(top_markers)
}

.sn_prepare_marker_table <- function(de_result,
                                     n_markers = 10,
                                     selection = c("specific", "top")) {
  selection <- match.arg(selection)
  marker_table <- .sn_prepare_marker_candidates(
    de_result = de_result,
    positive_only = identical(de_result$analysis, "markers")
  )
  group_col <- de_result$group_col

  if (nrow(marker_table) == 0 || is_null(group_col) || !group_col %in% colnames(marker_table)) {
    return(tibble::tibble())
  }

  ordered <- if (identical(selection, "specific")) {
    marker_table |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
      dplyr::arrange(.data$..specificity_freq, dplyr::desc(.data$..ranking_value), .by_group = TRUE) |>
      dplyr::slice_head(n = n_markers) |>
      dplyr::ungroup()
  } else {
    marker_table |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
      dplyr::slice_max(order_by = .data$..ranking_value, n = n_markers, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  dplyr::select(ordered, -dplyr::any_of(c("..ranking_value", "..specificity_freq")))
}

.sn_prepare_prediction_summary <- function(object, cluster_col = "seurat_clusters") {
  metadata <- object[[]]
  prediction_cols <- grep("(_predicted_labels|_majority_voting)$", colnames(metadata), value = TRUE)
  if (length(prediction_cols) == 0) {
    return(tibble::tibble())
  }

  summary_list <- lapply(prediction_cols, function(current_col) {
    dplyr::count(metadata, .data[[cluster_col]], .data[[current_col]], name = "n_cells") |>
      dplyr::rename(cluster = dplyr::all_of(cluster_col), label = dplyr::all_of(current_col)) |>
      dplyr::group_by(.data$cluster) |>
      dplyr::slice_max(order_by = .data$n_cells, n = 1, with_ties = FALSE) |>
      dplyr::transmute(
        cluster = .data$cluster,
        !!current_col := paste0(.data$label, " (", .data$n_cells, " cells)")
      )
  })

  Reduce(function(x, y) dplyr::full_join(x, y, by = "cluster"), summary_list)
}

.sn_prepare_annotation_qc_summary <- function(object, cluster_col = "seurat_clusters") {
  metadata <- object[[]]
  if (!cluster_col %in% colnames(metadata)) {
    stop(glue("Column '{cluster_col}' was not found in metadata."))
  }

  qc_cols <- c("percent.mt", "percent.ribo", "percent.hb")
  qc_cols <- qc_cols[qc_cols %in% colnames(metadata)]
  zero_cols <- grep("_zero_count$", colnames(metadata), value = TRUE)
  qc_flag_cols <- grep("_qc$", colnames(metadata), value = TRUE)
  doublet_col <- c("scDblFinder.class_corrected", "scDblFinder.class")[c("scDblFinder.class_corrected", "scDblFinder.class") %in% colnames(metadata)][1] %||% NULL

  rows <- lapply(split(metadata, metadata[[cluster_col]], drop = TRUE), function(current) {
    out <- data.frame(
      cluster = as.character(current[[cluster_col]][[1]]),
      n_cells = nrow(current),
      stringsAsFactors = FALSE
    )

    if ("nFeature_RNA" %in% colnames(current)) {
      out$median_nFeature_RNA <- stats::median(current$nFeature_RNA, na.rm = TRUE)
    }
    if ("nCount_RNA" %in% colnames(current)) {
      out$median_nCount_RNA <- stats::median(current$nCount_RNA, na.rm = TRUE)
    }
    for (col_name in qc_cols) {
      out[[paste0("median_", col_name)]] <- stats::median(current[[col_name]], na.rm = TRUE)
    }
    if (length(qc_flag_cols) > 0) {
      failed_flags <- vapply(qc_flag_cols, function(col_name) mean(current[[col_name]] == "Failed", na.rm = TRUE), numeric(1))
      out$max_failed_qc_fraction <- max(failed_flags, na.rm = TRUE)
    }
    if (!is.null(doublet_col)) {
      out$doublet_fraction <- mean(current[[doublet_col]] == "doublet", na.rm = TRUE)
    }
    if (length(zero_cols) > 0) {
      zero_rates <- vapply(zero_cols, function(col_name) mean(as.logical(current[[col_name]]), na.rm = TRUE), numeric(1))
      out$max_zero_count_fraction <- max(zero_rates, na.rm = TRUE)
    }
    out
  })

  tibble::as_tibble(do.call(rbind, rows))
}

.sn_prepare_cluster_enrichment_summary <- function(object,
                                                  enrichment_name,
                                                  n_terms = 5,
                                                  selection = c("specific", "top")) {
  selection <- match.arg(selection)
  if (is.null(enrichment_name)) {
    return(tibble::tibble())
  }

  stored <- .sn_get_misc_result(object = object, collection = "enrichment_results", store_name = enrichment_name)
  table <- tibble::as_tibble(stored$table)
  group_candidates <- c("Cluster", "cluster")[c("Cluster", "cluster") %in% colnames(table)]
  term_candidates <- c("Description", "ID")[c("Description", "ID") %in% colnames(table)]
  rank_candidates <- c("NES", "Count", "GeneRatio", "p.adjust", "pvalue")[c("NES", "Count", "GeneRatio", "p.adjust", "pvalue") %in% colnames(table)]
  group_col <- if (length(group_candidates) > 0) group_candidates[[1]] else NULL
  term_col <- if (length(term_candidates) > 0) term_candidates[[1]] else NULL
  rank_col <- if (length(rank_candidates) > 0) rank_candidates[[1]] else NULL

  if (is.null(group_col) || is.null(term_col)) {
    return(tibble::tibble())
  }

  working <- table
  if (!is.null(rank_col) && rank_col %in% c("NES")) {
    positive <- !is.na(working[[rank_col]]) & working[[rank_col]] > 0
    if (any(positive)) {
      working <- working[positive, , drop = FALSE]
    }
  }
  if (nrow(working) == 0L) {
    working <- table
  }
  working$..rank_value <- if (is.null(rank_col)) {
    0
  } else if (rank_col %in% c("p.adjust", "pvalue")) {
    -log10(pmax(working[[rank_col]], .Machine$double.xmin))
  } else {
    abs(working[[rank_col]])
  }
  working$..specificity_freq <- .sn_specificity_frequency(
    table = working,
    group_col = group_col,
    feature_col = term_col
  )[as.character(working[[term_col]])] %||% rep(NA_real_, nrow(working))
  working$..specificity_freq[is.na(working$..specificity_freq)] <- Inf

  ordered <- if (identical(selection, "specific")) {
    working |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
      dplyr::arrange(.data$..specificity_freq, dplyr::desc(.data$..rank_value), .by_group = TRUE) |>
      dplyr::slice_head(n = n_terms)
  } else if (is.null(rank_col)) {
    working |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
      dplyr::slice_head(n = n_terms)
  } else {
    working |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
      dplyr::slice_max(order_by = .data$..rank_value, n = n_terms, with_ties = FALSE)
  }

  ordered |>
    dplyr::summarise(
      top_functions = .sn_compact_value(.data[[term_col]], max_items = n_terms),
      .groups = "drop"
    ) |>
    dplyr::rename(cluster = dplyr::all_of(group_col))
}

.sn_prepare_cluster_geometry_summary <- function(object,
                                                 cluster_col = "seurat_clusters",
                                                 reduction = "umap",
                                                 n_neighbors = 3) {
  if (is.null(reduction) || !nzchar(reduction)) {
    return(tibble::tibble())
  }
  metadata <- object[[]]
  if (!cluster_col %in% colnames(metadata)) {
    stop(glue("Column '{cluster_col}' was not found in metadata."))
  }
  reduction_names <- names(object@reductions %||% list())
  if (!reduction %in% reduction_names) {
    return(tibble::tibble())
  }
  embeddings <- tryCatch(Seurat::Embeddings(object[[reduction]]), error = function(...) NULL)
  if (is.null(embeddings) || nrow(embeddings) == 0L) {
    return(tibble::tibble())
  }

  coords <- tibble::as_tibble(embeddings, rownames = "barcode")
  coords$cluster <- as.character(metadata[coords$barcode, cluster_col, drop = TRUE])
  coords <- coords[!is.na(coords$cluster), , drop = FALSE]
  coord_cols <- setdiff(colnames(coords), c("barcode", "cluster"))
  if (length(coord_cols) == 0L) {
    return(tibble::tibble())
  }

  centroids <- coords |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(coord_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  centroid_mat <- as.matrix(centroids[, coord_cols, drop = FALSE])
  rownames(centroid_mat) <- centroids$cluster
  if (nrow(centroid_mat) <= 1L) {
    return(tibble::tibble(cluster = rownames(centroid_mat), nearest_clusters = NA_character_))
  }

  dist_mat <- as.matrix(stats::dist(centroid_mat))
  nearest <- lapply(seq_len(nrow(dist_mat)), function(i) {
    distances <- dist_mat[i, ]
    distances <- distances[names(distances) != rownames(dist_mat)[i]]
    distances <- sort(distances, decreasing = FALSE, na.last = TRUE)
    if (length(distances) == 0L) {
      return(NA_character_)
    }
    keep <- utils::head(distances, n_neighbors)
    paste0(names(keep), " (distance=", formatC(unname(keep), format = "f", digits = 2), ")", collapse = "; ")
  })

  tibble::tibble(
    cluster = rownames(dist_mat),
    nearest_clusters = unlist(nearest, use.names = FALSE),
    geometry_reduction = reduction
  )
}

.sn_annotation_signature_catalog <- function(species = NULL) {
  list(
    ilc2 = c("IL7R", "KLRB1", "RORA", "GATA3", "PTGDR2", "HPGDS", "IL1RL1", "IL17RB"),
    kit_ilc = c("KIT", "IL1R1", "IL7R", "KLRB1", "RORA", "AHR", "TOX2"),
    ilc3 = c("RORC", "IL23R", "AHR", "KIT", "KLRB1", "NCR2"),
    t_cell = c("CD3D", "CD3E", "TRAC", "CD2", "CD5", "CD6", "LTB"),
    nk_cytotoxic = c("NKG7", "GNLY", "PRF1", "KLRD1", "NCR1", "FCGR3A", "TYROBP", "CST7"),
    b_cell = c("MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "HLA-DPA1"),
    dendritic_apc = c("FCER1A", "CD1C", "CLEC10A", "HLA-DRA", "CST3", "ZBTB46", "CLEC9A"),
    mast_basophil = c("KIT", "IL1R1", "MS4A2", "CCR3", "IL3RA", "GATA2", "HDC", "CLC", "TPSAB1", "TPSB2"),
    epithelial_contam = c("KRT8", "KRT18", "KRT19", "EPCAM", "CLDN3", "CEACAM6", "KRT17", "FXYD3", "TFF1"),
    erythroid_contam = c("HBB", "HBA1", "HBA2", "HBD", "HBM", "AHSP", "KLF1", "CA1", "CA2", "ALAS2", "BLVRB")
  )
}

.sn_select_annotation_hint <- function(score_row,
                                       marker_support = character(),
                                       feature_values = numeric()) {
  score_row <- sort(score_row, decreasing = TRUE, na.last = TRUE)
  top_name <- names(score_row)[1] %||% NA_character_
  top_value <- unname(score_row[1] %||% NA_real_)
  second_value <- unname(score_row[2] %||% NA_real_)
  margin <- top_value - (second_value %||% 0)
  support_text <- paste(marker_support, collapse = ", ")
  feature_values <- feature_values %||% numeric()
  fv <- function(name) {
    if (is.null(names(feature_values)) || !name %in% names(feature_values)) {
      return(0)
    }
    value <- feature_values[name][[1]]
    if (is.na(value)) 0 else value
  }
  t_cell_signal <- max(fv("CD3D"), fv("CD3E"), fv("TRAC"), fv("CD2"), fv("CD5"), fv("CD6"))
  nk_signal <- max(fv("NKG7"), fv("GNLY"), fv("PRF1"), fv("KLRD1"), fv("NCR1"))
  b_cell_signal <- max(fv("MS4A1"), fv("CD79A"), fv("CD79B"), fv("CD74"))
  epithelial_signal <- max(fv("KRT8"), fv("KRT18"), fv("KRT19"), fv("EPCAM"), fv("CLDN3"), fv("CEACAM6"))
  erythroid_signal <- max(fv("HBB"), fv("HBA1"), fv("HBA2"), fv("HBD"), fv("HBM"), fv("AHSP"), fv("KLF1"), fv("CA1"), fv("CA2"))
  ilc3_signal <- max(fv("RORC"), fv("IL23R"), fv("AHR"), fv("NCR2"))

  if (is.na(top_name) || is.na(top_value)) {
    return(list(
      hint = NA_character_,
      rationale = NA_character_
    ))
  }

  if (fv("HPGDS") >= 0.05 || fv("PTGDR2") >= 0.05 ||
      (fv("RORA") >= 8 && fv("IL7R") >= 5 && fv("GATA3") >= 4 && t_cell_signal < 1.5 && nk_signal < 1.5) ||
      length(intersect(marker_support, c("HPGDS", "PTGDR2", "RORA", "GATA3", "IL7R", "KLRB1"))) >= 2L) {
    return(list(
      hint = "ILC2-like",
      rationale = paste0("Canonical ILC2 program is strongest", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if ((fv("KIT") >= 1 || fv("IL1R1") >= 0.5) &&
      fv("IL7R") >= 4 &&
      t_cell_signal < 2 &&
      nk_signal < 2.5) {
    return(list(
      hint = "KIT+ ILC-like / ILC precursor-like",
      rationale = paste0("KIT-associated innate lymphoid program is strongest", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (epithelial_signal >= 2 && epithelial_signal > max(t_cell_signal, nk_signal)) {
    return(list(
      hint = "Epithelial contamination-like",
      rationale = paste0("Non-hematopoietic epithelial markers dominate", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (erythroid_signal >= 2 && erythroid_signal > max(t_cell_signal, nk_signal, epithelial_signal)) {
    return(list(
      hint = "Erythroid contamination-like",
      rationale = paste0("Hemoglobin/erythroid markers dominate", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (length(intersect(marker_support, c("MS4A1", "CD79A", "CD79B", "CD74"))) >= 2L) {
    return(list(
      hint = "B-cell-like",
      rationale = paste0("Canonical B-cell markers appear among the most specific markers", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (t_cell_signal >= 2.5 && nk_signal < 2) {
    return(list(
      hint = "T-cell-like",
      rationale = paste0("Canonical T-cell markers dominate", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (b_cell_signal >= 2 && t_cell_signal < 2 && nk_signal < 2) {
    return(list(
      hint = "B-cell-like",
      rationale = paste0("Canonical B-cell markers dominate", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (ilc3_signal >= 1.5 && nk_signal >= 2 && t_cell_signal < 2.5) {
    return(list(
      hint = "NK/ILC3 transitional-like",
      rationale = paste0("Type-3 ILC and cytotoxic NK-like programs coexist", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (nk_signal >= 2 && t_cell_signal >= 2) {
    return(list(
      hint = "T/NK mixed lymphoid-like",
      rationale = paste0("Both T-cell receptor and cytotoxic/NK programs are substantial", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (nk_signal >= 2.5 && t_cell_signal < 2) {
    return(list(
      hint = "NK/cytotoxic lymphocyte-like",
      rationale = paste0("Canonical cytotoxic/NK markers dominate", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (identical(top_name, "b_cell") && margin >= 0.1) {
    return(list(
      hint = "B-cell-like",
      rationale = paste0("Canonical B-cell program dominates", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (identical(top_name, "dendritic_apc") && margin >= 0.1) {
    return(list(
      hint = "Dendritic/APC-like",
      rationale = paste0("Antigen-presenting cell program dominates", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }
  if (identical(top_name, "ilc3") && top_value >= 0.75) {
    if (fv("RORC") >= 0.15 || fv("IL23R") >= 0.15 || fv("AHR") >= 1.5) {
      return(list(
        hint = "ILC3-like",
        rationale = paste0("Canonical ILC3-associated markers are enriched", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
      ))
    }
  }
  if (identical(top_name, "mast_basophil") && top_value >= 0.75) {
    return(list(
      hint = "Mast/basophil-like",
      rationale = paste0("KIT or basophil/mast-cell program is prominent", if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
    ))
  }

  list(
    hint = paste0(top_name, "-like"),
    rationale = paste0("Top canonical program is ", gsub("_", " ", top_name), if (nzchar(support_text)) paste0(" (", support_text, ")") else ".")
  )
}

.sn_prepare_annotation_lineage_hints <- function(object,
                                                 cluster_col = "seurat_clusters",
                                                 top_marker_table = tibble::tibble(),
                                                 species = NULL) {
  signatures <- .sn_annotation_signature_catalog(species = species)
  signature_genes <- unique(unlist(signatures, use.names = FALSE))
  available_genes <- intersect(signature_genes, rownames(object))
  if (length(available_genes) == 0L || !cluster_col %in% colnames(object[[]])) {
    return(tibble::tibble())
  }

  avg_expr <- Seurat::AverageExpression(
    object,
    assays = Seurat::DefaultAssay(object),
    group.by = cluster_col,
    features = available_genes,
    return.seurat = FALSE
  )[[Seurat::DefaultAssay(object)]]
  avg_expr <- as.matrix(avg_expr)
  if (nrow(avg_expr) == 0L || ncol(avg_expr) == 0L) {
    return(tibble::tibble())
  }
  cluster_levels <- unique(as.character(object[[cluster_col]][, 1]))
  prefixed_levels <- paste0("g", cluster_levels)
  if (all(colnames(avg_expr) %in% prefixed_levels)) {
    colnames(avg_expr) <- sub("^g", "", colnames(avg_expr))
  }

  score_mat <- sapply(signatures, function(genes) {
    present <- intersect(genes, rownames(avg_expr))
    if (length(present) == 0L) {
      rep(NA_real_, ncol(avg_expr))
    } else {
      colMeans(avg_expr[present, , drop = FALSE], na.rm = TRUE)
    }
  })
  if (is.null(dim(score_mat))) {
    score_mat <- matrix(score_mat, ncol = 1)
  }
  rownames(score_mat) <- colnames(avg_expr)
  score_tbl <- as.data.frame(score_mat, check.names = FALSE) |>
    tibble::rownames_to_column("cluster")

  top_marker_table <- tibble::as_tibble(top_marker_table)
  support_lookup <- if (nrow(top_marker_table) > 0 && all(c("cluster", "gene") %in% colnames(top_marker_table))) {
    split(as.character(top_marker_table$gene), as.character(top_marker_table$cluster))
  } else {
    list()
  }

  hint_rows <- lapply(seq_len(nrow(score_tbl)), function(i) {
    cluster <- as.character(score_tbl$cluster[[i]])
    row_scores <- unlist(score_tbl[i, setdiff(colnames(score_tbl), "cluster"), drop = TRUE], use.names = TRUE)
    expr_values <- avg_expr[, cluster, drop = TRUE]
    marker_support <- support_lookup[[cluster]] %||% character()
    selected <- .sn_select_annotation_hint(
      score_row = row_scores,
      marker_support = marker_support,
      feature_values = expr_values
    )
    top_scores <- sort(row_scores, decreasing = TRUE, na.last = TRUE)
    top_scores <- utils::head(top_scores, 3)
    tibble::tibble(
      cluster = cluster,
      heuristic_hint = selected$hint,
      heuristic_rationale = selected$rationale,
      heuristic_top_signatures = paste0(
        names(top_scores),
        "=",
        formatC(unname(top_scores), format = "f", digits = 2),
        collapse = "; "
      )
    )
  })

  dplyr::bind_rows(hint_rows)
}

.sn_prepare_annotation_canonical_snapshot <- function(object,
                                                      cluster_col = "seurat_clusters",
                                                      species = NULL) {
  signatures <- .sn_annotation_signature_catalog(species = species)
  snapshot_genes <- unique(unlist(signatures[c(
    "ilc2", "kit_ilc", "ilc3", "t_cell", "b_cell", "dendritic_apc",
    "nk_cytotoxic", "mast_basophil", "epithelial_contam", "erythroid_contam"
  )], use.names = FALSE))
  available_genes <- intersect(snapshot_genes, rownames(object))
  if (length(available_genes) == 0L || !cluster_col %in% colnames(object[[]])) {
    return(tibble::tibble())
  }

  avg_expr <- Seurat::AverageExpression(
    object,
    assays = Seurat::DefaultAssay(object),
    group.by = cluster_col,
    features = available_genes,
    return.seurat = FALSE
  )[[Seurat::DefaultAssay(object)]]
  avg_expr <- as.matrix(avg_expr)
  if (nrow(avg_expr) == 0L || ncol(avg_expr) == 0L) {
    return(tibble::tibble())
  }

  cluster_levels <- unique(as.character(object[[cluster_col]][, 1]))
  prefixed_levels <- paste0("g", cluster_levels)
  if (all(colnames(avg_expr) %in% prefixed_levels)) {
    colnames(avg_expr) <- sub("^g", "", colnames(avg_expr))
  }

  snapshot <- as.data.frame(t(avg_expr), check.names = FALSE)
  snapshot$cluster <- rownames(snapshot)
  snapshot <- tibble::as_tibble(snapshot[, c("cluster", setdiff(colnames(snapshot), "cluster")), drop = FALSE])
  numeric_cols <- setdiff(colnames(snapshot), "cluster")
  snapshot[numeric_cols] <- lapply(snapshot[numeric_cols], function(x) round(as.numeric(x), 3))
  snapshot
}

#' Store an enrichment result on a Seurat object
#'
#' This helper stores enrichment output inside
#' `object@misc$enrichment_results[[store_name]]` so interpretation and writing
#' helpers can reuse it later.
#'
#' @param object A \code{Seurat} object.
#' @param result An enrichment result object or data frame coercible with
#'   \code{as.data.frame()}.
#' @param store_name Name used under \code{object@misc$enrichment_results}.
#' @param analysis One of \code{"ora"} or \code{"gsea"}.
#' @param database Database used for enrichment, for example \code{"GOBP"}.
#' @param species Species label used in the enrichment run.
#' @param source_de_name Optional stored DE result name that produced the input
#'   ranked gene list or gene set.
#' @param gene_col Column containing gene symbols when the enrichment input came
#'   from a data frame.
#' @param score_col Column containing ranking scores for GSEA inputs.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#'
#' @return A \code{Seurat} object or a stored-result list.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(10 * 12, lambda = 1), nrow = 10, ncol = 12)
#'   rownames(counts) <- c(
#'     "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
#'     "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
#'   )
#'   colnames(counts) <- paste0("cell", 1:12)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   enrich_tbl <- tibble::tibble(
#'     ID = c("GO:0001", "GO:0002"),
#'     Description = c("immune response", "lymphocyte activation"),
#'     NES = c(2.1, 1.7),
#'     p.adjust = c(0.01, 0.03)
#'   )
#'   obj <- sn_store_enrichment(obj, enrich_tbl, store_name = "demo_gsea")
#'   names(obj@misc$enrichment_results)
#' }
#' @export
sn_store_enrichment <- function(object,
                                result,
                                store_name = "default",
                                analysis = c("ora", "gsea"),
                                database = "GOBP",
                                species = NULL,
                                source_de_name = NULL,
                                gene_col = "gene",
                                score_col = NULL,
                                return_object = TRUE) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  analysis <- match.arg(analysis)
  stored_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = .sn_as_enrichment_table(result),
    analysis = analysis,
    database = database,
    species = species,
    source_de_name = source_de_name,
    gene_col = gene_col,
    score_col = score_col
  )

  object <- .sn_store_misc_result(
    object = object,
    collection = "enrichment_results",
    store_name = store_name,
    result = stored_result
  )

  if (return_object) {
    return(.sn_log_seurat_command(object = object, name = "sn_store_enrichment"))
  }

  stored_result
}

#' Prepare cluster-annotation evidence from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param de_name Optional stored marker-result name in
#'   \code{object@misc$de_results}. When omitted, Shennong prefers
#'   \code{"default"}, then a single available result, and otherwise the most
#'   recent marker result.
#' @param cluster_col Metadata column containing cluster labels.
#' @param n_markers Number of top markers to retain per cluster.
#' @param marker_selection How to choose marker genes for annotation evidence:
#'   \code{"specific"} prefers genes that are relatively unique to one cluster,
#'   while \code{"top"} keeps the raw top-ranked genes.
#' @param enrichment_name Optional stored enrichment result used to add
#'   cluster-level functional evidence to the annotation prompt.
#' @param n_terms Number of enrichment terms to retain per cluster when
#'   \code{enrichment_name} is supplied.
#' @param enrichment_selection How to choose pathway/function terms for
#'   annotation evidence: \code{"specific"} prefers terms concentrated in fewer
#'   clusters, while \code{"top"} keeps the raw top-ranked terms.
#' @param include_qc Logical; whether to attach cluster-level QC summaries such
#'   as mitochondrial burden, failed-QC fractions, and doublet fractions when
#'   available in metadata.
#' @param reduction Optional dimensional reduction name used to summarize
#'   cluster neighborhood geometry, for example \code{"umap"}. Use
#'   \code{NULL} to disable geometry evidence.
#' @param n_neighbor_clusters Number of nearest clusters to report from the
#'   reduction centroid distances.
#'   Canonical lineage heuristic hints derived from known marker programs are
#'   included automatically when the required genes are present.
#'
#' @return A structured list ready for prompt construction.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(20 * 24, lambda = 1), nrow = 20, ncol = 24)
#'   rownames(counts) <- c(
#'     paste0("GENE", 1:14),
#'     "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "HLA-DRA"
#'   )
#'   colnames(counts) <- paste0("cell", 1:24)
#'   counts[c("CD3D", "CD3E", "TRAC"), 1:12] <-
#'     counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
#'   counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <-
#'     counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   evidence <- sn_prepare_annotation_evidence(
#'     obj,
#'     de_name = "celltype_markers",
#'     cluster_col = "cell_type"
#'   )
#'   names(evidence)
#' }
#' @export
sn_prepare_annotation_evidence <- function(object,
                                           de_name = NULL,
                                           cluster_col = "seurat_clusters",
                                           n_markers = 10,
                                           marker_selection = c("specific", "top"),
                                           enrichment_name = NULL,
                                           n_terms = 5,
                                           enrichment_selection = c("specific", "top"),
                                           include_qc = TRUE,
                                           reduction = "umap",
                                           n_neighbor_clusters = 3) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }
  marker_selection <- match.arg(marker_selection)
  enrichment_selection <- match.arg(enrichment_selection)

  de_name <- .sn_resolve_misc_result_name(
    object = object,
    collection = "de_results",
    store_name = de_name,
    preferred_analysis = "markers",
    arg_name = "de_name"
  )
  de_result <- .sn_get_misc_result(object = object, collection = "de_results", store_name = de_name)
  cluster_summary <- .sn_prepare_cluster_summary(object = object, cluster_col = cluster_col)
  marker_summary <- .sn_prepare_marker_summary(
    de_result = de_result,
    n_markers = n_markers,
    selection = marker_selection
  )
  marker_table <- .sn_prepare_marker_table(
    de_result = de_result,
    n_markers = n_markers,
    selection = marker_selection
  )
  prediction_summary <- .sn_prepare_prediction_summary(object = object, cluster_col = cluster_col)
  enrichment_summary <- .sn_prepare_cluster_enrichment_summary(
    object = object,
    enrichment_name = enrichment_name,
    n_terms = n_terms,
    selection = enrichment_selection
  )
  qc_summary <- if (isTRUE(include_qc)) {
    .sn_prepare_annotation_qc_summary(object = object, cluster_col = cluster_col)
  } else {
    tibble::tibble()
  }
  lineage_hints <- .sn_prepare_annotation_lineage_hints(
    object = object,
    cluster_col = cluster_col,
    top_marker_table = marker_table,
    species = tryCatch(sn_get_species(object), error = function(...) NULL)
  )
  canonical_marker_snapshot <- .sn_prepare_annotation_canonical_snapshot(
    object = object,
    cluster_col = cluster_col,
    species = tryCatch(sn_get_species(object), error = function(...) NULL)
  )
  geometry_summary <- .sn_prepare_cluster_geometry_summary(
    object = object,
    cluster_col = cluster_col,
    reduction = reduction,
    n_neighbors = n_neighbor_clusters
  )

  merged_summary <- dplyr::left_join(cluster_summary, marker_summary, by = "cluster")
  if (nrow(prediction_summary) > 0) {
    merged_summary <- dplyr::left_join(merged_summary, prediction_summary, by = "cluster")
  }
  if (nrow(enrichment_summary) > 0) {
    merged_summary <- dplyr::left_join(merged_summary, enrichment_summary, by = "cluster")
  }
  if (nrow(qc_summary) > 0) {
    qc_join <- dplyr::select(qc_summary, -dplyr::any_of("n_cells"))
    merged_summary <- dplyr::left_join(merged_summary, qc_join, by = "cluster")
  }
  if (nrow(lineage_hints) > 0) {
    merged_summary <- dplyr::left_join(merged_summary, lineage_hints, by = "cluster")
  }
  if (nrow(geometry_summary) > 0) {
    geometry_join <- dplyr::select(geometry_summary, -dplyr::any_of("geometry_reduction"))
    merged_summary <- dplyr::left_join(merged_summary, geometry_join, by = "cluster")
  }

  list(
    task = "annotation",
    cluster_col = cluster_col,
    source_de_name = de_name,
    source_enrichment_name = enrichment_name,
    analysis_method = de_result$method,
    species = tryCatch(sn_get_species(object), error = function(...) NULL),
    marker_selection = marker_selection,
    enrichment_selection = enrichment_selection,
    geometry_reduction = if (nrow(geometry_summary) > 0) reduction else NULL,
    cluster_summary = merged_summary,
    top_marker_table = marker_table,
    enrichment_summary = enrichment_summary,
    qc_summary = qc_summary,
    lineage_hints = lineage_hints,
    canonical_marker_snapshot = canonical_marker_snapshot,
    geometry_summary = geometry_summary,
    caveats = character()
  )
}

#' Prepare differential-expression evidence from a stored DE result
#'
#' @param object A \code{Seurat} object.
#' @param de_name Name of a stored DE result in \code{object@misc$de_results}.
#' @param n_genes Number of top genes to include per direction or group.
#'
#' @return A structured list ready for prompt construction.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(20 * 24, lambda = 1), nrow = 20, ncol = 24)
#'   rownames(counts) <- c(
#'     paste0("GENE", 1:14),
#'     "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "HLA-DRA"
#'   )
#'   colnames(counts) <- paste0("cell", 1:24)
#'   counts[c("CD3D", "CD3E", "TRAC"), 1:12] <-
#'     counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
#'   counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <-
#'     counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   evidence <- sn_prepare_de_evidence(obj, de_name = "celltype_markers", n_genes = 3)
#'   names(evidence)
#' }
#' @export
sn_prepare_de_evidence <- function(object, de_name, n_genes = 15) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  de_result <- .sn_get_misc_result(object = object, collection = "de_results", store_name = de_name)
  result_table <- tibble::as_tibble(de_result$table)
  rank_col <- de_result$rank_col
  if (is_null(rank_col) || !rank_col %in% colnames(result_table)) {
    stop("The stored DE result does not contain a ranking column.")
  }

  summary <- list(
    analysis = de_result$analysis,
    method = de_result$method,
    group_by = de_result$group_by,
    ident_1 = de_result$ident_1,
    ident_2 = de_result$ident_2,
    subset_by = de_result$subset_by
  )

  if (identical(de_result$analysis, "markers")) {
    top_table <- .sn_prepare_marker_summary(de_result = de_result, n_markers = n_genes)
    marker_table <- .sn_prepare_marker_table(de_result = de_result, n_markers = n_genes)
    return(list(
      task = "de",
      source_de_name = de_name,
      summary = summary,
      top_markers = top_table,
      top_marker_table = marker_table,
      caveats = character()
    ))
  }

  ordered <- result_table[order(-abs(result_table[[rank_col]])), , drop = FALSE]
  top_hits <- utils::head(ordered, n_genes)
  top_up <- ordered[ordered[[rank_col]] > 0, , drop = FALSE] |> utils::head(n_genes)
  top_down <- ordered[ordered[[rank_col]] < 0, , drop = FALSE] |> utils::head(n_genes)

  list(
    task = "de",
    source_de_name = de_name,
    summary = summary,
    top_hits = tibble::as_tibble(top_hits),
    top_up = tibble::as_tibble(top_up),
    top_down = tibble::as_tibble(top_down),
    caveats = character()
  )
}

#' Prepare enrichment evidence
#'
#' @param object Optional \code{Seurat} object containing stored enrichment
#'   results.
#' @param enrichment_name Name of a stored enrichment result.
#' @param result Optional enrichment result object supplied directly.
#' @param n_terms Number of top terms to keep.
#'
#' @return A structured list ready for prompt construction.
#'
#' @examples
#' enrich_tbl <- tibble::tibble(
#'   ID = c("GO:0001", "GO:0002"),
#'   Description = c("immune response", "lymphocyte activation"),
#'   NES = c(2.1, 1.7),
#'   p.adjust = c(0.01, 0.03)
#' )
#' evidence <- sn_prepare_enrichment_evidence(result = enrich_tbl, n_terms = 1)
#' evidence$top_terms
#' @export
sn_prepare_enrichment_evidence <- function(object = NULL,
                                           enrichment_name = NULL,
                                           result = NULL,
                                           n_terms = 10) {
  if (is_null(object) && is_null(result)) {
    stop("Supply either `object` + `enrichment_name` or `result`.")
  }

  stored <- if (!is_null(object)) {
    .sn_get_misc_result(object = object, collection = "enrichment_results", store_name = enrichment_name)
  } else {
    list(
      table = .sn_as_enrichment_table(result),
      analysis = NA_character_,
      database = NA_character_,
      species = NA_character_,
      source_de_name = NULL
    )
  }

  table <- tibble::as_tibble(stored$table)
  rank_candidates <- c("NES", "Count", "GeneRatio", "p.adjust", "pvalue")
  rank_col <- rank_candidates[rank_candidates %in% colnames(table)][1]
  ordered <- if (is_null(rank_col)) table else table[order(if (rank_col == "p.adjust" || rank_col == "pvalue") table[[rank_col]] else -abs(table[[rank_col]])), , drop = FALSE]

  list(
    task = "enrichment",
    source_enrichment_name = enrichment_name,
    analysis = stored$analysis,
    database = stored$database,
    species = stored$species,
    source_de_name = stored$source_de_name,
    top_terms = tibble::as_tibble(utils::head(ordered, n_terms)),
    full_term_table = tibble::as_tibble(table),
    caveats = character()
  )
}

#' Prepare manuscript-style results evidence
#'
#' @param object A \code{Seurat} object.
#' @param cluster_de_name Optional stored cluster-marker result.
#' @param contrast_de_name Optional stored contrast or pseudobulk result.
#' @param enrichment_name Optional stored enrichment result.
#' @param cluster_col Metadata column containing cluster labels.
#' @param n_markers Number of marker genes to retain per cluster.
#' @param n_terms Number of enrichment terms to retain.
#'
#' @return A structured list ready for prompt construction.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(10 * 24, lambda = 1), nrow = 10, ncol = 24)
#'   rownames(counts) <- c(
#'     "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
#'     "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
#'   )
#'   colnames(counts) <- paste0("cell", 1:24)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   obj <- sn_store_enrichment(
#'     obj,
#'     tibble::tibble(ID = "GO:0001", Description = "immune response", NES = 2, p.adjust = 0.01),
#'     store_name = "demo_gsea"
#'   )
#'   evidence <- sn_prepare_results_evidence(
#'     obj,
#'     cluster_de_name = "celltype_markers",
#'     enrichment_name = "demo_gsea",
#'     cluster_col = "cell_type"
#'   )
#'   names(evidence)
#' }
#' @export
sn_prepare_results_evidence <- function(object,
                                        cluster_de_name = NULL,
                                        contrast_de_name = NULL,
                                        enrichment_name = NULL,
                                        cluster_col = "seurat_clusters",
                                        n_markers = 5,
                                        n_terms = 10) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  evidence <- list(
    task = "results",
    dataset = list(
      n_cells = ncol(object),
      n_features = nrow(object),
      cluster_col = cluster_col,
      clusters = if (cluster_col %in% colnames(object[[]])) nlevels(factor(object[[cluster_col]][, 1])) else NULL
    ),
    cluster_summary = .sn_prepare_cluster_summary(object = object, cluster_col = cluster_col)
  )

  if (!is_null(cluster_de_name)) {
    evidence$cluster_markers <- sn_prepare_annotation_evidence(
      object = object,
      de_name = cluster_de_name,
      cluster_col = cluster_col,
      n_markers = n_markers
    )$cluster_summary
  }

  if (!is_null(contrast_de_name)) {
    evidence$de_summary <- sn_prepare_de_evidence(
      object = object,
      de_name = contrast_de_name,
      n_genes = n_markers
    )
  }

  if (!is_null(enrichment_name)) {
    evidence$enrichment_summary <- sn_prepare_enrichment_evidence(
      object = object,
      enrichment_name = enrichment_name,
      n_terms = n_terms
    )
  }

  evidence
}

#' Build an LLM prompt from structured Shennong evidence
#'
#' @param evidence A structured evidence list created by
#'   \code{sn_prepare_*_evidence()}.
#' @param task Interpretation task type.
#' @param style Optional style instruction, for example \code{"manuscript"}.
#' @param audience Intended audience such as \code{"scientist"}.
#' @param language Output language.
#' @param background Optional user-supplied study background or biological
#'   context to inject into the prompt.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable markdown brief.
#' @param include_json_schema Whether to request structured JSON output.
#'
#' @return A prompt bundle with \code{system}, \code{user}, and \code{messages}.
#'
#' @examples
#' evidence <- list(
#'   task = "annotation",
#'   cluster_summary = data.frame(cluster = "0", top_markers = "CD3D, TRAC")
#' )
#' prompt <- sn_build_prompt(evidence = evidence, task = "annotation")
#' names(prompt)
#' @export
sn_build_prompt <- function(evidence,
                            task = c("annotation", "de", "enrichment", "results", "figure_legend", "presentation_summary"),
                            style = NULL,
                            audience = c("scientist", "clinician", "general"),
                            language = c("en", "zh"),
                            background = NULL,
                            output_format = c("llm", "human"),
                            include_json_schema = FALSE) {
  task <- match.arg(task)
  audience <- match.arg(audience)
  language <- match.arg(language)
  output_format <- match.arg(output_format)

  instruction <- .sn_interpretation_task_instructions(task = task, evidence = evidence)
  evidence_max_rows <- if (identical(task, "annotation")) {
    200L
  } else {
    8L
  }

  if (identical(output_format, "human")) {
    return(.sn_human_readable_prompt(
      task = task,
      evidence = evidence,
      background = background,
      instruction = instruction
    ))
  }

  system_prompt <- paste(
    .sn_render_template(file.path("interpretation", "system_prompt.txt")),
    collapse = " "
  )

  style_line <- if (!is_null(style)) paste0("- Target style: ", style) else NULL
  json_line <- if (isTRUE(include_json_schema)) {
    if (identical(task, "annotation")) {
      .sn_annotation_json_schema_text()
    } else {
      "Return a structured JSON object followed by a brief narrative explanation."
    }
  } else {
    NULL
  }
  metadata_lines <- c(
    paste0("- Task: ", task),
    paste0("- Audience: ", audience),
    paste0("- Language: ", language),
    style_line
  )

  user_sections <- c(
    "# Interpretation Request",
    "## Task Metadata",
    paste(metadata_lines[nzchar(metadata_lines)], collapse = "\n"),
    "## Task Instructions",
    instruction
  )

  if (!is_null(background) && nzchar(background)) {
    user_sections <- c(
      user_sections,
      "## Background Context",
      as.character(background)
    )
  }
  if (!is_null(json_line) && nzchar(json_line)) {
    user_sections <- c(
      user_sections,
      "## Output Contract",
      json_line
    )
  }
  user_sections <- c(
    user_sections,
    "## Evidence",
    .sn_render_evidence_markdown(evidence, max_rows = evidence_max_rows)
  )
  user_prompt <- paste(user_sections, collapse = "\n\n")

  list(
    output_format = output_format,
    task = task,
    system = system_prompt,
    user = user_prompt,
    messages = .sn_build_messages(system_prompt = system_prompt, user_prompt = user_prompt),
    evidence = evidence
  )
}

#' Run an LLM provider on a prepared message list
#'
#' @param messages A message list, typically from \code{sn_build_prompt()}.
#' @param provider A user-supplied function that accepts \code{messages} and
#'   returns text or a list containing \code{text}.
#' @param model Optional model identifier passed through to the provider.
#' @param structured_type Optional structured-output schema or type object
#'   forwarded to providers that support typed responses.
#' @param tools Optional list of tool definitions forwarded to providers that
#'   support tool registration or tool calling.
#' @param ... Additional arguments passed to \code{provider}.
#' @importFrom utils modifyList
#'
#' @return A list containing at least \code{text}.
#'
#' @examples
#' provider <- function(messages, model = NULL, ...) {
#'   list(text = paste("received", length(messages), "messages"))
#' }
#' sn_run_llm(
#'   messages = list(list(role = "user", content = "Summarize this result.")),
#'   provider = provider
#' )
#' @export
sn_run_llm <- function(messages, provider, model = NULL, structured_type = NULL, tools = NULL, ...) {
  if (!is.function(provider)) {
    stop("`provider` must be a function.")
  }

  provider_args <- list(
    messages = messages,
    model = model,
    structured_type = structured_type,
    tools = tools,
    ...
  )
  provider_formals <- names(formals(provider) %||% alist(... = ))
  if (!"..." %in% provider_formals) {
    provider_args <- provider_args[intersect(names(provider_args), provider_formals)]
  }

  response <- do.call(provider, provider_args)
  if (is.character(response) && length(response) == 1) {
    return(list(text = response, model = model))
  }
  if (is.list(response) && any(c("text", "structured") %in% names(response))) {
    if (!"text" %in% names(response) && "structured" %in% names(response)) {
      response$text <- tryCatch(
        jsonlite::toJSON(response$structured, auto_unbox = TRUE, null = "null", dataframe = "rows"),
        error = function(...) NULL
      )
    }
    response$text <- .sn_text_scalar(response$text)
    return(response)
  }

  stop("`provider` must return either a single string or a list containing `text` and/or `structured`.")
}

.sn_text_scalar <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.character(x) && length(x) >= 1L) {
    return(as.character(x[[1]]))
  }

  values <- unlist(x, recursive = TRUE, use.names = FALSE)
  values <- values[!is.na(values)]
  if (length(values) == 0L) {
    return(NULL)
  }
  as.character(values[[1]])
}

.sn_format_openai_base_url <- function(base_url,
                                       require_v1 = TRUE) {
  if (!nzchar(base_url)) {
    stop("`base_url` must be a non-empty URL.", call. = FALSE)
  }

  base_url <- sub("/+$", "", base_url)
  if (grepl("/(responses|chat/completions)$", base_url)) {
    return(base_url)
  }
  if (require_v1 && !grepl("/v1$", base_url)) {
    return(paste0(base_url, "/v1"))
  }
  base_url
}

.sn_extract_openai_response_text <- function(parsed,
                                             wire_api = c("responses", "chat_completions")) {
  wire_api <- match.arg(wire_api)

  if (identical(wire_api, "chat_completions")) {
    choices <- parsed$choices %||% NULL
    choice <- if (is.data.frame(choices)) {
      as.list(choices[1, , drop = FALSE])
    } else if (is.list(choices) && length(choices) > 0) {
      choices[[1]]
    } else {
      NULL
    }
    message <- choice$message %||% NULL
    text <- if (is.list(message)) {
      message$content %||% NULL
    } else if (is.data.frame(message) && "content" %in% colnames(message)) {
      message$content[[1]]
    } else {
      choice$content %||% choice$text %||% NULL
    }
    return(text)
  }

  output_items <- parsed$output %||% list()
  output_text <- NULL
  if (is.list(output_items) && length(output_items) > 0L) {
    for (item in output_items) {
      contents <- item$content %||% list()
      if (is.list(contents) && length(contents) > 0L) {
        for (content_item in contents) {
          if (identical(content_item$type %||% NULL, "output_text") && !is.null(content_item$text)) {
            output_text <- content_item$text
            break
          }
        }
      }
      if (!is.null(output_text)) {
        break
      }
      if (!is.null(item$text)) {
        output_text <- item$text
        break
      }
    }
  }

  parsed$output_text %||%
    output_text %||%
    parsed$text %||%
    NULL
}

#' Create an \pkg{ellmer}-backed provider for Shennong interpretation helpers
#'
#' This adapter is useful when you want to manage transport, streaming, or
#' future structured-output features through \pkg{ellmer} while keeping
#' Shennong's interpretation API unchanged.
#'
#' @param api_key API key. Defaults to \code{OPENAI_API_KEY}.
#' @param base_url Base URL of the API. Defaults to \code{OPENAI_BASE_URL},
#'   then \code{"https://api.openai.com/v1"}.
#' @param model Default model identifier.
#' @param provider_type One of \code{"openai_compatible"} or \code{"openai"}.
#' @param echo Echo mode forwarded to \pkg{ellmer}.
#' @param reasoning_effort Optional reasoning effort forwarded to compatible
#'   GPT-5 chat-completions providers.
#' @param api_args Optional named list appended to each API request.
#' @param retries Number of retry attempts for retryable upstream proxy / HTML
#'   response failures.
#' @param retry_delay_sec Delay between retry attempts in seconds.
#'
#' @return A provider function suitable for \code{sn_run_llm()}.
#'
#' @examples
#' \dontrun{
#' provider <- sn_make_ellmer_provider(
#'   api_key = Sys.getenv("OPENAI_API_KEY"),
#'   base_url = "https://api.catplot.org/v1",
#'   model = "gpt-5.4"
#' )
#' }
#' @export
sn_make_ellmer_provider <- function(api_key = .sn_default_llm_api_key(),
                                    base_url = .sn_default_llm_base_url(),
                                    model = NULL,
                                    provider_type = c("openai_compatible", "openai"),
                                    echo = c("none", "output", "all"),
                                    reasoning_effort = NULL,
                                    api_args = list(),
                                    retries = 2L,
                                    retry_delay_sec = 1) {
  provider_type <- match.arg(provider_type)
  echo <- match.arg(echo)

  check_installed("ellmer")
  if (!nzchar(api_key)) {
    stop("`api_key` is required. Set `OPENAI_API_KEY` or pass it explicitly.", call. = FALSE)
  }

  default_model <- model
  formatted_base <- .sn_format_openai_base_url(base_url)
  provider_api_args <- api_args
  if (!is.null(reasoning_effort) && nzchar(reasoning_effort)) {
    provider_api_args$reasoning_effort <- reasoning_effort
  }

  function(messages, model = NULL, structured_type = NULL, tools = NULL, ...) {
    system_messages <- vapply(
      Filter(function(message) identical(message$role, "system"), messages),
      function(message) as.character(message$content %||% ""),
      character(1)
    )
    system_prompt <- if (length(system_messages) > 0L) {
      paste(system_messages, collapse = "\n\n")
    } else {
      NULL
    }
    non_system_messages <- Filter(function(message) !identical(message$role, "system"), messages)
    prompt_text <- paste(
      vapply(
        non_system_messages,
        function(message) paste0(toupper(message$role), ":\n", as.character(message$content %||% "")),
        character(1)
      ),
      collapse = "\n\n"
    )

    credentials <- function() api_key
    attempts <- max(1L, as.integer(retries %||% 1L))
    last_error <- NULL

    for (attempt in seq_len(attempts)) {
      result <- tryCatch({
        chat <- if (identical(provider_type, "openai")) {
          ellmer::chat_openai(
            system_prompt = system_prompt,
            base_url = formatted_base,
            credentials = credentials,
            model = model %||% default_model,
            api_args = modifyList(provider_api_args, list(...)),
            echo = echo
          )
        } else {
          ellmer::chat_openai_compatible(
            base_url = formatted_base,
            name = "Shennong provider",
            system_prompt = system_prompt,
            credentials = credentials,
            model = model %||% default_model,
            api_args = modifyList(provider_api_args, list(...)),
            echo = echo
          )
        }

        if (length(tools %||% list()) > 0L && is.null(structured_type)) {
          for (current_tool in tools) {
            chat$register_tool(current_tool)
          }
        }

        if (!is.null(structured_type)) {
          structured <- chat$chat_structured(
            prompt_text,
            type = structured_type,
            echo = echo
          )
          return(list(
            text = jsonlite::toJSON(structured, auto_unbox = TRUE, null = "null", dataframe = "rows"),
            structured = structured,
            model = model %||% default_model,
            raw = NULL
          ))
        }

        text <- chat$chat(prompt_text)
        list(
          text = as.character(text),
          model = model %||% default_model,
          raw = NULL
        )
      }, error = identity)

      if (!inherits(result, "error")) {
        return(result)
      }

      last_error <- result
      if (attempt < attempts && .sn_is_retryable_llm_error(conditionMessage(result))) {
        Sys.sleep(retry_delay_sec)
        next
      }
      break
    }

    stop(.sn_harden_llm_error(conditionMessage(last_error)), call. = FALSE)
  }
}

#' Test whether an ellmer-backed LLM provider is reachable and usable
#'
#' @param name Optional label used in the returned summary table.
#' @param model Optional override model used for the test request.
#' @param prompt Prompt used for the connectivity check.
#' @param provider Optional provider function. When supplied, Shennong tests
#'   this function directly; otherwise it builds an ellmer-backed provider from
#'   environment variables.
#'
#' @return A one-row tibble summarizing the result.
#'
#' @examples
#' test_provider <- function(messages, model = NULL, ...) list(text = "OK", model = model %||% "demo")
#' sn_test_llm_provider(provider = test_provider)
#' @export
sn_test_llm_provider <- function(name = NULL,
                                 model = NULL,
                                 prompt = "Reply with exactly OK.",
                                 provider = NULL) {
  provider_name <- name %||% "ellmer"
  provider <- provider %||% .sn_get_default_ellmer_provider()

  started_at <- Sys.time()
  result <- tryCatch(
    sn_run_llm(
      messages = list(list(role = "user", content = prompt)),
      provider = provider,
      model = model
    ),
    error = identity
  )
  elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))

  if (inherits(result, "error")) {
    return(tibble::tibble(
      provider = provider_name,
      ok = FALSE,
      model = model %||% NA_character_,
      elapsed_sec = elapsed,
      text = NA_character_,
      error = conditionMessage(result)
    ))
  }

  tibble::tibble(
    provider = provider_name,
    ok = TRUE,
    model = result$model %||% model %||% NA_character_,
    elapsed_sec = elapsed,
    text = result$text %||% NA_character_,
    error = NA_character_
  )
}

#' Interpret cluster markers for cell-type annotation
#'
#' @param object A \code{Seurat} object.
#' @param de_name Optional stored marker-result name. When omitted, Shennong
#'   prefers \code{"default"}, then a single available result, and otherwise
#'   the most recent marker result.
#' @param cluster_col Metadata column containing cluster labels.
#' @param n_markers Number of top markers per cluster.
#' @param marker_selection How to choose marker genes for annotation evidence:
#'   \code{"specific"} prefers genes that are relatively unique to one cluster,
#'   while \code{"top"} keeps the raw top-ranked genes.
#' @param enrichment_name Optional stored enrichment result used to add
#'   cluster-level functional evidence.
#' @param n_terms Number of enrichment terms per cluster when
#'   \code{enrichment_name} is supplied.
#' @param enrichment_selection How to choose pathway/function terms for
#'   annotation evidence: \code{"specific"} prefers terms concentrated in fewer
#'   clusters, while \code{"top"} keeps the raw top-ranked terms.
#' @param include_qc Logical; whether to include cluster-level QC summaries in
#'   the evidence bundle.
#' @param reduction Optional dimensional reduction name used to summarize
#'   cluster neighborhood geometry, for example \code{"umap"}. Use
#'   \code{NULL} to disable geometry evidence.
#' @param n_neighbor_clusters Number of nearest clusters to report from the
#'   reduction centroid distances.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param annotation_mode Annotation workflow mode. \code{"single_pass"} sends
#'   one compact prompt. \code{"agentic"} runs a broad-pass annotation followed
#'   by a focused refinement pass on ambiguous or ILC-relevant clusters, then
#'   reconciles the result against canonical lineage hints.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function. When left \code{NULL},
#'   Shennong will try to construct an \pkg{ellmer}-backed provider from
#'   \code{OPENAI_API_KEY} plus optional \code{OPENAI_BASE_URL} and
#'   \code{OPENAI_MODEL} environment variables.
#' @param model Optional model identifier.
#' @param reasoning_effort Optional reasoning effort forwarded to compatible
#'   GPT-5 chat-completions providers, for example \code{"minimal"},
#'   \code{"low"}, \code{"medium"}, \code{"high"}, or \code{"xhigh"} when the
#'   upstream endpoint supports it.
#' @param include_json_schema Logical; whether to request structured JSON output
#'   from the provider. Defaults to \code{TRUE} for annotation workflows.
#' @param apply_metadata Logical; if \code{TRUE} and a structured annotation
#'   response is returned, map the cluster labels back onto each cell in the
#'   Seurat metadata.
#' @param metadata_prefix Prefix used for metadata columns written back to the
#'   Seurat object when \code{apply_metadata = TRUE}.
#' @param metadata_fields Annotation fields to write back into Seurat metadata
#'   when \code{apply_metadata = TRUE}. Defaults to core fields only:
#'   \code{primary_label}, \code{broad_label}, \code{confidence},
#'   \code{status}, and \code{risk_flags}. Detailed evidence remains in the
#'   stored interpretation table under \code{object@misc}.
#' @param label_candidates Optional vector of candidate cell-type labels or
#'   broad lineages that should be treated as annotation priors for sorted or
#'   enriched datasets. These are injected into the prompt as constraints, but
#'   the model may still return a broader label when evidence is weak.
#' @param label_style Naming style used to normalize returned cell-type labels.
#'   One of \code{"title"}, \code{"snake"}, or \code{"asis"}.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param store_name Name used under \code{object@misc$interpretation_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#' @param show_progress Logical; if \code{TRUE}, emit step-wise progress logs
#'   and, when \pkg{cli} is available, a console progress bar while waiting for
#'   the LLM response.
#' @param ... Additional arguments forwarded to \code{provider}.
#'
#' @return A prompt bundle, response, or updated \code{Seurat} object.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(20 * 24, lambda = 1), nrow = 20, ncol = 24)
#'   rownames(counts) <- c(
#'     paste0("GENE", 1:14),
#'     "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "HLA-DRA"
#'   )
#'   colnames(counts) <- paste0("cell", 1:24)
#'   counts[c("CD3D", "CD3E", "TRAC"), 1:12] <-
#'     counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
#'   counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <-
#'     counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_interpret_annotation(
#'     obj,
#'     de_name = "celltype_markers",
#'     cluster_col = "cell_type",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_interpret_annotation <- function(object,
                                    de_name = NULL,
                                    cluster_col = "seurat_clusters",
                                    n_markers = 10,
                                    marker_selection = c("specific", "top"),
                                    enrichment_name = NULL,
                                    n_terms = 5,
                                    enrichment_selection = c("specific", "top"),
                                    include_qc = TRUE,
                                    reduction = "umap",
                                    n_neighbor_clusters = 3,
                                    background = NULL,
                                    annotation_mode = c("single_pass", "agentic"),
                                    output_format = c("llm", "human"),
                                    provider = NULL,
                                    model = NULL,
                                    reasoning_effort = NULL,
                                    include_json_schema = TRUE,
                                    apply_metadata = TRUE,
                                    metadata_prefix = "sn_annotation",
                                    metadata_fields = c("primary_label", "broad_label", "confidence", "status", "risk_flags"),
                                    label_candidates = NULL,
                                    label_style = c("title", "snake", "asis"),
                                    return_prompt = FALSE,
                                    store_name = "default",
                                    return_object = TRUE,
                                    show_progress = interactive(),
                                    ...) {
  output_format <- match.arg(output_format)
  marker_selection <- match.arg(marker_selection)
  enrichment_selection <- match.arg(enrichment_selection)
  label_style <- match.arg(label_style)
  annotation_mode <- match.arg(annotation_mode)
  progress_state <- .sn_interpret_progress_start(
    task = "interpret_annotation",
    enabled = show_progress,
    total_steps = if (identical(annotation_mode, "agentic")) 8L else 5L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing annotation evidence")
  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = de_name,
    cluster_col = cluster_col,
    n_markers = n_markers,
    marker_selection = marker_selection,
    enrichment_name = enrichment_name,
    n_terms = n_terms,
    enrichment_selection = enrichment_selection,
    include_qc = include_qc,
    reduction = reduction,
    n_neighbor_clusters = n_neighbor_clusters
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Building annotation prompt")

  if (identical(annotation_mode, "agentic")) {
    broad_evidence <- .sn_prepare_annotation_compact_evidence(
      evidence = evidence,
      include_marker_table = FALSE,
      include_enrichment_table = FALSE,
      include_qc = include_qc,
      include_canonical_snapshot = FALSE,
      include_geometry = !is.null(reduction)
    )
    broad_prompt <- sn_build_prompt(
      evidence = broad_evidence,
      task = "annotation",
      audience = "scientist",
      style = "broad cell type annotation",
      background = .sn_build_annotation_stage_background(
        background = background,
        label_candidates = label_candidates,
        stage = "broad_pass"
      ),
      output_format = output_format,
      include_json_schema = include_json_schema
    )

    return(.sn_finish_annotation_agentic(
      object = object,
      evidence = evidence,
      broad_prompt = broad_prompt,
      provider = provider,
      model = model,
      store_name = store_name,
      cluster_col = cluster_col,
      metadata_prefix = metadata_prefix,
      metadata_fields = metadata_fields,
      label_style = label_style,
      apply_metadata = apply_metadata,
      label_candidates = label_candidates,
      background = background,
      return_prompt = return_prompt,
      return_object = return_object,
      progress_state = progress_state,
      reasoning_effort = reasoning_effort,
      ...
    ))
  }

  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "annotation",
    audience = "scientist",
    style = "cell type annotation",
    background = .sn_build_annotation_stage_background(
      background = background,
      label_candidates = label_candidates,
      stage = "single_pass"
    ),
    output_format = output_format,
    include_json_schema = include_json_schema
  )

  .sn_finish_interpretation(
    object = object,
    task = "interpret_annotation",
    evidence = evidence,
    prompt = prompt,
    provider = provider,
    model = model,
    store_name = store_name,
    cluster_col = cluster_col,
    metadata_prefix = metadata_prefix,
    metadata_fields = metadata_fields,
    label_style = label_style,
    apply_metadata = apply_metadata,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    reasoning_effort = reasoning_effort,
    ...
  )
}

#' Interpret a stored differential-expression result
#'
#' @param object A \code{Seurat} object.
#' @param de_name Name of a stored DE result.
#' @param n_genes Number of top genes to retain.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param store_name Name used under \code{object@misc$interpretation_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#' @param show_progress Logical; if \code{TRUE}, emit step-wise progress logs
#'   and, when \pkg{cli} is available, a console progress bar while waiting for
#'   the LLM response.
#' @param ... Additional arguments forwarded to \code{provider}.
#'
#' @return A prompt bundle, response, or updated \code{Seurat} object.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(20 * 24, lambda = 1), nrow = 20, ncol = 24)
#'   rownames(counts) <- c(
#'     paste0("GENE", 1:14),
#'     "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "HLA-DRA"
#'   )
#'   colnames(counts) <- paste0("cell", 1:24)
#'   counts[c("CD3D", "CD3E", "TRAC"), 1:12] <-
#'     counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
#'   counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <-
#'     counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_interpret_de(obj, de_name = "celltype_markers", return_prompt = TRUE)
#'   prompt$task
#' }
#' @export
sn_interpret_de <- function(object,
                            de_name,
                            n_genes = 15,
                            background = NULL,
                            output_format = c("llm", "human"),
                            provider = NULL,
                            model = NULL,
                            return_prompt = FALSE,
                            store_name = "default",
                            return_object = TRUE,
                            show_progress = interactive(),
                            ...) {
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "interpret_de",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing DE evidence")
  evidence <- sn_prepare_de_evidence(object = object, de_name = de_name, n_genes = n_genes)
  progress_state <- .sn_interpret_progress_step(progress_state, "Building interpretation prompt")
  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "de",
    audience = "scientist",
    style = "differential expression interpretation",
    background = background,
    output_format = output_format
  )

  .sn_finish_interpretation(
    object = object,
    task = "interpret_de",
    evidence = evidence,
    prompt = prompt,
    provider = provider,
    model = model,
    store_name = store_name,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}

#' Interpret a stored enrichment result
#'
#' @param object A \code{Seurat} object.
#' @param enrichment_name Name of a stored enrichment result.
#' @param n_terms Number of top enrichment terms to retain.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param store_name Name used under \code{object@misc$interpretation_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#' @param show_progress Logical; if \code{TRUE}, emit step-wise progress logs
#'   and, when \pkg{cli} is available, a console progress bar while waiting for
#'   the LLM response.
#' @param ... Additional arguments forwarded to \code{provider}.
#'
#' @return A prompt bundle, response, or updated \code{Seurat} object.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(10 * 12, lambda = 1), nrow = 10, ncol = 12)
#'   rownames(counts) <- c(
#'     "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
#'     "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
#'   )
#'   colnames(counts) <- paste0("cell", 1:12)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj <- sn_store_enrichment(
#'     obj,
#'     tibble::tibble(ID = "GO:0001", Description = "immune response", NES = 2, p.adjust = 0.01),
#'     store_name = "demo_gsea"
#'   )
#'   prompt <- sn_interpret_enrichment(
#'     obj,
#'     enrichment_name = "demo_gsea",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_interpret_enrichment <- function(object,
                                    enrichment_name,
                                    n_terms = 10,
                                    background = NULL,
                                    output_format = c("llm", "human"),
                                    provider = NULL,
                                    model = NULL,
                                    return_prompt = FALSE,
                                    store_name = "default",
                                    return_object = TRUE,
                                    show_progress = interactive(),
                                    ...) {
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "interpret_enrichment",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing enrichment evidence")
  evidence <- sn_prepare_enrichment_evidence(
    object = object,
    enrichment_name = enrichment_name,
    n_terms = n_terms
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Building interpretation prompt")
  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "enrichment",
    audience = "scientist",
    style = "pathway interpretation",
    background = background,
    output_format = output_format
  )

  .sn_finish_interpretation(
    object = object,
    task = "interpret_enrichment",
    evidence = evidence,
    prompt = prompt,
    provider = provider,
    model = model,
    store_name = store_name,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}

#' Write a manuscript-style results summary from stored analysis outputs
#'
#' @param object A \code{Seurat} object.
#' @param cluster_de_name Optional stored cluster-marker result.
#' @param contrast_de_name Optional stored contrast or pseudobulk result.
#' @param enrichment_name Optional stored enrichment result.
#' @param cluster_col Metadata column containing cluster labels.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param store_name Name used under \code{object@misc$interpretation_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#' @param show_progress Logical; if \code{TRUE}, emit step-wise progress logs
#'   and, when \pkg{cli} is available, a console progress bar while waiting for
#'   the LLM response.
#' @param ... Additional arguments forwarded to \code{provider}.
#'
#' @return A prompt bundle, response, or updated \code{Seurat} object.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(10 * 24, lambda = 1), nrow = 10, ncol = 24)
#'   rownames(counts) <- c(
#'     "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
#'     "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
#'   )
#'   colnames(counts) <- paste0("cell", 1:24)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   obj <- sn_store_enrichment(
#'     obj,
#'     tibble::tibble(ID = "GO:0001", Description = "immune response", NES = 2, p.adjust = 0.01),
#'     store_name = "demo_gsea"
#'   )
#'   prompt <- sn_write_results(
#'     obj,
#'     cluster_de_name = "celltype_markers",
#'     enrichment_name = "demo_gsea",
#'     cluster_col = "cell_type",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_write_results <- function(object,
                             cluster_de_name = NULL,
                             contrast_de_name = NULL,
                             enrichment_name = NULL,
                             cluster_col = "seurat_clusters",
                             background = NULL,
                             output_format = c("llm", "human"),
                             provider = NULL,
                             model = NULL,
                             return_prompt = FALSE,
                             store_name = "default",
                             return_object = TRUE,
                             show_progress = interactive(),
                             ...) {
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "write_results",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing results evidence")
  evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_name = cluster_de_name,
    contrast_de_name = contrast_de_name,
    enrichment_name = enrichment_name,
    cluster_col = cluster_col
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Building writing prompt")
  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "results",
    style = "manuscript-style Results section",
    audience = "scientist",
    background = background,
    output_format = output_format
  )

  .sn_finish_interpretation(
    object = object,
    task = "write_results",
    evidence = evidence,
    prompt = prompt,
    provider = provider,
    model = model,
    store_name = store_name,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}

#' Write a figure legend from stored analysis outputs
#'
#' @param object A \code{Seurat} object.
#' @param cluster_de_name Optional stored cluster-marker result.
#' @param enrichment_name Optional stored enrichment result.
#' @param cluster_col Metadata column containing cluster labels.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param store_name Name used under \code{object@misc$interpretation_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#' @param show_progress Logical; if \code{TRUE}, emit step-wise progress logs
#'   and, when \pkg{cli} is available, a console progress bar while waiting for
#'   the LLM response.
#' @param ... Additional arguments forwarded to \code{provider}.
#'
#' @return A prompt bundle, response, or updated \code{Seurat} object.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(10 * 24, lambda = 1), nrow = 10, ncol = 24)
#'   rownames(counts) <- c(
#'     "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
#'     "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
#'   )
#'   colnames(counts) <- paste0("cell", 1:24)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_write_figure_legend(
#'     obj,
#'     cluster_de_name = "celltype_markers",
#'     cluster_col = "cell_type",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_write_figure_legend <- function(object,
                                   cluster_de_name = NULL,
                                   enrichment_name = NULL,
                                   cluster_col = "seurat_clusters",
                                   background = NULL,
                                   output_format = c("llm", "human"),
                                   provider = NULL,
                                   model = NULL,
                                   return_prompt = FALSE,
                                   store_name = "default",
                                   return_object = TRUE,
                                   show_progress = interactive(),
                                   ...) {
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "write_figure_legend",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing legend evidence")
  evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_name = cluster_de_name,
    enrichment_name = enrichment_name,
    cluster_col = cluster_col
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Building legend prompt")
  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "figure_legend",
    style = "figure legend",
    audience = "scientist",
    background = background,
    output_format = output_format
  )

  .sn_finish_interpretation(
    object = object,
    task = "write_figure_legend",
    evidence = evidence,
    prompt = prompt,
    provider = provider,
    model = model,
    store_name = store_name,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}

#' Write a presentation-style summary from stored analysis outputs
#'
#' @param object A \code{Seurat} object.
#' @param cluster_de_name Optional stored cluster-marker result.
#' @param contrast_de_name Optional stored contrast or pseudobulk result.
#' @param enrichment_name Optional stored enrichment result.
#' @param cluster_col Metadata column containing cluster labels.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param store_name Name used under \code{object@misc$interpretation_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#' @param show_progress Logical; if \code{TRUE}, emit step-wise progress logs
#'   and, when \pkg{cli} is available, a console progress bar while waiting for
#'   the LLM response.
#' @param ... Additional arguments forwarded to \code{provider}.
#'
#' @return A prompt bundle, response, or updated \code{Seurat} object.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(10 * 24, lambda = 1), nrow = 10, ncol = 24)
#'   rownames(counts) <- c(
#'     "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
#'     "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
#'   )
#'   colnames(counts) <- paste0("cell", 1:24)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_write_presentation_summary(
#'     obj,
#'     cluster_de_name = "celltype_markers",
#'     cluster_col = "cell_type",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_write_presentation_summary <- function(object,
                                          cluster_de_name = NULL,
                                          contrast_de_name = NULL,
                                          enrichment_name = NULL,
                                          cluster_col = "seurat_clusters",
                                          background = NULL,
                                          output_format = c("llm", "human"),
                                          provider = NULL,
                                          model = NULL,
                                          return_prompt = FALSE,
                                          store_name = "default",
                                          return_object = TRUE,
                                          show_progress = interactive(),
                                          ...) {
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "write_presentation_summary",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing presentation evidence")
  evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_name = cluster_de_name,
    contrast_de_name = contrast_de_name,
    enrichment_name = enrichment_name,
    cluster_col = cluster_col
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Building presentation prompt")
  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "presentation_summary",
    style = "presentation slide summary",
    audience = "scientist",
    background = background,
    output_format = output_format
  )

  .sn_finish_interpretation(
    object = object,
    task = "write_presentation_summary",
    evidence = evidence,
    prompt = prompt,
    provider = provider,
    model = model,
    store_name = store_name,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}
