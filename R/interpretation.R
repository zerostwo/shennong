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
    .sn_compact_collection_summary(object, "deconvolution_results", "deconvolution"),
    .sn_compact_collection_summary(object, "milo_results", "milo"),
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
    return(paste(utils::capture.output(print(utils::head(x, max_rows))), collapse = "\n"))
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

.sn_parse_annotation_response <- function(response) {
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
                                          metadata_prefix = "sn_annotation") {
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

  metadata_to_add <- data.frame(row.names = colnames(object))
  for (source_col in names(field_map)) {
    if (source_col %in% colnames(merged)) {
      metadata_to_add[[field_map[[source_col]]]] <- merged[colnames(object), source_col, drop = TRUE]
    }
  }
  metadata_to_add[[paste0(metadata_prefix, "_cluster")]] <- merged[colnames(object), "cluster", drop = TRUE]
  SeuratObject::AddMetaData(object, metadata = metadata_to_add)
}

.sn_annotation_json_schema_text <- function() {
  paste(
    "Return valid JSON with keys `cluster_annotations` and `narrative_summary`.",
    "`cluster_annotations` must be an array of objects with:",
    "`cluster`, `primary_label`, `broad_label`, `confidence`, `status`,",
    "`alternatives`, `supporting_markers`, `supporting_functions`, `risk_flags`,",
    "`note`, and `recommended_checks`.",
    "`confidence` should use one of: high, medium, low.",
    "`status` should use one of: confident, ambiguous, possible_contamination, possible_low_quality, possible_transition.",
    "`risk_flags` should be an array that may include: contamination, low_quality, transitional_state, doublet_like, proliferating_state."
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

.sn_finish_interpretation <- function(object,
                                      task,
                                      evidence,
                                      prompt,
                                      provider = NULL,
                                      model = NULL,
                                      store_name = "default",
                                      cluster_col = NULL,
                                      metadata_prefix = "sn_annotation",
                                      label_style = c("title", "snake", "asis"),
                                      apply_metadata = FALSE,
                                      return_prompt = FALSE,
                                      return_object = TRUE,
                                      ...) {
  label_style <- match.arg(label_style)
  if (return_prompt || identical(prompt$output_format, "human")) {
    return(prompt)
  }

  if (is_null(provider)) {
    provider <- tryCatch(
      sn_get_llm_provider(),
      error = function(...) NULL
    )
  }
  if (is_null(provider)) {
    stop("`provider` must be supplied unless `return_prompt = TRUE`, or you must configure a default provider under `~/.shennong`.", call. = FALSE)
  }

  response <- sn_run_llm(
    messages = prompt$messages,
    provider = provider,
    model = model,
    ...
  )

  parsed_annotation <- NULL
  if (identical(task, "interpret_annotation")) {
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
          metadata_prefix = metadata_prefix
        )
      }
    }
  }

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
    return(.sn_log_seurat_command(object = object, name = paste0("sn_", task)))
  }

  response
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

.sn_prepare_marker_summary <- function(de_result, n_markers = 10) {
  marker_table <- de_result$table
  group_col <- de_result$group_col
  rank_col <- de_result$rank_col
  p_col <- de_result$p_col

  if (nrow(marker_table) == 0) {
    return(tibble::tibble(cluster = character(), top_markers = character()))
  }
  if (is_null(group_col) || !group_col %in% colnames(marker_table)) {
    return(tibble::tibble(cluster = character(), top_markers = character()))
  }
  if (is_null(rank_col) || !rank_col %in% colnames(marker_table)) {
    return(tibble::tibble(cluster = character(), top_markers = character()))
  }

  working <- tibble::as_tibble(marker_table)
  if (!is_null(p_col) && p_col %in% colnames(working)) {
    working <- working[working[[p_col]] <= de_result$p_val_cutoff, , drop = FALSE]
  }
  if (nrow(working) == 0) {
    working <- tibble::as_tibble(marker_table)
  }

  ranking <- working[[rank_col]]
  if (!identical(de_result$analysis, "markers")) {
    ranking <- abs(ranking)
  }
  working$..ranking_value <- ranking

  top_markers <- working |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
    dplyr::slice_max(order_by = .data$..ranking_value, n = n_markers, with_ties = FALSE) |>
    dplyr::summarise(
      top_markers = .sn_compact_value(.data$gene, max_items = n_markers),
      .groups = "drop"
    ) |>
    dplyr::rename(cluster = dplyr::all_of(group_col))

  tibble::as_tibble(top_markers)
}

.sn_prepare_marker_table <- function(de_result, n_markers = 10) {
  marker_table <- tibble::as_tibble(de_result$table)
  group_col <- de_result$group_col
  rank_col <- de_result$rank_col

  if (nrow(marker_table) == 0 || is_null(group_col) || !group_col %in% colnames(marker_table)) {
    return(tibble::tibble())
  }

  ranking <- if (!is_null(rank_col) && rank_col %in% colnames(marker_table)) {
    if (identical(de_result$analysis, "markers")) marker_table[[rank_col]] else abs(marker_table[[rank_col]])
  } else {
    rep(0, nrow(marker_table))
  }
  marker_table$..ranking_value <- ranking

  marker_table |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
    dplyr::slice_max(order_by = .data$..ranking_value, n = n_markers, with_ties = FALSE) |>
    dplyr::ungroup()
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
                                                  n_terms = 5) {
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

  if (is.null(rank_col)) {
    ordered <- table
  } else if (rank_col %in% c("p.adjust", "pvalue")) {
    ordered <- table[order(table[[rank_col]], decreasing = FALSE), , drop = FALSE]
  } else {
    ordered <- table[order(abs(table[[rank_col]]), decreasing = TRUE), , drop = FALSE]
  }

  ordered |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) |>
    dplyr::slice_head(n = n_terms) |>
    dplyr::summarise(
      top_functions = .sn_compact_value(.data[[term_col]], max_items = n_terms),
      .groups = "drop"
    ) |>
    dplyr::rename(cluster = dplyr::all_of(group_col))
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
#' @param de_name Name of a stored marker result in \code{object@misc$de_results}.
#' @param cluster_col Metadata column containing cluster labels.
#' @param n_markers Number of top markers to retain per cluster.
#' @param enrichment_name Optional stored enrichment result used to add
#'   cluster-level functional evidence to the annotation prompt.
#' @param n_terms Number of enrichment terms to retain per cluster when
#'   \code{enrichment_name} is supplied.
#' @param include_qc Logical; whether to attach cluster-level QC summaries such
#'   as mitochondrial burden, failed-QC fractions, and doublet fractions when
#'   available in metadata.
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
                                           de_name,
                                           cluster_col = "seurat_clusters",
                                           n_markers = 10,
                                           enrichment_name = NULL,
                                           n_terms = 5,
                                           include_qc = TRUE) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  de_result <- .sn_get_misc_result(object = object, collection = "de_results", store_name = de_name)
  cluster_summary <- .sn_prepare_cluster_summary(object = object, cluster_col = cluster_col)
  marker_summary <- .sn_prepare_marker_summary(de_result = de_result, n_markers = n_markers)
  marker_table <- .sn_prepare_marker_table(de_result = de_result, n_markers = n_markers)
  prediction_summary <- .sn_prepare_prediction_summary(object = object, cluster_col = cluster_col)
  enrichment_summary <- .sn_prepare_cluster_enrichment_summary(
    object = object,
    enrichment_name = enrichment_name,
    n_terms = n_terms
  )
  qc_summary <- if (isTRUE(include_qc)) {
    .sn_prepare_annotation_qc_summary(object = object, cluster_col = cluster_col)
  } else {
    tibble::tibble()
  }

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

  list(
    task = "annotation",
    cluster_col = cluster_col,
    source_de_name = de_name,
    source_enrichment_name = enrichment_name,
    analysis_method = de_result$method,
    species = tryCatch(sn_get_species(object), error = function(...) NULL),
    cluster_summary = merged_summary,
    top_marker_table = marker_table,
    enrichment_summary = enrichment_summary,
    qc_summary = qc_summary,
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

  style_line <- if (!is_null(style)) paste0("Target style: ", style, ".") else ""
  json_line <- if (isTRUE(include_json_schema)) {
    if (identical(task, "annotation")) {
      .sn_annotation_json_schema_text()
    } else {
      "Return a structured JSON object followed by a brief narrative explanation."
    }
  } else {
    ""
  }
  background_line <- if (!is_null(background) && nzchar(background)) paste0("Background context:\n", background) else ""

  user_prompt <- paste(
    paste0("Task: ", task, "."),
    paste0("Audience: ", audience, "."),
    paste0("Language: ", language, "."),
    style_line,
    paste0("Task instructions: ", instruction),
    background_line,
    json_line,
    "Evidence:",
    .sn_render_prompt_value(evidence),
    sep = "\n\n"
  )

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
sn_run_llm <- function(messages, provider, model = NULL, ...) {
  if (!is.function(provider)) {
    stop("`provider` must be a function.")
  }

  response <- provider(messages = messages, model = model, ...)
  if (is.character(response) && length(response) == 1) {
    return(list(text = response, model = model))
  }
  if (is.list(response) && "text" %in% names(response)) {
    response$text <- .sn_text_scalar(response$text)
    return(response)
  }

  stop("`provider` must return either a single string or a list containing `text`.")
}

# Keep local provider configuration under ~/.shennong so package workflows can
# reuse a default endpoint without forcing users to edit .Renviron by hand.
.sn_llm_config_default <- function() {
  list(
    schema_version = "1.0.0",
    updated_at = NULL,
    default_provider = NULL,
    providers = list()
  )
}

.sn_llm_root_dir <- function(config_dir = "~/.shennong", create = TRUE) {
  path <- path.expand(config_dir)
  if (create && !dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

.sn_llm_config_path <- function(config_dir = "~/.shennong") {
  file.path(.sn_llm_root_dir(config_dir = config_dir), "llm-providers.json")
}

.sn_llm_history_dir <- function(config_dir = "~/.shennong", create = TRUE) {
  path <- file.path(.sn_llm_root_dir(config_dir = config_dir, create = create), "llm-history")
  if (create && !dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

.sn_coerce_named_list <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (is.list(x) && !is.data.frame(x)) {
    return(x)
  }
  as.list(x)
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

.sn_read_llm_config <- function(config_dir = "~/.shennong") {
  path <- .sn_llm_config_path(config_dir = config_dir)
  if (!file.exists(path)) {
    return(.sn_llm_config_default())
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(...) .sn_llm_config_default()
  )
  parsed <- modifyList(.sn_llm_config_default(), parsed)
  parsed$providers <- .sn_coerce_named_list(parsed$providers)
  parsed
}

.sn_write_llm_config <- function(config,
                                 config_dir = "~/.shennong") {
  path <- .sn_llm_config_path(config_dir = config_dir)
  config <- modifyList(.sn_llm_config_default(), config)
  config$providers <- .sn_coerce_named_list(config$providers)
  config$updated_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)

  writeLines(
    jsonlite::toJSON(config, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    con = path,
    useBytes = TRUE
  )
  suppressWarnings(Sys.chmod(path, mode = "600"))
  invisible(path)
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

.sn_openai_endpoint <- function(base_url,
                                wire_api = c("responses", "chat_completions")) {
  wire_api <- match.arg(wire_api)
  formatted_base <- .sn_format_openai_base_url(base_url)
  suffix <- if (identical(wire_api, "responses")) "responses" else "chat/completions"
  if (grepl(paste0("/", suffix, "$"), formatted_base, fixed = FALSE)) {
    return(formatted_base)
  }
  paste0(formatted_base, "/", suffix)
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

.sn_record_llm_history <- function(provider_name,
                                   base_url,
                                   wire_api,
                                   backend,
                                   messages,
                                   response,
                                   requested_model = NULL,
                                   config_dir = "~/.shennong") {
  history_dir <- .sn_llm_history_dir(config_dir = config_dir)
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
  history_path <- file.path(
    history_dir,
    paste0(stamp, "-", sprintf("%06d", sample.int(999999, 1)), ".json")
  )
  payload <- list(
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    provider_name = provider_name,
    base_url = base_url,
    wire_api = wire_api,
    backend = backend,
    requested_model = requested_model,
    resolved_model = response$model %||% requested_model,
    messages = messages,
    response = response
  )
  writeLines(
    jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    con = history_path,
    useBytes = TRUE
  )
  suppressWarnings(Sys.chmod(history_path, mode = "600"))
  invisible(history_path)
}

.sn_wrap_provider_with_history <- function(provider,
                                           provider_name,
                                           base_url,
                                           wire_api,
                                           backend,
                                           config_dir = "~/.shennong",
                                           history = TRUE) {
  force(provider)
  force(provider_name)
  force(base_url)
  force(wire_api)
  force(backend)
  force(config_dir)
  force(history)

  function(messages, model = NULL, ...) {
    response <- provider(messages = messages, model = model, ...)
    if (isTRUE(history)) {
      .sn_record_llm_history(
        provider_name = provider_name,
        base_url = base_url,
        wire_api = wire_api,
        backend = backend,
        messages = messages,
        response = response,
        requested_model = model,
        config_dir = config_dir
      )
    }
    response
  }
}

.sn_resolve_llm_provider_record <- function(name = NULL,
                                            config_dir = "~/.shennong") {
  config <- .sn_read_llm_config(config_dir = config_dir)
  providers <- config$providers %||% list()

  if (length(providers) == 0L) {
    stop("No LLM providers are configured under `~/.shennong`.", call. = FALSE)
  }

  selected_name <- name %||% config$default_provider
  if (is.null(selected_name) || !nzchar(selected_name)) {
    stop("No default LLM provider is configured.", call. = FALSE)
  }
  if (!selected_name %in% names(providers)) {
    stop(glue("Provider '{selected_name}' was not found in the local Shennong config."), call. = FALSE)
  }

  record <- providers[[selected_name]]
  record$name <- selected_name
  record
}

#' Configure a reusable local LLM provider for Shennong
#'
#' Stores provider metadata under \code{~/.shennong/llm-providers.json} so
#' Shennong interpretation helpers can reuse a default endpoint without forcing
#' every call site to pass an explicit provider function.
#'
#' @param name Provider name used in the local registry.
#' @param base_url Base URL of the upstream API. A bare host like
#'   \code{"https://api.catplot.org"} is accepted.
#' @param api_key API key to store in the local provider config. If omitted,
#'   Shennong will fall back to \code{Sys.getenv(api_key_env)} at runtime.
#' @param api_key_env Environment variable name used as a fallback secret source.
#' @param provider_type Provider family. Currently \code{"openai_compatible"}
#'   and \code{"openai"} are supported.
#' @param wire_api Upstream API style: \code{"responses"} or
#'   \code{"chat_completions"}.
#' @param default_model Default model used when no explicit \code{model} is
#'   supplied.
#' @param review_model Optional second model name kept for local bookkeeping.
#' @param reasoning_effort Optional reasoning level stored alongside the
#'   provider config.
#' @param backend_preference Preferred Shennong backend. Use
#'   \code{"native"} for the built-in HTTP client or \code{"ellmer"} for the
#'   optional \pkg{ellmer}-backed adapter.
#' @param disable_response_storage Logical; if \code{TRUE}, disable provider-side
#'   response storage when the upstream API supports it.
#' @param enable_history Logical; whether configured providers should write call
#'   history into \code{~/.shennong/llm-history/}.
#' @param set_default Logical; whether to mark this provider as the default one.
#' @param supports_websockets Logical flag stored for user inspection.
#' @param requires_openai_auth Logical flag stored for user inspection.
#' @param config_dir Local Shennong config root.
#'
#' @return Invisibly returns the updated provider table.
#'
#' @examples
#' cfg_dir <- tempfile("shennong-config-")
#' sn_configure_llm_provider(
#'   name = "demo",
#'   base_url = "https://example.invalid",
#'   api_key = "demo-key",
#'   default_model = "gpt-4.1",
#'   config_dir = cfg_dir
#' )
#' sn_list_llm_providers(config_dir = cfg_dir)
#' @export
sn_configure_llm_provider <- function(name = "default",
                                      base_url,
                                      api_key = NULL,
                                      api_key_env = "OPENAI_API_KEY",
                                      provider_type = c("openai_compatible", "openai"),
                                      wire_api = c("responses", "chat_completions"),
                                      default_model = NULL,
                                      review_model = NULL,
                                      reasoning_effort = NULL,
                                      backend_preference = c("native", "ellmer", "auto"),
                                      disable_response_storage = TRUE,
                                      enable_history = TRUE,
                                      set_default = TRUE,
                                      supports_websockets = FALSE,
                                      requires_openai_auth = TRUE,
                                      config_dir = "~/.shennong") {
  provider_type <- match.arg(provider_type)
  wire_api <- match.arg(wire_api)
  backend_preference <- match.arg(backend_preference)

  if (!nzchar(name)) {
    stop("`name` must be a non-empty provider name.", call. = FALSE)
  }

  config <- .sn_read_llm_config(config_dir = config_dir)
  providers <- config$providers %||% list()
  providers[[name]] <- list(
    provider_type = provider_type,
    base_url = .sn_format_openai_base_url(base_url, require_v1 = TRUE),
    api_key = if (!is.null(api_key) && nzchar(api_key)) api_key else NULL,
    api_key_env = api_key_env,
    wire_api = wire_api,
    default_model = default_model,
    review_model = review_model,
    reasoning_effort = reasoning_effort,
    backend_preference = backend_preference,
    disable_response_storage = isTRUE(disable_response_storage),
    enable_history = isTRUE(enable_history),
    supports_websockets = isTRUE(supports_websockets),
    requires_openai_auth = isTRUE(requires_openai_auth)
  )

  config$providers <- providers
  if (isTRUE(set_default) || is.null(config$default_provider)) {
    config$default_provider <- name
  }
  .sn_write_llm_config(config = config, config_dir = config_dir)
  invisible(sn_list_llm_providers(config_dir = config_dir))
}

#' List locally configured LLM providers
#'
#' @param config_dir Local Shennong config root.
#'
#' @return A tibble with one row per configured provider.
#'
#' @examples
#' cfg_dir <- tempfile("shennong-config-")
#' sn_configure_llm_provider(
#'   name = "demo",
#'   base_url = "https://example.invalid",
#'   api_key = "demo-key",
#'   default_model = "gpt-4.1",
#'   config_dir = cfg_dir
#' )
#' sn_list_llm_providers(config_dir = cfg_dir)
#' @export
sn_list_llm_providers <- function(config_dir = "~/.shennong") {
  config <- .sn_read_llm_config(config_dir = config_dir)
  providers <- config$providers %||% list()
  provider_names <- names(providers)

  if (length(provider_names) == 0L) {
    return(tibble::tibble(
      name = character(0),
      is_default = logical(0),
      provider_type = character(0),
      backend_preference = character(0),
      base_url = character(0),
      wire_api = character(0),
      default_model = character(0),
      review_model = character(0),
      reasoning_effort = character(0),
      history_enabled = logical(0),
      response_storage_disabled = logical(0),
      api_key_source = character(0)
    ))
  }

  tibble::tibble(
    name = provider_names,
    is_default = provider_names %in% (config$default_provider %||% ""),
    provider_type = vapply(provider_names, function(x) providers[[x]]$provider_type %||% NA_character_, character(1)),
    backend_preference = vapply(provider_names, function(x) providers[[x]]$backend_preference %||% "native", character(1)),
    base_url = vapply(provider_names, function(x) providers[[x]]$base_url %||% NA_character_, character(1)),
    wire_api = vapply(provider_names, function(x) providers[[x]]$wire_api %||% NA_character_, character(1)),
    default_model = vapply(provider_names, function(x) providers[[x]]$default_model %||% NA_character_, character(1)),
    review_model = vapply(provider_names, function(x) providers[[x]]$review_model %||% NA_character_, character(1)),
    reasoning_effort = vapply(provider_names, function(x) providers[[x]]$reasoning_effort %||% NA_character_, character(1)),
    history_enabled = vapply(provider_names, function(x) isTRUE(providers[[x]]$enable_history), logical(1)),
    response_storage_disabled = vapply(provider_names, function(x) isTRUE(providers[[x]]$disable_response_storage), logical(1)),
    api_key_source = vapply(
      provider_names,
      function(x) if (!is.null(providers[[x]]$api_key) && nzchar(providers[[x]]$api_key)) "config" else providers[[x]]$api_key_env %||% "unset",
      character(1)
    )
  )
}

#' Create an OpenAI-style provider for Shennong interpretation helpers
#'
#' @param api_key API key.
#' @param base_url Base URL of the API.
#' @param model Default model identifier used when no explicit model is
#'   supplied.
#' @param wire_api Upstream API style: \code{"responses"} or
#'   \code{"chat_completions"}.
#' @param reasoning_effort Optional reasoning effort sent with the request.
#' @param temperature Optional temperature. For the responses API this is only
#'   included when supplied.
#' @param disable_response_storage Logical; if \code{TRUE}, request that the
#'   upstream API does not retain the response.
#'
#' @return A provider function suitable for \code{sn_run_llm()}.
#'
#' @examples
#' \dontrun{
#' provider <- sn_make_openai_provider(
#'   api_key = Sys.getenv("OPENAI_API_KEY"),
#'   base_url = "https://api.openai.com/v1",
#'   model = "gpt-4.1",
#'   wire_api = "responses"
#' )
#' }
#' @export
sn_make_openai_provider <- function(api_key = Sys.getenv("OPENAI_API_KEY"),
                                    base_url = Sys.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1"),
                                    model = Sys.getenv("OPENAI_MODEL", "gpt-4.1"),
                                    wire_api = c("responses", "chat_completions"),
                                    reasoning_effort = NULL,
                                    temperature = NULL,
                                    disable_response_storage = TRUE) {
  wire_api <- match.arg(wire_api)

  if (!nzchar(api_key)) {
    stop("`api_key` is required. Set `OPENAI_API_KEY` or pass it explicitly.", call. = FALSE)
  }

  endpoint <- .sn_openai_endpoint(base_url = base_url, wire_api = wire_api)
  default_model <- model

  function(messages, model = NULL, ...) {
    resolved_model <- model %||% default_model
    extra_args <- list(...)

    payload <- if (identical(wire_api, "responses")) {
      base_payload <- list(
        model = resolved_model,
        input = lapply(messages, function(message) {
          list(
            role = message$role,
            content = list(list(
              type = "input_text",
              text = as.character(message$content %||% "")
            ))
          )
        })
      )
      if (!is.null(reasoning_effort) && nzchar(reasoning_effort)) {
        base_payload$reasoning <- list(effort = reasoning_effort)
      }
      if (!is.null(temperature)) {
        base_payload$temperature <- temperature
      }
      if (isTRUE(disable_response_storage)) {
        base_payload$store <- FALSE
      }
      modifyList(base_payload, extra_args)
    } else {
      base_payload <- list(
        model = resolved_model,
        messages = messages
      )
      if (!is.null(temperature)) {
        base_payload$temperature <- temperature
      }
      modifyList(base_payload, extra_args)
    }

    response <- curl::curl_fetch_memory(
      url = endpoint,
      handle = curl::new_handle(
        post = TRUE,
        postfields = jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"),
        httpheader = c(
          "Content-Type: application/json",
          paste("Authorization: Bearer", api_key),
          paste("x-api-key:", api_key)
        )
      )
    )

    parsed <- jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
    text <- .sn_text_scalar(.sn_extract_openai_response_text(parsed = parsed, wire_api = wire_api))
    if (is.null(text)) {
      stop(glue("OpenAI-style response from '{endpoint}' did not contain extractable text."), call. = FALSE)
    }

    list(
      text = text,
      model = parsed$model %||% resolved_model,
      raw = parsed
    )
  }
}

#' Create a Sub2API-compatible chat provider for Shennong interpretation helpers
#'
#' This helper returns a provider function suitable for \code{sn_run_llm()} and
#' the high-level \code{sn_interpret_*()} wrappers. It assumes an
#' OpenAI-compatible chat-completions endpoint such as Sub2API's proxy layer.
#'
#' @param api_key API key. Defaults to \code{Sys.getenv("SUB2API_API_KEY")}.
#' @param base_url Chat-completions endpoint URL. Defaults to
#'   \code{Sys.getenv("SUB2API_BASE_URL", "https://api.sub2api.com/v1/chat/completions")}.
#' @param model Default model identifier used when the caller does not supply
#'   one.
#' @param temperature Sampling temperature.
#'
#' @return A provider function that can be passed to \code{sn_run_llm()} or
#'   \code{sn_interpret_annotation()}.
#'
#' @examples
#' \dontrun{
#' provider <- sn_make_sub2api_provider()
#' response <- sn_run_llm(
#'   messages = list(list(role = "user", content = "Say hello.")),
#'   provider = provider
#' )
#' }
#' @export
sn_make_sub2api_provider <- function(api_key = Sys.getenv("SUB2API_API_KEY"),
                                     base_url = Sys.getenv("SUB2API_BASE_URL", "https://api.sub2api.com/v1/chat/completions"),
                                     model = Sys.getenv("SUB2API_MODEL", "gpt-4o-mini"),
                                     temperature = 0.2) {
  sn_make_openai_provider(
    api_key = api_key,
    base_url = base_url,
    model = model,
    wire_api = "chat_completions",
    temperature = temperature
  )
}

#' Create an \pkg{ellmer}-backed provider for Shennong interpretation helpers
#'
#' This adapter is useful when you want to manage transport, streaming, or
#' future structured-output features through \pkg{ellmer} while keeping
#' Shennong's interpretation API unchanged.
#'
#' @param api_key API key. Defaults to \code{Sys.getenv("OPENAI_API_KEY")}.
#' @param base_url Base URL of the API.
#' @param model Default model identifier.
#' @param provider_type One of \code{"openai_compatible"} or \code{"openai"}.
#' @param echo Echo mode forwarded to \pkg{ellmer}.
#' @param api_args Optional named list appended to each API request.
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
sn_make_ellmer_provider <- function(api_key = Sys.getenv("OPENAI_API_KEY"),
                                    base_url,
                                    model = NULL,
                                    provider_type = c("openai_compatible", "openai"),
                                    echo = c("none", "output", "all"),
                                    api_args = list()) {
  provider_type <- match.arg(provider_type)
  echo <- match.arg(echo)

  if (!requireNamespace("ellmer", quietly = TRUE)) {
    stop("Package `ellmer` is required for `sn_make_ellmer_provider()`. Install it first.", call. = FALSE)
  }
  if (!nzchar(api_key)) {
    stop("`api_key` is required. Set `OPENAI_API_KEY` or pass it explicitly.", call. = FALSE)
  }

  default_model <- model
  formatted_base <- .sn_format_openai_base_url(base_url)

  function(messages, model = NULL, ...) {
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
    chat <- if (identical(provider_type, "openai")) {
      ellmer::chat_openai(
        system_prompt = system_prompt,
        base_url = formatted_base,
        credentials = credentials,
        model = model %||% default_model,
        api_args = modifyList(api_args, list(...)),
        echo = echo
      )
    } else {
      ellmer::chat_openai_compatible(
        base_url = formatted_base,
        name = "Shennong provider",
        system_prompt = system_prompt,
        credentials = credentials,
        model = model %||% default_model,
        api_args = modifyList(api_args, list(...)),
        echo = echo
      )
    }

    text <- chat$chat(prompt_text)
    list(
      text = as.character(text),
      model = model %||% default_model,
      raw = NULL
    )
  }
}

#' Construct a provider from local Shennong LLM configuration
#'
#' @param name Optional provider name. Defaults to the configured default
#'   provider.
#' @param backend Backend to use: \code{"auto"}, \code{"native"}, or
#'   \code{"ellmer"}.
#' @param config_dir Local Shennong config root.
#'
#' @return A provider function suitable for \code{sn_run_llm()}.
#'
#' @examples
#' cfg_dir <- tempfile("shennong-config-")
#' sn_configure_llm_provider(
#'   name = "demo",
#'   base_url = "https://example.invalid",
#'   api_key = "demo-key",
#'   default_model = "gpt-4.1",
#'   config_dir = cfg_dir
#' )
#' provider <- sn_get_llm_provider(config_dir = cfg_dir)
#' @export
sn_get_llm_provider <- function(name = NULL,
                                backend = c("auto", "native", "ellmer"),
                                config_dir = "~/.shennong") {
  backend <- match.arg(backend)
  record <- .sn_resolve_llm_provider_record(name = name, config_dir = config_dir)

  resolved_backend <- if (identical(backend, "auto")) {
    preference <- record$backend_preference %||% "native"
    if (identical(preference, "ellmer") && requireNamespace("ellmer", quietly = TRUE)) {
      "ellmer"
    } else {
      "native"
    }
  } else {
    backend
  }

  api_key <- record$api_key %||% Sys.getenv(record$api_key_env %||% "OPENAI_API_KEY")
  if (!nzchar(api_key)) {
    stop(glue("Provider '{record$name}' does not have an API key in config and `{record$api_key_env %||% 'OPENAI_API_KEY'}` is empty."), call. = FALSE)
  }

  provider <- if (identical(resolved_backend, "ellmer")) {
    sn_make_ellmer_provider(
      api_key = api_key,
      base_url = record$base_url,
      model = record$default_model %||% NULL,
      provider_type = record$provider_type %||% "openai_compatible"
    )
  } else {
    sn_make_openai_provider(
      api_key = api_key,
      base_url = record$base_url,
      model = record$default_model %||% "gpt-4.1",
      wire_api = record$wire_api %||% "responses",
      reasoning_effort = record$reasoning_effort %||% NULL,
      disable_response_storage = isTRUE(record$disable_response_storage)
    )
  }

  .sn_wrap_provider_with_history(
    provider = provider,
    provider_name = record$name,
    base_url = record$base_url,
    wire_api = record$wire_api %||% "responses",
    backend = resolved_backend,
    config_dir = config_dir,
    history = isTRUE(record$enable_history)
  )
}

#' Test whether a configured LLM provider is reachable and usable
#'
#' @param name Optional provider name. Defaults to the configured default
#'   provider.
#' @param model Optional override model used for the test request.
#' @param prompt Prompt used for the connectivity check.
#' @param backend Backend to use: \code{"auto"}, \code{"native"}, or
#'   \code{"ellmer"}.
#' @param config_dir Local Shennong config root.
#' @param provider Optional provider function. When supplied, Shennong skips the
#'   config lookup and tests this function directly.
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
                                 backend = c("auto", "native", "ellmer"),
                                 config_dir = "~/.shennong",
                                 provider = NULL) {
  backend <- match.arg(backend)
  provider_name <- name %||% .sn_read_llm_config(config_dir = config_dir)$default_provider %||% "manual"
  provider <- provider %||% sn_get_llm_provider(
    name = name,
    backend = backend,
    config_dir = config_dir
  )

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
#' @param de_name Name of a stored marker result.
#' @param cluster_col Metadata column containing cluster labels.
#' @param n_markers Number of top markers per cluster.
#' @param enrichment_name Optional stored enrichment result used to add
#'   cluster-level functional evidence.
#' @param n_terms Number of enrichment terms per cluster when
#'   \code{enrichment_name} is supplied.
#' @param include_qc Logical; whether to include cluster-level QC summaries in
#'   the evidence bundle.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function. When left \code{NULL},
#'   Shennong will try to use the default provider configured under
#'   \code{~/.shennong}.
#' @param model Optional model identifier.
#' @param include_json_schema Logical; whether to request structured JSON output
#'   from the provider. Defaults to \code{TRUE} for annotation workflows.
#' @param apply_metadata Logical; if \code{TRUE} and a structured annotation
#'   response is returned, map the cluster labels back onto each cell in the
#'   Seurat metadata.
#' @param metadata_prefix Prefix used for metadata columns written back to the
#'   Seurat object when \code{apply_metadata = TRUE}.
#' @param label_style Naming style used to normalize returned cell-type labels.
#'   One of \code{"title"}, \code{"snake"}, or \code{"asis"}.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param store_name Name used under \code{object@misc$interpretation_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
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
                                    de_name,
                                    cluster_col = "seurat_clusters",
                                    n_markers = 10,
                                    enrichment_name = NULL,
                                    n_terms = 5,
                                    include_qc = TRUE,
                                    background = NULL,
                                    output_format = c("llm", "human"),
                                    provider = NULL,
                                    model = NULL,
                                    include_json_schema = TRUE,
                                    apply_metadata = TRUE,
                                    metadata_prefix = "sn_annotation",
                                    label_style = c("title", "snake", "asis"),
                                    return_prompt = FALSE,
                                    store_name = "default",
                                    return_object = TRUE,
                                    ...) {
  output_format <- match.arg(output_format)
  label_style <- match.arg(label_style)
  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = de_name,
    cluster_col = cluster_col,
    n_markers = n_markers,
    enrichment_name = enrichment_name,
    n_terms = n_terms,
    include_qc = include_qc
  )
  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "annotation",
    audience = "scientist",
    style = "cell type annotation",
    background = background,
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
    label_style = label_style,
    apply_metadata = apply_metadata,
    return_prompt = return_prompt,
    return_object = return_object,
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
                            ...) {
  output_format <- match.arg(output_format)
  evidence <- sn_prepare_de_evidence(object = object, de_name = de_name, n_genes = n_genes)
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
                                    ...) {
  output_format <- match.arg(output_format)
  evidence <- sn_prepare_enrichment_evidence(
    object = object,
    enrichment_name = enrichment_name,
    n_terms = n_terms
  )
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
                             ...) {
  output_format <- match.arg(output_format)
  evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_name = cluster_de_name,
    contrast_de_name = contrast_de_name,
    enrichment_name = enrichment_name,
    cluster_col = cluster_col
  )
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
                                   ...) {
  output_format <- match.arg(output_format)
  evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_name = cluster_de_name,
    enrichment_name = enrichment_name,
    cluster_col = cluster_col
  )
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
                                          ...) {
  output_format <- match.arg(output_format)
  evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_name = cluster_de_name,
    contrast_de_name = contrast_de_name,
    enrichment_name = enrichment_name,
    cluster_col = cluster_col
  )
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
    ...
  )
}
