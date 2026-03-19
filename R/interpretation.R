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

.sn_build_messages <- function(system_prompt, user_prompt) {
  list(
    list(role = "system", content = system_prompt),
    list(role = "user", content = user_prompt)
  )
}

.sn_finish_interpretation <- function(object,
                                      task,
                                      evidence,
                                      prompt,
                                      provider = NULL,
                                      model = NULL,
                                      store_name = "default",
                                      return_prompt = FALSE,
                                      return_object = TRUE,
                                      ...) {
  if (return_prompt) {
    return(prompt)
  }

  if (is_null(provider)) {
    stop("`provider` must be supplied unless `return_prompt = TRUE`.")
  }

  response <- sn_run_llm(
    messages = prompt$messages,
    provider = provider,
    model = model,
    ...
  )

  interpretation_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    task = task,
    evidence = evidence,
    prompt = prompt,
    response = response,
    model_info = list(model = model)
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
#'   rownames(counts) <- c("CD3D", "CD3E", "TRAC", "LTB", "MS4A1", "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1")
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
#'   counts[c("CD3D", "CD3E", "TRAC"), 1:12] <- counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
#'   counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <- counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   evidence <- sn_prepare_annotation_evidence(obj, de_name = "celltype_markers", cluster_col = "cell_type")
#'   names(evidence)
#' }
#' @export
sn_prepare_annotation_evidence <- function(object,
                                           de_name,
                                           cluster_col = "seurat_clusters",
                                           n_markers = 10) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  de_result <- .sn_get_misc_result(object = object, collection = "de_results", store_name = de_name)
  cluster_summary <- .sn_prepare_cluster_summary(object = object, cluster_col = cluster_col)
  marker_summary <- .sn_prepare_marker_summary(de_result = de_result, n_markers = n_markers)
  prediction_summary <- .sn_prepare_prediction_summary(object = object, cluster_col = cluster_col)

  merged_summary <- dplyr::left_join(cluster_summary, marker_summary, by = "cluster")
  if (nrow(prediction_summary) > 0) {
    merged_summary <- dplyr::left_join(merged_summary, prediction_summary, by = "cluster")
  }

  list(
    task = "annotation",
    cluster_col = cluster_col,
    source_de_name = de_name,
    analysis_method = de_result$method,
    species = tryCatch(sn_get_species(object), error = function(...) NULL),
    cluster_summary = merged_summary,
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
#'   counts[c("CD3D", "CD3E", "TRAC"), 1:12] <- counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
#'   counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <- counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
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
    return(list(
      task = "de",
      source_de_name = de_name,
      summary = summary,
      top_markers = top_table,
      caveats = character()
    ))
  }

  ordered <- result_table[order(-abs(result_table[[rank_col]])), , drop = FALSE]
  top_hits <- utils::head(ordered, n_genes)

  list(
    task = "de",
    source_de_name = de_name,
    summary = summary,
    top_hits = tibble::as_tibble(top_hits),
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
#'   rownames(counts) <- c("CD3D", "CD3E", "TRAC", "LTB", "MS4A1", "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1")
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
#'   evidence <- sn_prepare_results_evidence(obj, cluster_de_name = "celltype_markers", enrichment_name = "demo_gsea", cluster_col = "cell_type")
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
#' @param include_json_schema Whether to request structured JSON output.
#'
#' @return A prompt bundle with \code{system}, \code{user}, and \code{messages}.
#'
#' @examples
#' evidence <- list(task = "annotation", cluster_summary = data.frame(cluster = "0", top_markers = "CD3D, TRAC"))
#' prompt <- sn_build_prompt(evidence = evidence, task = "annotation")
#' names(prompt)
#' @export
sn_build_prompt <- function(evidence,
                            task = c("annotation", "de", "enrichment", "results", "figure_legend", "presentation_summary"),
                            style = NULL,
                            audience = c("scientist", "clinician", "general"),
                            language = c("en", "zh"),
                            include_json_schema = FALSE) {
  task <- match.arg(task)
  audience <- match.arg(audience)
  language <- match.arg(language)

  system_prompt <- paste(
    "You are assisting with interpretation of single-cell transcriptomics analyses produced by the Shennong R package.",
    "Use only the supplied evidence. Distinguish direct evidence from inference.",
    "Keep explanations biologically grounded and concise."
  )

  style_line <- if (!is_null(style)) paste0("Target style: ", style, ".") else ""
  json_line <- if (isTRUE(include_json_schema)) "Return a structured JSON object followed by a brief narrative explanation." else ""

  user_prompt <- paste(
    paste0("Task: ", task, "."),
    paste0("Audience: ", audience, "."),
    paste0("Language: ", language, "."),
    style_line,
    json_line,
    "Evidence:",
    .sn_render_prompt_value(evidence),
    sep = "\n\n"
  )

  list(
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
    return(response)
  }

  stop("`provider` must return either a single string or a list containing `text`.")
}

#' Interpret cluster markers for cell-type annotation
#'
#' @param object A \code{Seurat} object.
#' @param de_name Name of a stored marker result.
#' @param cluster_col Metadata column containing cluster labels.
#' @param n_markers Number of top markers per cluster.
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
#'   counts[c("CD3D", "CD3E", "TRAC"), 1:12] <- counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
#'   counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <- counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_interpret_annotation(obj, de_name = "celltype_markers", cluster_col = "cell_type", return_prompt = TRUE)
#'   prompt$task
#' }
#' @export
sn_interpret_annotation <- function(object,
                                    de_name,
                                    cluster_col = "seurat_clusters",
                                    n_markers = 10,
                                    provider = NULL,
                                    model = NULL,
                                    return_prompt = FALSE,
                                    store_name = "default",
                                    return_object = TRUE,
                                    ...) {
  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = de_name,
    cluster_col = cluster_col,
    n_markers = n_markers
  )
  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "annotation",
    audience = "scientist"
  )

  .sn_finish_interpretation(
    object = object,
    task = "interpret_annotation",
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

#' Interpret a stored differential-expression result
#'
#' @param object A \code{Seurat} object.
#' @param de_name Name of a stored DE result.
#' @param n_genes Number of top genes to retain.
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
#'   counts[c("CD3D", "CD3E", "TRAC"), 1:12] <- counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
#'   counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <- counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
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
                            provider = NULL,
                            model = NULL,
                            return_prompt = FALSE,
                            store_name = "default",
                            return_object = TRUE,
                            ...) {
  evidence <- sn_prepare_de_evidence(object = object, de_name = de_name, n_genes = n_genes)
  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "de",
    audience = "scientist"
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
#'   rownames(counts) <- c("CD3D", "CD3E", "TRAC", "LTB", "MS4A1", "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1")
#'   colnames(counts) <- paste0("cell", 1:12)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj <- sn_store_enrichment(
#'     obj,
#'     tibble::tibble(ID = "GO:0001", Description = "immune response", NES = 2, p.adjust = 0.01),
#'     store_name = "demo_gsea"
#'   )
#'   prompt <- sn_interpret_enrichment(obj, enrichment_name = "demo_gsea", return_prompt = TRUE)
#'   prompt$task
#' }
#' @export
sn_interpret_enrichment <- function(object,
                                    enrichment_name,
                                    n_terms = 10,
                                    provider = NULL,
                                    model = NULL,
                                    return_prompt = FALSE,
                                    store_name = "default",
                                    return_object = TRUE,
                                    ...) {
  evidence <- sn_prepare_enrichment_evidence(
    object = object,
    enrichment_name = enrichment_name,
    n_terms = n_terms
  )
  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "enrichment",
    audience = "scientist"
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
#'   rownames(counts) <- c("CD3D", "CD3E", "TRAC", "LTB", "MS4A1", "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1")
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
#'   prompt <- sn_write_results(obj, cluster_de_name = "celltype_markers", enrichment_name = "demo_gsea", cluster_col = "cell_type", return_prompt = TRUE)
#'   prompt$task
#' }
#' @export
sn_write_results <- function(object,
                             cluster_de_name = NULL,
                             contrast_de_name = NULL,
                             enrichment_name = NULL,
                             cluster_col = "seurat_clusters",
                             provider = NULL,
                             model = NULL,
                             return_prompt = FALSE,
                             store_name = "default",
                             return_object = TRUE,
                             ...) {
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
    audience = "scientist"
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
#'   rownames(counts) <- c("CD3D", "CD3E", "TRAC", "LTB", "MS4A1", "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1")
#'   colnames(counts) <- paste0("cell", 1:24)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_write_figure_legend(obj, cluster_de_name = "celltype_markers", cluster_col = "cell_type", return_prompt = TRUE)
#'   prompt$task
#' }
#' @export
sn_write_figure_legend <- function(object,
                                   cluster_de_name = NULL,
                                   enrichment_name = NULL,
                                   cluster_col = "seurat_clusters",
                                   provider = NULL,
                                   model = NULL,
                                   return_prompt = FALSE,
                                   store_name = "default",
                                   return_object = TRUE,
                                   ...) {
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
    audience = "scientist"
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
#'   rownames(counts) <- c("CD3D", "CD3E", "TRAC", "LTB", "MS4A1", "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1")
#'   colnames(counts) <- paste0("cell", 1:24)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
#'     layer = "data", min_pct = 0, logfc_threshold = 0,
#'     store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_write_presentation_summary(obj, cluster_de_name = "celltype_markers", cluster_col = "cell_type", return_prompt = TRUE)
#'   prompt$task
#' }
#' @export
sn_write_presentation_summary <- function(object,
                                          cluster_de_name = NULL,
                                          contrast_de_name = NULL,
                                          enrichment_name = NULL,
                                          cluster_col = "seurat_clusters",
                                          provider = NULL,
                                          model = NULL,
                                          return_prompt = FALSE,
                                          store_name = "default",
                                          return_object = TRUE,
                                          ...) {
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
    audience = "scientist"
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
