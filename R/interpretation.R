.sn_misc_result_registry <- function() {
  artifact_collections <- c(
    "sn_run_cluster",
    "integration",
    "integration_comparison",
    "mmochi",
    "bpcells_layers",
    "infercnvpy",
    "input_source",
    "qc",
    "hvg_selection",
    "rare_feature_selection",
    "label_transfer",
    "label_transfer_reference",
    "scdesign3",
    "coralysis",
    "scarches",
    "scpoli",
    "cellphonedb",
    "cell2location",
    "tangram",
    "squidpy",
    "spatialdata",
    "stlearn"
  )
  artifact_registry <- tibble::tibble(
    collection = artifact_collections,
    type = paste0(artifact_collections, "_artifact"),
    schema_version = rep(NA_character_, length(artifact_collections)),
    required_fields = I(rep(list(character()), length(artifact_collections))),
    contract_scope = rep("artifact", length(artifact_collections)),
    listable = rep(FALSE, length(artifact_collections)),
    table_required = rep(FALSE, length(artifact_collections)),
    reader = rep(NA_character_, length(artifact_collections)),
    writer = rep(NA_character_, length(artifact_collections))
  )
  artifact_registry$type[artifact_registry$collection == "sn_run_cluster"] <-
    "clustering_stage_cache"
  artifact_registry$type[artifact_registry$collection == "input_source"] <-
    "input_source"
  artifact_registry
}


.sn_resolve_stored_result_id <- function(object,
                                         type,
                                         result_id = NULL,
                                         preferred_analysis = NULL) {
  collection_data <- .sn_result_store(object)[[type]] %||% list()

  if (!is_null(result_id) && nzchar(result_id)) {
    if (!result_id %in% names(collection_data)) {
      stop(glue("No result with `result_id = \"{result_id}\"` was found for analysis type '{type}'."))
    }
    return(result_id)
  }

  if (length(collection_data) == 0L) {
    stop(glue("No stored results were found for analysis type '{type}'."), call. = FALSE)
  }

  available_names <- names(collection_data)
  latest_name <- function(candidates) {
    if (length(candidates) == 0L) {
      return(NULL)
    }
    created_at <- vapply(
      candidates,
      function(candidate) {
        collection_data[[candidate]]$provenance$timestamp %||% ""
      },
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

  .sn_log_info("`result_id` was not supplied; using '{resolved_name}' for analysis type '{type}'.")
  resolved_name
}

.sn_result_n_rows <- function(result) {
  tables <- result[["tables"]] %||% list()
  if (is.data.frame(tables[["primary"]])) {
    return(nrow(tables[["primary"]]))
  }
  0L
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


.sn_as_enrichment_table <- function(result) {
  table <- if (is.data.frame(result)) {
    tibble::as_tibble(result)
  } else {
    tibble::as_tibble(as.data.frame(result))
  }
  if (ncol(table) == 0L) {
    # A no-hit enrichment run is still a valid analytical result. Preserve a
    # typed, discoverable zero-row schema instead of an unstructured tibble()
    # that cannot satisfy the unified result contract.
    table <- tibble::tibble(
      ID = character(),
      Description = character(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalue = numeric()
    )
  }
  table
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
                                      result_id = "default",
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
    schema_version = .sn_analysis_result_schema_version(),
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

  object <- sn_store_result(
    object = object,
    type = "interpretation",
    result_id = result_id,
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
                                          result_id = "default",
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
    schema_version = .sn_analysis_result_schema_version(),
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

  object <- sn_store_result(
    object = object,
    type = "interpretation",
    result_id = result_id,
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

#' Interpret cluster markers for cell-type annotation
#'
#' @param object A \code{Seurat} object.
#' @param de_result_id Optional stored marker-result name. When omitted, Shennong
#'   prefers \code{"default"}, then a single available result, and otherwise
#'   the most recent marker result.
#' @param cluster_by Metadata column containing cluster labels.
#' @param n_markers Number of top markers per cluster.
#' @param marker_selection How to choose marker genes for annotation evidence:
#'   \code{"specific"} prefers genes that are relatively unique to one cluster,
#'   while \code{"top"} keeps the raw top-ranked genes.
#' @param enrichment_result_id Optional stored enrichment result used to add
#'   cluster-level functional evidence.
#' @param n_terms Number of enrichment terms per cluster when
#'   \code{enrichment_result_id} is supplied.
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
#' @param result_id Stable identifier for the stored interpretation result.
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
#'     result_id = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_interpret_annotation(
#'     obj,
#'     de_result_id = "celltype_markers",
#'     cluster_by = "cell_type",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_interpret_annotation <- function(object,
                                    de_result_id = NULL,
                                    cluster_by = NULL,
                                    n_markers = 10,
                                    marker_selection = c("specific", "top"),
                                    enrichment_result_id = NULL,
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
                                    result_id = "default",
                                    return_object = TRUE,
                                    show_progress = interactive(),
                                    ...) {
  result_id <- .sn_validate_result_id(result_id)
  cluster_by <- cluster_by %||% "seurat_clusters"
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
    de_result_id = de_result_id,
    cluster_by = cluster_by,
    n_markers = n_markers,
    marker_selection = marker_selection,
    enrichment_result_id = enrichment_result_id,
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
      result_id = result_id,
      cluster_col = cluster_by,
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
    result_id = result_id,
    cluster_col = cluster_by,
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
#' @param de_result_id Name of a stored DE result.
#' @param n_genes Number of top genes to retain.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param result_id Stable identifier for the stored interpretation result.
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
#'     result_id = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_interpret_de(obj, de_result_id = "celltype_markers", return_prompt = TRUE)
#'   prompt$task
#' }
#' @export
sn_interpret_de <- function(object,
                            de_result_id,
                            n_genes = 15,
                            background = NULL,
                            output_format = c("llm", "human"),
                            provider = NULL,
                            model = NULL,
                            return_prompt = FALSE,
                            result_id = "default",
                            return_object = TRUE,
                            show_progress = interactive(),
                            ...) {
  result_id <- .sn_validate_result_id(result_id)
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "interpret_de",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing DE evidence")
  evidence <- sn_prepare_de_evidence(object = object, de_result_id = de_result_id, n_genes = n_genes)
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
    result_id = result_id,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}

#' Interpret a stored enrichment result
#'
#' @param object A \code{Seurat} object.
#' @param enrichment_result_id Name of a stored enrichment result.
#' @param n_terms Number of top enrichment terms to retain.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param result_id Stable identifier for the stored interpretation result.
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
#'     result_id = "demo_gsea"
#'   )
#'   prompt <- sn_interpret_enrichment(
#'     obj,
#'     enrichment_result_id = "demo_gsea",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_interpret_enrichment <- function(object,
                                    enrichment_result_id,
                                    n_terms = 10,
                                    background = NULL,
                                    output_format = c("llm", "human"),
                                    provider = NULL,
                                    model = NULL,
                                    return_prompt = FALSE,
                                    result_id = "default",
                                    return_object = TRUE,
                                    show_progress = interactive(),
                                    ...) {
  result_id <- .sn_validate_result_id(result_id)
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "interpret_enrichment",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing enrichment evidence")
  evidence <- sn_prepare_enrichment_evidence(
    object = object,
    enrichment_result_id = enrichment_result_id,
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
    result_id = result_id,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}

#' Write a manuscript-style results summary from stored analysis outputs
#'
#' @param object A \code{Seurat} object.
#' @param cluster_de_result_id Optional stored cluster-marker result.
#' @param contrast_de_result_id Optional stored contrast or pseudobulk result.
#' @param enrichment_result_id Optional stored enrichment result.
#' @param cluster_by Metadata column containing cluster labels.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param result_id Stable identifier for the stored interpretation result.
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
#'     result_id = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   obj <- sn_store_enrichment(
#'     obj,
#'     tibble::tibble(ID = "GO:0001", Description = "immune response", NES = 2, p.adjust = 0.01),
#'     result_id = "demo_gsea"
#'   )
#'   prompt <- sn_write_results(
#'     obj,
#'     cluster_de_result_id = "celltype_markers",
#'     enrichment_result_id = "demo_gsea",
#'     cluster_by = "cell_type",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_write_results <- function(object,
                             cluster_de_result_id = NULL,
                             contrast_de_result_id = NULL,
                             enrichment_result_id = NULL,
                             cluster_by = NULL,
                             background = NULL,
                             output_format = c("llm", "human"),
                             provider = NULL,
                             model = NULL,
                             return_prompt = FALSE,
                             result_id = "default",
                             return_object = TRUE,
                             show_progress = interactive(),
                             ...) {
  result_id <- .sn_validate_result_id(result_id)
  cluster_by <- cluster_by %||% "seurat_clusters"
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "write_results",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing results evidence")
  evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_result_id = cluster_de_result_id,
    contrast_de_result_id = contrast_de_result_id,
    enrichment_result_id = enrichment_result_id,
    cluster_by = cluster_by
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
    result_id = result_id,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}

#' Write a figure legend from stored analysis outputs
#'
#' @param object A \code{Seurat} object.
#' @param cluster_de_result_id Optional stored cluster-marker result.
#' @param enrichment_result_id Optional stored enrichment result.
#' @param cluster_by Metadata column containing cluster labels.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param result_id Stable identifier for the stored interpretation result.
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
#'     result_id = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_write_figure_legend(
#'     obj,
#'     cluster_de_result_id = "celltype_markers",
#'     cluster_by = "cell_type",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_write_figure_legend <- function(object,
                                   cluster_de_result_id = NULL,
                                   enrichment_result_id = NULL,
                                   cluster_by = NULL,
                                   background = NULL,
                                   output_format = c("llm", "human"),
                                   provider = NULL,
                                   model = NULL,
                                   return_prompt = FALSE,
                                   result_id = "default",
                                   return_object = TRUE,
                                   show_progress = interactive(),
                                   ...) {
  result_id <- .sn_validate_result_id(result_id)
  cluster_by <- cluster_by %||% "seurat_clusters"
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "write_figure_legend",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing legend evidence")
  evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_result_id = cluster_de_result_id,
    enrichment_result_id = enrichment_result_id,
    cluster_by = cluster_by
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
    result_id = result_id,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}

#' Write a presentation-style summary from stored analysis outputs
#'
#' @param object A \code{Seurat} object.
#' @param cluster_de_result_id Optional stored cluster-marker result.
#' @param contrast_de_result_id Optional stored contrast or pseudobulk result.
#' @param enrichment_result_id Optional stored enrichment result.
#' @param cluster_by Metadata column containing cluster labels.
#' @param background Optional study-specific background information to provide
#'   additional interpretation context.
#' @param output_format One of \code{"llm"} for a model-ready prompt bundle or
#'   \code{"human"} for a human-readable summary.
#' @param provider Optional model provider function.
#' @param model Optional model identifier.
#' @param return_prompt If \code{TRUE}, return the prompt bundle without calling
#'   the provider.
#' @param result_id Stable identifier for the stored interpretation result.
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
#'     result_id = "celltype_markers", return_object = TRUE, verbose = FALSE
#'   )
#'   prompt <- sn_write_presentation_summary(
#'     obj,
#'     cluster_de_result_id = "celltype_markers",
#'     cluster_by = "cell_type",
#'     return_prompt = TRUE
#'   )
#'   prompt$task
#' }
#' @export
sn_write_presentation_summary <- function(object,
                                          cluster_de_result_id = NULL,
                                          contrast_de_result_id = NULL,
                                          enrichment_result_id = NULL,
                                          cluster_by = NULL,
                                          background = NULL,
                                          output_format = c("llm", "human"),
                                          provider = NULL,
                                          model = NULL,
                                          return_prompt = FALSE,
                                          result_id = "default",
                                          return_object = TRUE,
                                          show_progress = interactive(),
                                          ...) {
  result_id <- .sn_validate_result_id(result_id)
  cluster_by <- cluster_by %||% "seurat_clusters"
  output_format <- match.arg(output_format)
  progress_state <- .sn_interpret_progress_start(
    task = "write_presentation_summary",
    enabled = show_progress,
    total_steps = 4L
  )
  progress_state <- .sn_interpret_progress_step(progress_state, "Preparing presentation evidence")
  evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_result_id = cluster_de_result_id,
    contrast_de_result_id = contrast_de_result_id,
    enrichment_result_id = enrichment_result_id,
    cluster_by = cluster_by
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
    result_id = result_id,
    return_prompt = return_prompt,
    return_object = return_object,
    progress_state = progress_state,
    ...
  )
}
