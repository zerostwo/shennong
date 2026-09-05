.sn_plot_method_registry <- function() {
  list(
    annotation = list(
      default = "confidence_cluster",
      views = c("confidence_cluster", "confidence_cell", "confusion", "projection"),
      input = "Seurat or annotation result"
    ),
    differential_abundance = list(
      default = "effect", views = "effect",
      input = "Seurat or differential-abundance result"
    ),
    de = list(
      default = "volcano", views = c("volcano", "ma", "effect", "heatmap"),
      input = "Seurat, DE result, or compatible table"
    ),
    enrichment = list(
      default = "dot", views = c("dot", "bar", "ridge", "network", "emap", "gsea"),
      input = "Seurat, enrichment result, compatible table, or clusterProfiler result"
    ),
    bulk_qc = list(
      default = "library_size",
      views = c(
        "library_size", "detected_features", "mean_correlation", "pca",
        "correlation_heatmap", "sample_scatter"
      ),
      input = "Seurat or bulk-QC result"
    ),
    bulk_de = list(
      default = "volcano", views = "volcano",
      input = "Seurat or bulk-DE result"
    ),
    bulk_network = list(
      default = "traits", views = c("traits", "modules"),
      input = "Seurat or bulk-network result"
    ),
    bulk_survival = list(
      default = "forest",
      views = c("forest", "km", "risk_table", "ph", "ph_test", "cumulative_hazard"),
      input = "Seurat or bulk-survival result"
    ),
    cell_communication = list(
      default = "bubble",
      views = c("bubble", "heatmap", "network", "chord", "river", "ligand_target", "comparison"),
      input = "Seurat, cell-communication result, or compatible table"
    ),
    cnv = list(
      default = "heatmap", views = c("heatmap", "umap", "score", "sample", "association"),
      input = "Seurat or CNV result"
    ),
    grn = list(
      default = "network", views = c("network", "activity", "specificity"),
      input = "Seurat or GRN result"
    ),
    metabolism = list(
      default = "activity", views = c("activity", "heatmap", "differential", "sample"),
      input = "Seurat or metabolism result"
    ),
    program_discovery = list(
      default = "weights", views = c("weights", "activity", "stability"),
      input = "Seurat or program-discovery result"
    ),
    program_scoring = list(
      default = "activity", views = c("activity", "heatmap"),
      input = "Seurat or program-scoring result"
    ),
    state_priority = list(
      default = "ranking", views = "ranking",
      input = "Seurat or state-priority result"
    ),
    scissor = list(
      default = "states", views = c("states", "cells", "samples", "correlations", "reliability"),
      input = "Seurat or Scissor result"
    ),
    trajectory = list(
      default = "embedding",
      views = c(
        "embedding", "pseudotime", "lineage_probability", "dynamic_heatmap",
        "gene_trend", "branch_comparison"
      ),
      input = "Seurat or trajectory result"
    ),
    velocity = list(
      default = "embedding", views = "embedding",
      input = "Seurat or velocity result"
    ),
    fate = list(
      default = "probability", views = "probability",
      input = "Seurat or fate result"
    ),
    spatial_domains = list(
      default = "domain", views = "domain",
      input = "Seurat or spatial-domain result"
    ),
    spatial_features = list(
      default = "ranked", views = "ranked",
      input = "Seurat or spatial-feature result"
    ),
    spatial_neighborhood = list(
      default = "enrichment", views = c("enrichment", "cooccurrence"),
      input = "Seurat or spatial-neighborhood result"
    ),
    spatial_communication = list(
      default = "interactions", views = "interactions",
      input = "Seurat or spatial-communication result"
    ),
    qc = list(
      default = "qc_score", views = c("qc_score", "n_cells", "retention_fraction"),
      input = "QC assessment result or by-sample table"
    ),
    resolution_sweep = list(
      default = "summary", views = "summary",
      input = "resolution-sweep result or summary table"
    ),
    integration_assessment = list(
      default = "summary", views = "summary",
      input = "integration assessment or summary table"
    )
  )
}

.sn_normalize_plot_analysis_type <- function(analysis_type) {
  if (is_null(analysis_type)) return(NULL)
  if (!is.character(analysis_type) || length(analysis_type) != 1L ||
      is.na(analysis_type) || !nzchar(analysis_type)) {
    stop("`analysis_type` must be one non-empty name.", call. = FALSE)
  }
  aliases <- c(
    abundance = "differential_abundance",
    communication = "cell_communication",
    differential_expression = "de",
    wgcna = "bulk_network",
    survival = "bulk_survival",
    program = "program_scoring",
    spatial_domain = "spatial_domains",
    spatial_svg = "spatial_features",
    integration = "integration_assessment"
  )
  normalized <- tolower(analysis_type)
  if (normalized %in% names(aliases)) unname(aliases[[normalized]]) else normalized
}

.sn_plot_registry_table <- function() {
  registry <- .sn_plot_method_registry()
  dplyr::bind_rows(lapply(names(registry), function(analysis_type) {
    entry <- registry[[analysis_type]]
    required <- unname(vapply(entry$views, function(view) {
      switch(
        paste(analysis_type, view, sep = "/"),
        "annotation/confusion" = "truth",
        "annotation/projection" = "Seurat input",
        "trajectory/lineage_probability" = "lineage",
        "trajectory/gene_trend" = "features",
        ""
      )
    }, character(1)))
    tibble::tibble(
      analysis_type = analysis_type,
      view = entry$views,
      default = entry$views == entry$default,
      required_parameters = required,
      accepted_input = entry$input
    )
  }))
}

#' List result-aware plot methods
#'
#' Returns the canonical analysis-type and view vocabulary accepted by
#' [sn_plot_result()]. The table is intended for discovery and programmatic UI
#' generation; one row represents one valid analysis/view pair.
#'
#' @param analysis_type Optional analysis type or documented alias used to
#'   filter the returned table.
#'
#' @return A tibble with `analysis_type`, `view`, `default`,
#'   `required_parameters`, and `accepted_input` columns.
#'
#' @examples
#' sn_list_plot_methods("de")
#' @export
sn_list_plot_methods <- function(analysis_type = NULL) {
  table <- .sn_plot_registry_table()
  if (!is_null(analysis_type)) {
    analysis_type <- .sn_normalize_plot_analysis_type(analysis_type)
    table <- table[table$analysis_type == analysis_type, , drop = FALSE]
    if (nrow(table) == 0L) {
      stop("Unknown plot analysis type '", analysis_type, "'.", call. = FALSE)
    }
  }
  table
}

.sn_infer_stored_plot_type <- function(object, analysis_type, result_id) {
  store <- .sn_result_store(object)
  available_types <- names(store)[vapply(store, length, integer(1)) > 0L]
  if (!is_null(analysis_type)) return(analysis_type)
  candidates <- if (!is_null(result_id)) {
    available_types[vapply(store[available_types], function(results) {
      result_id %in% names(results)
    }, logical(1))]
  } else if (length(available_types) == 1L) {
    available_types
  } else {
    character()
  }
  if (length(candidates) == 1L) return(candidates[[1L]])
  detail <- if (length(available_types) > 0L) {
    paste0(" Available types: ", paste(available_types, collapse = ", "), ".")
  } else {
    " The object contains no stored Shennong results."
  }
  stop(
    "`analysis_type` is required because the stored result type is ambiguous.",
    detail,
    call. = FALSE
  )
}

.sn_select_stored_plot_id <- function(object, analysis_type, result_id) {
  results <- .sn_result_store(object)[[analysis_type]] %||% list()
  ids <- names(results)
  if (length(ids) == 0L) {
    stop("No stored results were found for analysis type '", analysis_type, "'.", call. = FALSE)
  }
  if (!is_null(result_id)) {
    .sn_validate_result_id(result_id)
    if (!result_id %in% ids) {
      stop(
        "No result with `result_id = \"", result_id, "\"` was found for analysis type '",
        analysis_type, "'. Available IDs: ", paste(ids, collapse = ", "), ".",
        call. = FALSE
      )
    }
    return(result_id)
  }
  if ("default" %in% ids) return("default")
  if (length(ids) == 1L) return(ids[[1L]])
  stop(
    "`result_id` is required because multiple '", analysis_type,
    "' results are stored. Available IDs: ", paste(ids, collapse = ", "), ".",
    call. = FALSE
  )
}

.sn_resolve_canonical_plot_input <- function(object, analysis_type, result_id) {
  analysis_type <- .sn_normalize_plot_analysis_type(analysis_type)
  if (inherits(object, "Seurat")) {
    analysis_type <- .sn_infer_stored_plot_type(object, analysis_type, result_id)
    result_id <- .sn_select_stored_plot_id(object, analysis_type, result_id)
    result <- sn_get_result(object, analysis_type, result_id)
    return(list(
      result = result, analysis_type = analysis_type, result_id = result_id,
      source_object = object
    ))
  }
  if (is.list(object) && !is.data.frame(object) &&
      is.character(object$analysis_type) && length(object$analysis_type) == 1L) {
    sn_validate_result(object)
    inferred <- .sn_normalize_plot_analysis_type(object$analysis_type)
    if (!is_null(analysis_type) && !identical(analysis_type, inferred)) {
      stop(
        "`analysis_type = \"", analysis_type, "\"` does not match result analysis type '",
        inferred, "'.",
        call. = FALSE
      )
    }
    if (!is_null(result_id) && !identical(result_id, object$result_id)) {
      stop(
        "`result_id = \"", result_id, "\"` does not match direct result ID '",
        object$result_id, "'.",
        call. = FALSE
      )
    }
    return(list(
      result = object, analysis_type = inferred,
      result_id = object$result_id %||% result_id, source_object = NULL
    ))
  }
  if (is_null(analysis_type)) {
    stop(
      "`analysis_type` is required for a table or legacy assessment input.",
      call. = FALSE
    )
  }
  list(result = object, analysis_type = analysis_type, result_id = result_id, source_object = NULL)
}

.sn_plot_call <- function(fun, fixed, dots) {
  duplicates <- intersect(names(fixed), names(dots))
  if (length(duplicates) > 0L) {
    stop(
      "Argument(s) controlled by `sn_plot_result()` were also supplied in `...`: ",
      paste(duplicates, collapse = ", "), ". Use `view` for the plot mode.",
      call. = FALSE
    )
  }
  do.call(fun, c(fixed, dots))
}

.sn_dispatch_plot_result <- function(input, view, dots) {
  result <- input$result
  analysis_type <- input$analysis_type
  result_id <- input$result_id
  source_object <- input$source_object
  call <- function(fun, ...) .sn_plot_call(fun, list(...), dots)

  if (identical(analysis_type, "annotation")) {
    if (identical(view, "confidence_cluster")) return(call(sn_plot_annotation_confidence, x = result, level = "cluster"))
    if (identical(view, "confidence_cell")) return(call(sn_plot_annotation_confidence, x = result, level = "cell"))
    if (identical(view, "confusion")) return(call(sn_plot_annotation_confusion, x = result))
    if (is_null(source_object)) stop("The annotation projection view requires Seurat input.", call. = FALSE)
    return(call(sn_plot_reference_projection, object = source_object, result_id = result_id))
  }
  if (identical(analysis_type, "differential_abundance")) return(call(sn_plot_abundance, x = result))
  if (identical(analysis_type, "de")) return(call(sn_plot_de, result = result, type = view))
  if (identical(analysis_type, "enrichment")) {
    if (identical(view, "gsea")) return(call(sn_plot_gsea, result = result))
    return(call(sn_plot_enrichment, result = result, type = view))
  }
  if (identical(analysis_type, "bulk_qc")) {
    if (view %in% c("library_size", "detected_features", "mean_correlation")) {
      return(call(sn_plot_bulk_qc, x = result, metric = view))
    }
    if (identical(view, "pca")) return(call(sn_plot_bulk_pca, x = result))
    if (identical(view, "correlation_heatmap")) return(call(sn_plot_sample_correlation, x = result, view = "heatmap"))
    return(call(sn_plot_sample_correlation, x = result, view = "scatter"))
  }
  if (identical(analysis_type, "bulk_de")) return(call(sn_plot_bulk_de, x = result))
  if (identical(analysis_type, "bulk_network")) return(call(sn_plot_wgcna, x = result, type = view))
  if (identical(analysis_type, "bulk_survival")) return(call(sn_plot_survival, x = result, view = view))
  if (identical(analysis_type, "cell_communication")) {
    if (identical(view, "ligand_target")) return(call(sn_plot_ligand_target, x = result))
    if (identical(view, "comparison")) return(call(sn_plot_communication_comparison, x = result))
    return(call(sn_plot_communication, x = result, type = view))
  }
  if (identical(analysis_type, "cnv")) return(call(sn_plot_cnv, x = result, type = view))
  if (identical(analysis_type, "grn")) return(call(sn_plot_regulon, x = result, type = view))
  if (identical(analysis_type, "metabolism")) return(call(sn_plot_metabolism, x = result, type = view))
  if (identical(analysis_type, "program_discovery")) return(call(sn_plot_discovered_programs, x = result, type = view))
  if (identical(analysis_type, "program_scoring")) {
    if (identical(view, "heatmap")) return(call(sn_plot_program_heatmap, x = result, result_id = result_id %||% "default"))
    return(call(sn_plot_program_activity, x = result, result_id = result_id %||% "default"))
  }
  if (identical(analysis_type, "state_priority")) return(call(sn_plot_state_priority, x = result))
  if (identical(analysis_type, "scissor")) return(call(sn_plot_scissor, x = result, type = view))
  if (identical(analysis_type, "trajectory")) {
    if (identical(view, "embedding")) return(call(sn_plot_trajectory, x = result))
    if (identical(view, "pseudotime")) return(call(sn_plot_pseudotime, x = result))
    if (identical(view, "lineage_probability")) return(call(sn_plot_lineage_probability, x = result))
    if (identical(view, "dynamic_heatmap")) return(call(sn_plot_dynamic_heatmap, x = result))
    if (identical(view, "gene_trend")) return(call(sn_plot_gene_trend, x = result))
    return(call(sn_plot_branch_comparison, x = result))
  }
  if (identical(analysis_type, "velocity")) return(call(sn_plot_velocity, x = result))
  if (identical(analysis_type, "fate")) return(call(sn_plot_fate, x = result))
  if (identical(analysis_type, "spatial_domains")) return(call(sn_plot_spatial_domain, x = result))
  if (identical(analysis_type, "spatial_features")) return(call(sn_plot_spatial_svg, x = result))
  if (identical(analysis_type, "spatial_neighborhood")) return(call(sn_plot_spatial_neighborhood, x = result, type = view))
  if (identical(analysis_type, "spatial_communication")) return(call(sn_plot_spatial_communication, x = result))
  if (identical(analysis_type, "qc")) return(call(sn_plot_qc, x = result, metric = view))
  if (identical(analysis_type, "resolution_sweep")) return(call(sn_plot_resolution_sweep, x = result))
  if (identical(analysis_type, "integration_assessment")) return(call(sn_plot_integration, x = result))
  stop("No plot dispatcher is registered for analysis type '", analysis_type, "'.", call. = FALSE)
}

#' Plot a stored or direct Shennong analysis result
#'
#' `sn_plot_result()` is the canonical result-aware visualization entry point.
#' It uses one stable interface across analysis domains: `object`,
#' `analysis_type`, `result_id`, and `view`. Existing domain-specific
#' `sn_plot_*()` functions remain available as compatibility and advanced
#' interfaces.
#'
#' A direct Shennong result supplies its own `analysis_type`. For Seurat input,
#' the function retrieves a stored result. If `analysis_type` or `result_id` is
#' omitted, it is selected only when unambiguous; a stored result named
#' `"default"` is preferred within a selected analysis type.
#'
#' @param object A Seurat object with stored results, a direct Shennong result,
#'   or a compatible table/legacy assessment supported by the chosen analysis
#'   type.
#' @param analysis_type Optional analysis type. Required for ambiguous Seurat
#'   stores and for table/legacy inputs. See [sn_list_plot_methods()].
#' @param result_id Optional stored result identifier. A unique result or one
#'   named `"default"` is selected automatically when omitted.
#' @param view Optional registered plot view. The analysis-specific default is
#'   used when omitted.
#' @param ... View-specific parameters forwarded to the specialized plotter.
#'
#' @return Usually a `ggplot` object; heatmap or multi-panel backends may return
#'   a compatible composed plot object.
#'
#' @examples
#' de_table <- data.frame(
#'   gene = c("G1", "G2", "G3"),
#'   log2_fold_change = c(2, -1.5, 0.2),
#'   adjusted_p_value = c(0.01, 0.03, 0.8)
#' )
#' sn_plot_result(de_table, analysis_type = "de", view = "volcano")
#' \dontrun{
#' sn_list_results(seurat_object, type = "de")
#' sn_plot_result(
#'   seurat_object,
#'   analysis_type = "de",
#'   result_id = "cluster_markers",
#'   view = "volcano"
#' )
#' }
#' @export
sn_plot_result <- function(object,
                           analysis_type = NULL,
                           result_id = NULL,
                           view = NULL,
                           ...) {
  input <- .sn_resolve_canonical_plot_input(object, analysis_type, result_id)
  registry <- .sn_plot_method_registry()
  entry <- registry[[input$analysis_type]]
  if (is_null(entry)) {
    stop(
      "No plot methods are registered for analysis type '", input$analysis_type,
      "'. Use `sn_list_plot_methods()` to inspect supported types.",
      call. = FALSE
    )
  }
  if (is_null(view)) view <- entry$default
  if (!is.character(view) || length(view) != 1L || is.na(view) || !view %in% entry$views) {
    stop(
      "Unknown view for analysis type '", input$analysis_type, "'. Choose one of: ",
      paste(entry$views, collapse = ", "), ".",
      call. = FALSE
    )
  }
  .sn_dispatch_plot_result(input, view, list(...))
}
