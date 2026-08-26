# Evidence preparation and prompt building for LLM interpretation.
#
# Extracted from interpretation.R: annotation-response parsing and label normalization,
# cluster/marker/enrichment/geometry evidence summaries, lineage hints and canonical
# snapshots, public sn_prepare_*_evidence builders, and prompt rendering helpers.

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
  .sn_validate_seurat_object(object)
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


#' Prepare cluster-annotation evidence from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param de_name Optional stored marker-result name in
#'   \code{object@misc$de_results}. When omitted, Shennong prefers
#'   \code{"default"}, then a single available result, and otherwise the most
#'   recent marker result.
#' @param cluster_by Metadata column containing cluster labels.
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
#'     cluster_by = "cell_type"
#'   )
#'   names(evidence)
#' }
#' @export
sn_prepare_annotation_evidence <- function(object,
                                           de_name = NULL,
                                           cluster_by = NULL,
                                           n_markers = 10,
                                           marker_selection = c("specific", "top"),
                                           enrichment_name = NULL,
                                           n_terms = 5,
                                           enrichment_selection = c("specific", "top"),
                                           include_qc = TRUE,
                                           reduction = "umap",
                                           n_neighbor_clusters = 3) {
  .sn_validate_seurat_object(object)
  cluster_by <- cluster_by %||% "seurat_clusters"
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
  cluster_summary <- .sn_prepare_cluster_summary(object = object, cluster_col = cluster_by)
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
  prediction_summary <- .sn_prepare_prediction_summary(object = object, cluster_col = cluster_by)
  enrichment_summary <- .sn_prepare_cluster_enrichment_summary(
    object = object,
    enrichment_name = enrichment_name,
    n_terms = n_terms,
    selection = enrichment_selection
  )
  qc_summary <- if (isTRUE(include_qc)) {
    .sn_prepare_annotation_qc_summary(object = object, cluster_col = cluster_by)
  } else {
    tibble::tibble()
  }
  lineage_hints <- .sn_prepare_annotation_lineage_hints(
    object = object,
    cluster_col = cluster_by,
    top_marker_table = marker_table,
    species = tryCatch(sn_get_species(object), error = function(...) NULL)
  )
  canonical_marker_snapshot <- .sn_prepare_annotation_canonical_snapshot(
    object = object,
    cluster_col = cluster_by,
    species = tryCatch(sn_get_species(object), error = function(...) NULL)
  )
  geometry_summary <- .sn_prepare_cluster_geometry_summary(
    object = object,
    cluster_col = cluster_by,
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
    cluster_col = cluster_by,
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
  .sn_validate_seurat_object(object)

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
#' @param cluster_by Metadata column containing cluster labels.
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
#'     cluster_by = "cell_type"
#'   )
#'   names(evidence)
#' }
#' @export
sn_prepare_results_evidence <- function(object,
                                        cluster_de_name = NULL,
                                        contrast_de_name = NULL,
                                        enrichment_name = NULL,
                                        cluster_by = NULL,
                                        n_markers = 5,
                                        n_terms = 10) {
  .sn_validate_seurat_object(object)
  cluster_by <- cluster_by %||% "seurat_clusters"

  evidence <- list(
    task = "results",
    dataset = list(
      n_cells = ncol(object),
      n_features = nrow(object),
      cluster_col = cluster_by,
      clusters = if (cluster_by %in% colnames(object[[]])) nlevels(factor(object[[cluster_by]][, 1])) else NULL
    ),
    cluster_summary = .sn_prepare_cluster_summary(object = object, cluster_col = cluster_by)
  )

  if (!is_null(cluster_de_name)) {
    evidence$cluster_markers <- sn_prepare_annotation_evidence(
      object = object,
      de_name = cluster_de_name,
      cluster_by = cluster_by,
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
