.sn_validate_annotation_evidence <- function(evidence) {
  if (!is.data.frame(evidence)) {
    stop("`evidence` must be a data frame.", call. = FALSE)
  }
  required <- c("entity", "label", "score", "method")
  missing <- setdiff(required, colnames(evidence))
  if (length(missing) > 0L) {
    stop("`evidence` is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  evidence <- tibble::as_tibble(evidence)
  evidence$entity <- as.character(evidence$entity)
  evidence$label <- as.character(evidence$label)
  evidence$method <- as.character(evidence$method)
  evidence$score <- suppressWarnings(as.numeric(evidence$score))
  evidence <- evidence[
    !is.na(evidence$entity) & nzchar(evidence$entity) &
      !is.na(evidence$label) & nzchar(evidence$label) &
      is.finite(evidence$score),
    ,
    drop = FALSE
  ]
  if (nrow(evidence) == 0L) {
    stop("`evidence` contains no complete finite scores.", call. = FALSE)
  }
  evidence
}

.sn_normalize_annotation_scores <- function(evidence) {
  groups <- split(seq_len(nrow(evidence)), paste(evidence$entity, evidence$method, sep = "\r"))
  evidence$normalized_score <- NA_real_
  for (indices in groups) {
    values <- evidence$score[indices]
    value_range <- range(values, finite = TRUE)
    if (length(values) == 1L) {
      normalized <- min(1, max(0, values))
    } else if (diff(value_range) <= sqrt(.Machine$double.eps)) {
      normalized <- rep(0.5, length(values))
    } else {
      normalized <- (values - value_range[[1]]) / diff(value_range)
    }
    evidence$normalized_score[indices] <- normalized
  }
  evidence
}

#' Calibrate confidence from annotation evidence
#'
#' Scores are normalized within each entity/backend before weighted consensus.
#' The output retains the best and second-best labels, their margin, method
#' support, and a conservative low-confidence flag.
#'
#' @param evidence Long data frame with \code{entity}, \code{label},
#'   \code{score}, and \code{method} columns.
#' @param weights Optional named numeric vector of backend weights.
#' @param low_confidence_threshold Minimum calibrated prediction score.
#' @param margin_threshold Minimum best-versus-second-best margin.
#'
#' @return A tibble with one calibrated prediction per entity.
#'
#' @examples
#' evidence <- data.frame(
#'   entity = rep(c("0", "1"), each = 2),
#'   label = rep(c("B cells", "T cells"), 2),
#'   score = c(2, 0.2, 0.1, 3), method = "markers"
#' )
#' sn_annotation_confidence(evidence)
#'
#' @export
sn_annotation_confidence <- function(evidence,
                                     weights = NULL,
                                     low_confidence_threshold = 0.55,
                                     margin_threshold = 0.1) {
  evidence <- .sn_normalize_annotation_scores(.sn_validate_annotation_evidence(evidence))
  methods <- unique(evidence$method)
  if (is_null(weights)) {
    weights <- stats::setNames(rep(1, length(methods)), methods)
  }
  if (is.null(names(weights)) || any(!methods %in% names(weights))) {
    stop("`weights` must be named for every method present in `evidence`.", call. = FALSE)
  }
  weights <- as.numeric(weights[methods]) |>
    stats::setNames(methods)
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("All annotation `weights` must be positive finite values.", call. = FALSE)
  }

  entity_groups <- split(evidence, evidence$entity)
  out <- lapply(names(entity_groups), function(entity) {
    current <- entity_groups[[entity]]
    available_methods <- unique(current$method)
    denominator <- sum(weights[available_methods])
    labels <- unique(current$label)
    label_scores <- vapply(labels, function(label) {
      rows <- current[current$label == label, , drop = FALSE]
      sum(rows$normalized_score * weights[rows$method]) / denominator
    }, numeric(1))
    ordering <- order(label_scores, decreasing = TRUE, na.last = TRUE)
    best_index <- ordering[[1]]
    second_index <- if (length(ordering) >= 2L) ordering[[2]] else NA_integer_
    best_score <- unname(label_scores[[best_index]])
    second_score <- if (is.na(second_index)) 0 else unname(label_scores[[second_index]])
    margin <- best_score - second_score
    calibrated <- min(1, max(0, 0.7 * best_score + 0.3 * margin))
    tibble::tibble(
      entity = entity,
      prediction = labels[[best_index]],
      second_best_label = if (is.na(second_index)) NA_character_ else labels[[second_index]],
      prediction_score = calibrated,
      margin = margin,
      method_count = length(available_methods),
      methods = paste(sort(available_methods), collapse = ";"),
      low_confidence = calibrated < low_confidence_threshold || margin < margin_threshold
    )
  })
  dplyr::bind_rows(out)
}

.sn_annotation_collapse_field <- function(x) {
  x <- unique(trimws(unlist(strsplit(as.character(x[!is.na(x)]), ";", fixed = TRUE))))
  x <- x[nzchar(x)]
  paste(x, collapse = ";")
}

#' Build consensus annotation labels from evidence
#'
#' @param evidence Long evidence table accepted by
#'   \code{sn_annotation_confidence()}. Optional columns include
#'   \code{parent_label}, \code{supporting_markers},
#'   \code{conflicting_markers}, and \code{reference_coverage}.
#' @param weights Optional named backend weights.
#' @param ontology Logical; map final labels to the bundled Cell Ontology
#'   snapshot.
#' @param ontology_mapping Optional custom ontology mapping.
#' @param low_confidence_threshold,margin_threshold Confidence thresholds.
#'
#' @return A traceable annotation tibble with hierarchical labels, evidence,
#'   confidence, and ontology identifiers.
#'
#' @examples
#' evidence <- data.frame(
#'   entity = rep(c("0", "1"), each = 2),
#'   label = rep(c("B cells", "T cells"), 2),
#'   score = c(3, 0.1, 0.2, 2.5), method = "markers",
#'   parent_label = rep(c("B cells", "T cells"), 2)
#' )
#' sn_annotation_consensus(evidence)
#'
#' @export
sn_annotation_consensus <- function(evidence,
                                    weights = NULL,
                                    ontology = TRUE,
                                    ontology_mapping = NULL,
                                    low_confidence_threshold = 0.55,
                                    margin_threshold = 0.1) {
  evidence <- .sn_validate_annotation_evidence(evidence)
  confidence <- sn_annotation_confidence(
    evidence,
    weights = weights,
    low_confidence_threshold = low_confidence_threshold,
    margin_threshold = margin_threshold
  )
  selected <- lapply(seq_len(nrow(confidence)), function(index) {
    entity <- confidence$entity[[index]]
    label <- confidence$prediction[[index]]
    rows <- evidence[evidence$entity == entity & evidence$label == label, , drop = FALSE]
    field <- function(name, default = NA_character_) {
      if (!name %in% colnames(rows)) return(default)
      value <- .sn_annotation_collapse_field(rows[[name]])
      if (nzchar(value)) value else default
    }
    coverage <- if ("reference_coverage" %in% colnames(rows)) {
      values <- suppressWarnings(as.numeric(rows$reference_coverage))
      if (all(is.na(values))) NA_real_ else max(values, na.rm = TRUE)
    } else {
      NA_real_
    }
    tibble::tibble(
      entity = entity,
      level_1 = field("parent_label", NA_character_),
      level_2 = label,
      level_3 = label,
      supporting_markers = field("supporting_markers"),
      conflicting_markers = field("conflicting_markers"),
      reference_coverage = coverage
    )
  }) |>
    dplyr::bind_rows()
  out <- dplyr::left_join(confidence, selected, by = "entity")
  if (isTRUE(ontology)) {
    mapped <- sn_map_cell_ontology(out$prediction, ontology = ontology_mapping)
    out$ontology_id <- mapped$ontology_id
    out$ontology_label <- mapped$ontology_label
  } else {
    out$ontology_id <- NA_character_
    out$ontology_label <- NA_character_
  }
  out
}
