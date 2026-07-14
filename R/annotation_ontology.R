.sn_cell_ontology_path <- function() {
  installed <- system.file("ontology", "cell_ontology.json", package = "Shennong")
  candidates <- c(
    if (nzchar(installed)) installed else character(),
    file.path(getwd(), "inst", "ontology", "cell_ontology.json")
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) {
    stop("The bundled Cell Ontology snapshot was not found.", call. = FALSE)
  }
  candidates[[1]]
}

.sn_cell_ontology_terms <- function(ontology = NULL) {
  if (is_null(ontology)) {
    payload <- jsonlite::fromJSON(.sn_cell_ontology_path(), simplifyVector = FALSE)
    ontology <- payload$terms
  } else if (is.character(ontology) && length(ontology) == 1L) {
    payload <- jsonlite::fromJSON(ontology, simplifyVector = FALSE)
    ontology <- payload$terms %||% payload
  }

  if (is.data.frame(ontology)) {
    required <- c("id", "label")
    missing <- setdiff(required, colnames(ontology))
    if (length(missing) > 0L) {
      stop("`ontology` is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    return(lapply(seq_len(nrow(ontology)), function(index) {
      list(
        id = as.character(ontology$id[[index]]),
        label = as.character(ontology$label[[index]]),
        aliases = if ("aliases" %in% colnames(ontology)) ontology$aliases[[index]] else character()
      )
    }))
  }
  if (!is.list(ontology)) {
    stop("`ontology` must be NULL, a JSON path, a data frame, or a list of terms.", call. = FALSE)
  }
  ontology
}

.sn_normalize_cell_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[_/-]+", " ", x)
  x <- gsub("[^[:alnum:]+ ]", "", x)
  gsub("[[:space:]]+", " ", x)
}

#' Map cell labels to Cell Ontology identifiers
#'
#' Uses a small, versioned Cell Ontology snapshot shipped with Shennong. A
#' project-specific mapping can be supplied as a JSON file, data frame, or list
#' with \code{id}, \code{label}, and optional \code{aliases} fields.
#'
#' @param labels Character vector of cell-type labels.
#' @param ontology Optional custom ontology mapping.
#' @param strict If \code{TRUE}, fail when any label is unmapped.
#'
#' @return A tibble with input label, ontology identifier, canonical ontology
#'   label, and the alias that matched.
#'
#' @examples
#' sn_map_cell_ontology(c("B cells", "T cells", "unknown"))
#'
#' @export
sn_map_cell_ontology <- function(labels, ontology = NULL, strict = FALSE) {
  if (!is.character(labels)) {
    stop("`labels` must be a character vector.", call. = FALSE)
  }
  terms <- .sn_cell_ontology_terms(ontology)
  lookup <- list()
  for (term in terms) {
    if (is_null(term$id) || is_null(term$label)) next
    aliases <- unique(c(term$label, unlist(term$aliases %||% character(), use.names = FALSE)))
    for (alias in aliases) {
      key <- .sn_normalize_cell_label(alias)
      if (!nzchar(key) || !is_null(lookup[[key]])) next
      lookup[[key]] <- list(
        ontology_id = as.character(term$id),
        ontology_label = as.character(term$label),
        matched_alias = as.character(alias)
      )
    }
  }

  rows <- lapply(labels, function(label) {
    match <- lookup[[.sn_normalize_cell_label(label)]]
    tibble::tibble(
      input_label = label,
      ontology_id = match$ontology_id %||% NA_character_,
      ontology_label = match$ontology_label %||% NA_character_,
      matched_alias = match$matched_alias %||% NA_character_
    )
  })
  out <- dplyr::bind_rows(rows)
  if (isTRUE(strict) && anyNA(out$ontology_id)) {
    stop(
      "No Cell Ontology mapping was found for: ",
      paste(unique(out$input_label[is.na(out$ontology_id)]), collapse = ", "),
      call. = FALSE
    )
  }
  out
}
