.sn_supported_signature_categories <- function() {
  names(shennong_signatures$human)
}

#' Retrieve bundled Shennong signature genes by preset category names
#'
#' Shennong ships a package-owned snapshot of the signature categories required
#' by its core workflows. The snapshot is generated from \pkg{SignatuR} during
#' package development and stored inside the package so runtime behavior is
#' stable and does not depend on the user's installed \pkg{SignatuR} version.
#'
#' @param species One of \code{"human"} or \code{"mouse"}.
#' @param category A character vector of predefined signature categories.
#'
#' @return A unique character vector of signature gene symbols.
#'
#' @examples
#' block_genes <- sn_get_signatures(
#'   species = "human",
#'   category = c("mito", "ribo", "tcr", "immunoglobulins")
#' )
#' @export
sn_get_signatures <- function(species = "human",
                              category = NULL) {
  species <- rlang::arg_match(species, c("human", "mouse"))
  if (is_null(category) || length(category) == 0) {
    stop("`category` must contain at least one signature category.")
  }
  category <- unique(category)
  gene_list <- shennong_signatures[[species]]
  supported_categories <- names(gene_list)
  unknown_categories <- setdiff(category, supported_categories)

  if (length(unknown_categories) > 0) {
    warning(
      glue(
        "Unknown signature categories: {paste(unknown_categories, collapse = ', ')}. ",
        "Supported categories are: {paste(supported_categories, collapse = ', ')}."
      ),
      call. = FALSE
    )
  }

  selected_categories <- intersect(category, supported_categories)
  unique(unlist(gene_list[selected_categories], use.names = FALSE))
}
