.sn_is_seurat_object <- function(object) {
  methods::is(object, "Seurat")
}

.sn_extract_feature_names <- function(object) {
  if (.sn_is_seurat_object(object)) {
    assay <- SeuratObject::DefaultAssay(object)
    return(rownames(object[[assay]]))
  }

  if (inherits(object, c("matrix", "Matrix", "data.frame"))) {
    return(rownames(object))
  }

  if (is.character(object)) {
    return(object)
  }

  NULL
}

.sn_infer_species_from_features <- function(features) {
  features <- unique(stats::na.omit(as.character(features)))
  features <- features[nzchar(features)]
  if (length(features) == 0) {
    return(NULL)
  }

  homology_env <- new.env(parent = emptyenv())
  utils::data("hom_genes", package = "Shennong", envir = homology_env)
  homology_table <- get("hom_genes", envir = homology_env, inherits = FALSE)

  human_matches <- sum(features %in% homology_table$human)
  mouse_matches <- sum(features %in% homology_table$mouse)

  human_mt <- sum(grepl("^MT-", features))
  mouse_mt <- sum(grepl("^mt-", features))
  human_like <- human_matches + human_mt
  mouse_like <- mouse_matches + mouse_mt

  min_confident_matches <- max(3L, ceiling(length(features) * 0.01))
  if (max(human_like, mouse_like) < min_confident_matches) {
    return(NULL)
  }

  if (human_like == mouse_like) {
    return(NULL)
  }

  if (human_like > mouse_like) "human" else "mouse"
}

#' Retrieve or infer species information
#'
#' This helper returns the explicit \code{species} argument when supplied,
#' otherwise it tries the species stored in a Seurat object and then falls back
#' to feature-name based inference using the packaged \code{hom_genes} mapping
#' plus common mitochondrial naming patterns.
#'
#' @param object A Seurat object, matrix-like object, or character vector of
#'   gene symbols.
#' @param species Optional explicit species label. If provided, it is returned
#'   directly.
#'
#' @return A character string indicating the inferred or explicit species.
#' @export
sn_get_species <- function(object, species = NULL) {
  if (!is_null(species)) {
    return(arg_match(species, c("human", "mouse")))
  }

  if (.sn_is_seurat_object(object)) {
    species <- Seurat::Misc(object, slot = "species")
    if (!is_null(species)) {
      return(species)
    }
  }

  inferred_species <- .sn_infer_species_from_features(.sn_extract_feature_names(object))
  if (is_null(inferred_species)) {
    stop(
      paste(
        "Species information is required but could not be inferred from feature names.",
        "Please provide `species = \"human\"` or `species = \"mouse\"`."
      ),
      call. = FALSE
    )
  }

  if (.sn_is_seurat_object(object)) {
    Seurat::Misc(object = object, slot = "species") <- inferred_species
  }

  inferred_species
}
