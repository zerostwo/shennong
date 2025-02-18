#' Retrieve Species Information from a Seurat Object
#'
#' This function extracts the species information from a Seurat object. If the species is not found,
#' it throws an error.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param species (Optional) A character string indicating the species. If provided, it is returned directly.
#'
#' @return A character string indicating the species.
#' @export
sn_get_species <- function(object, species = NULL) {
  if (!is.null(species)) {
    return(species)
  }

  species <- Seurat::Misc(object, slot = "species")
  if (is.null(species)) {
    stop("Species information is required but not found in the Seurat object. Please provide a valid species name.")
  }

  return(species)
}
