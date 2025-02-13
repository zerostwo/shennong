#' @export
sn_add_data_from_anndata <- function(object, umap_path = NULL, metadata_path = NULL) {
  if (!is_null(x = umap_path)) {
    umap <- read.csv(file = umap_path, row.names = 1) |>
      Matrix::as.matrix()
    object <- object[, rownames(umap)]
    object[["umap"]] <- Seurat::CreateDimReducObject(umap, key = "umap")
  }
  if (!is_null(x = metadata_path)) {
    metadata <- read.csv(file = metadata_path, row.names = 1)
    object@meta.data <- metadata
  }
  object
}

#' @export
sn_write_h5ad <- function(object, path, mode = "w") {
  dir_path <- dirname(path = path)
  if (!dir.exists(paths = dir_path)) {
    dir.create(path = dir_path, recursive = TRUE)
  }

  if (!inherits(x = object, what = "SingleCellExperiment")) {
    object <- rlang::try_fetch(
      Seurat::as.SingleCellExperiment(x = object),
      error = function(cnd) {
        abort(
          "Failed to convert object to SingleCellExperiment.",
          parent = cnd
        )
      }
    )
  }

  rlang::try_fetch(
    anndataR::write_h5ad(object = object, path = path, mode = mode),
    error = function(cnd) {
      abort(
        paste("Failed to write h5ad file to:", path),
        parent = cnd
      )
    }
  )

  rlang::inform(paste("Successfully wrote h5ad file to:", path))
}
