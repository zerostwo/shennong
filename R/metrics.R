#' Calculate LISI score
#'
#' This function calculates the Local Intrinsic Dimensionality-based Outlier score (Lisi score) for a Seurat object.
#'
#' @param object A Seurat object.
#' @param reduction The dimensionality reduction method used to generate the embeddings (default is "pca").
#' @param label The column name in object@meta.data that specifies the sample labels (default is "sample").
#'
#' @return A data frame with the Lisi score for each cell, along with the cell ID.
#'
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' # Calculate Lisi score for a Seurat object
#' sn_calculate_lisi(x = seurat_obj)
#'
#' @export
sn_calculate_lisi <-
  function(x,
           reduction = "pca",
           label = "sample") {
    rlang::check_installed("lisi", action = \(pkg, ...) pak::pak("immunogenomics/lisi"))
    rlang::check_installed("SeuratObject")
    lisi_score <- lisi::compute_lisi(
      X = SeuratObject::Embeddings(object = x, reduction = reduction),
      meta_data = x@meta.data,
      label_colnames = label
    ) |> rownames_to_column("cell_id")
    return(lisi_score)
  }

#' Calculate ROGUE score for Seurat Object
#'
#' This function calculates ROGUE score based on Seurat object.
#'
#' @param x A Seurat object.
#' @param cluster Column name in metadata specifying cluster labels.
#' @param sample Column name in metadata specifying sample labels.
#' @param span The span parameter for rogue estimation.
#'
#' @return A data.frame of ROGUE score per cluster/sample or per cluster.
#' @export
sn_calculate_rogue <- function(x,
  cluster = NULL,
  sample = NULL,
  span = 0.9) {
rlang::check_installed("ROGUE", action = \(pkg, ...) pak::pak("PaulingLiu/ROGUE"))

if (!inherits(x, "Seurat")) {
stop("Input x must be a Seurat object.")
}

if (!"counts" %in% SeuratObject::LayerNames(x)) {
stop("The 'counts' layer is not found in the Seurat object.")
}

counts <- SeuratObject::LayerData(x, layer = "counts")
counts <- Matrix::as.matrix(counts)
metadata <- x@meta.data

counts <- ROGUE::matr.filter(counts, min.cells = 10, min.genes = 10)

message("Calculating entropy...")
entropy <- ROGUE::SE_fun(counts)

message("Calculating ROGUE score...")
rogue_result <- ROGUE::CalculateRogue(entropy, platform = "UMI")

if (!is.null(cluster) && !is.null(sample)) {
if (!cluster %in% colnames(metadata)) {
stop("Specified cluster column not found in metadata.")
}
if (!sample %in% colnames(metadata)) {
stop("Specified sample column not found in metadata.")
}

rogue_result <- ROGUE::rogue(
counts,
labels = as.character(metadata[[cluster]]),
samples = as.character(metadata[[sample]]),
platform = "UMI",
span = span
)

rogue_result <- as.data.frame(rogue_result)
rogue_result$cluster <- rownames(rogue_result)
rownames(rogue_result) <- NULL
}

return(rogue_result)
}