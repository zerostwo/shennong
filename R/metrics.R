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
#' muf_calculate_lisi(x = seurat_obj)
#'
#' @export
sn_calculate_lisi <-
  function(x,
           reduction = "pca",
           label = "sample") {
    check_installed("lisi", action = \(pkg, ...) pak::pak("immunogenomics/lisi"))
    check_installed("SeuratObject")
    lisi_score <- lisi::compute_lisi(
      X = SeuratObject::Embeddings(object = x, reduction = reduction),
      meta_data = x@meta.data,
      label_colnames = label
    ) |> rownames_to_column("cell_id")
    return(lisi_score)
  }

#' @export
sn_calculate_rogue <-
  function(x,
           cluster = NULL,
           sample = NULL,
           span = 0.9) {
    check_installed("ROGUE", action = \(pkg, ...) pak::pak("PaulingLiu/ROGUE"))
    counts <- SeuratObject::LayerData(object = x, layer = "counts")
    counts <- Matrix::as.matrix(counts)
    metadata <- x@meta.data
    counts <-
      ROGUE::matr.filter(counts, min.cells = 10, min.genes = 10)
    entropy <- ROGUE::SE_fun(counts)
    rogue <- ROGUE::CalculateRogue(entropy, platform = "UMI")
    if (!is.null(cluster) && !is.null(sample)) {
      rogue <-
        ROGUE::rogue(
          counts,
          labels = as.character(metadata[[cluster]]),
          samples = as.character(metadata[[sample]]),
          platform = "UMI",
          span = span
        )
    }
    return(rogue)
  }
