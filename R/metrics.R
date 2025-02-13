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
