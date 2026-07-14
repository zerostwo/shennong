.sn_resolve_cnv_result <- function(x, name = NULL) {
  if (inherits(x, "Seurat")) {
    if (is_null(name) || !nzchar(name)) stop("`name` is required when `x` is a Seurat object.", call. = FALSE)
    result <- sn_get_result(x, "cnv", name)
  } else {
    result <- x
    sn_validate_result(result)
  }
  if (!identical(result$analysis_type, "cnv")) stop("Expected a CNV result.", call. = FALSE)
  result
}

#' Plot CNV, malignancy, and subclone results
#'
#' @param x A Seurat object or unified CNV result.
#' @param name Stored result name when `x` is a Seurat object.
#' @param type Plot type: chromosome heatmap, CNV UMAP, score distribution,
#'   sample summary, or CNV-expression association.
#' @param n Maximum cells or features shown.
#' @return A `ggplot` object.
#' @examples
#' \dontrun{sn_plot_cnv(object, "cnv", type = "heatmap")}
#' @export
sn_plot_cnv <- function(x,
                        name = NULL,
                        type = c("heatmap", "umap", "score", "sample", "association"),
                        n = 100L) {
  type <- match.arg(type)
  result <- .sn_resolve_cnv_result(x, name)
  primary <- tibble::as_tibble(result$tables$primary)
  if (identical(type, "score")) {
    return(ggplot2::ggplot(primary, ggplot2::aes(x = .data$malignant_call, y = .data$malignant_score, fill = .data$malignant_call)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.15, size = 0.6, alpha = 0.5) +
      ggplot2::labs(x = NULL, y = "Reference-scaled malignancy score", fill = NULL) +
      ggplot2::theme_bw())
  }
  if (identical(type, "sample")) {
    data <- tibble::as_tibble(result$tables$sample_summary)
    if (nrow(data) == 0L) stop("CNV sample summary is empty.", call. = FALSE)
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$sample, y = .data$malignant_fraction, fill = .data$mean_cnv_score)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(x = NULL, y = "Malignant fraction", fill = "Mean CNV score") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
  }
  if (identical(type, "association")) {
    data <- tibble::as_tibble(result$tables$expression_association)
    data <- data[is.finite(data$correlation), , drop = FALSE]
    data <- utils::head(data[order(abs(data$correlation), decreasing = TRUE), , drop = FALSE], as.integer(n))
    if (nrow(data) == 0L) stop("CNV-expression association table is empty.", call. = FALSE)
    data$feature <- factor(data$feature, levels = rev(data$feature))
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$correlation, y = .data$feature, colour = .data$adjusted_p_value)) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$correlation, yend = .data$feature), colour = "grey70") +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_viridis_c(direction = -1, trans = "reverse") +
      ggplot2::labs(x = "Spearman correlation with malignancy score", y = NULL, colour = "BH adjusted p") +
      ggplot2::theme_bw())
  }
  if (identical(type, "umap")) {
    embedding <- result$embeddings$cnv_umap
    if (is_null(embedding)) stop("CNV UMAP is unavailable; rerun the backend with UMAP enabled.", call. = FALSE)
    embedding <- as.data.frame(embedding)
    if (ncol(embedding) < 2L) stop("CNV UMAP must have at least two dimensions.", call. = FALSE)
    data <- tibble::tibble(
      cell = rownames(embedding), x = embedding[[1]], y = embedding[[2]]
    ) |>
      dplyr::left_join(primary[, c("cell", "malignant_score", "subclone")], by = "cell")
    return(ggplot2::ggplot(data, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$malignant_score)) +
      ggplot2::geom_point(size = 0.8, alpha = 0.85) +
      ggplot2::scale_colour_viridis_c() +
      ggplot2::labs(x = "CNV UMAP 1", y = "CNV UMAP 2", colour = "Malignancy") +
      ggplot2::theme_bw())
  }
  chromosomes <- tibble::as_tibble(result$tables$chromosome)
  if (nrow(chromosomes) == 0L) stop("Chromosome-level CNV table is empty.", call. = FALSE)
  selected <- utils::head(primary[order(primary$malignant_score, decreasing = TRUE), "cell", drop = TRUE], as.integer(n))
  chromosomes <- chromosomes[chromosomes$cell %in% selected & is.finite(chromosomes$cnv), , drop = FALSE]
  chromosomes$cell <- factor(chromosomes$cell, levels = rev(selected))
  ggplot2::ggplot(chromosomes, ggplot2::aes(x = .data$chromosome, y = .data$cell, fill = .data$cnv)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0) +
    ggplot2::labs(x = "Chromosome", y = NULL, fill = "CNV") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
