#' Run multimodal clustering through the unified clustering workflow
#'
#' `sn_run_multimodal()` is the explicit multimodal entry point for the
#' CITE-seq workflows implemented by [sn_run_cluster()]. It forwards the
#' requested WNN, totalVI, Coralysis, or MMoCHi method without changing the
#' clustering return contract.
#'
#' @inheritParams sn_run_cluster
#' @param modality Multimodal assay type. CITE-seq is currently supported.
#' @param method Multimodal method: `"wnn"`, `"totalvi"`, `"coralysis"`, or
#'   `"mmochi"`.
#' @param ... Additional arguments passed to [sn_run_cluster()].
#'
#' @return The clustered Seurat object, or the cluster vector when
#'   `return_cluster = TRUE`.
#' @export
sn_run_multimodal <- function(object,
                              modality = "cite_seq",
                              method = c("wnn", "totalvi", "coralysis", "mmochi"),
                              ...) {
  modality <- match.arg(modality, "cite_seq")
  method <- match.arg(method)
  sn_run_cluster(
    object = object,
    modality = modality,
    multimodal_method = method,
    ...
  )
}
