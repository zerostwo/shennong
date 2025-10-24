sn_find_de <- function(
  object,
  treatment = NULL,
  control = NULL,
  group_by = NULL,
  assay = NULL,
  features = NULL,
  method = "wilcox",
  only_pos = NULL,
  logfc_threshold = 0.1,
  min_pct = 0.25,
  de_p_val = 0.05,
  de_logfc = 0.25,
  verbose = TRUE,
  ...
) {
  if (all(is_null(c(treatment, control, group_by)))) {
    stop("At least one of treatment, control, group_by must be provided.")
  }

  if (is_null(treatment) && is_null(control) && !is_null(group_by)) {
    only_pos <- only_pos %||% TRUE
    adata <- Seurat::FindAllMarkers(
      object = object, assay = assay, features = features,
      group.by = group_by, test.use = method, only_pos = TRUE,
      min.pct = min_pct, logfc.threshold = logfc_threshold, verbose = verbose
    )
    return(adata)
  }
  if (!(is_null(treatment) || is_null(control))) {
    only_pos <- only_pos %||% FALSE
    deg <- FindMarkers(
      object = seurat_obj,
      ident.1 = treatment,
      ident.2 = control,
      min.pct = min_pct,
      group.by = group_by,
      assay = assay,
      only.pos = only_pos,
      verbose = verbose,
      logfc.threshold = logfc_threshold
    )
  }
}
