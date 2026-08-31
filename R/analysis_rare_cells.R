# Rare-cell detection (GapClust / ScCAD / gini / local-marker backends).
#
# Extracted from analysis_clustering.R: the detection subsystem behind
# sn_detect_rare_cells().

.sn_gini_coefficient <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x >= 0]
  if (length(x) == 0 || sum(x) <= 0) {
    return(0)
  }
  x <- sort(x, decreasing = FALSE)
  n <- length(x)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - ((n + 1) / n)
}

.sn_detect_rare_features_gini <- function(expr,
                                          nfeatures = 200,
                                          min_cells = 3,
                                          max_fraction = 0.1) {
  stopifnot(nrow(expr) > 0, ncol(expr) > 1)

  prevalence <- Matrix::rowSums(expr > 0)
  prevalence_fraction <- prevalence / ncol(expr)
  keep <- prevalence >= min_cells & prevalence_fraction <= max_fraction

  if (!any(keep)) {
    return(data.frame(
      feature = character(0),
      score = numeric(0),
      prevalence = integer(0),
      prevalence_fraction = numeric(0),
      mean_expression = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  kept_features <- rownames(expr)[keep]
  gini_score <- vapply(kept_features, function(feature) {
    .sn_gini_coefficient(expr[feature, ])
  }, numeric(1))
  mean_expression <- Matrix::rowMeans(expr[kept_features, , drop = FALSE])

  out <- data.frame(
    feature = kept_features,
    score = gini_score,
    prevalence = as.integer(prevalence[keep]),
    prevalence_fraction = as.numeric(prevalence_fraction[keep]),
    mean_expression = as.numeric(mean_expression),
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$score, out$prevalence_fraction, -out$mean_expression, out$feature), , drop = FALSE]
  utils::head(out, nfeatures)
}

.sn_make_temporary_grouping <- function(object,
                                        features,
                                        npcs = 20,
                                        dims = 1:10,
                                        resolution = 0.2) {
  temp_object <- .sn_with_default_seurat_acceleration(
    suppressWarnings(
      Seurat::ScaleData(
        object = object,
        features = features,
        verbose = FALSE
      )
    ),
    object = object
  )
  temp_object <- .sn_with_default_seurat_acceleration(
    suppressWarnings(
      Seurat::RunPCA(
        object = temp_object,
        features = features,
        npcs = npcs,
        verbose = FALSE,
        seed.use = 717
      )
    ),
    object = temp_object
  )
  temp_dims <- dims[dims <= npcs]
  temp_object <- Seurat::FindNeighbors(
    object = temp_object,
    reduction = "pca",
    dims = temp_dims,
    verbose = FALSE
  )
  temp_object <- Seurat::FindClusters(
    object = temp_object,
    resolution = resolution,
    random.seed = 717,
    verbose = FALSE
  )
  as.character(temp_object$seurat_clusters)
}

.sn_resolve_rare_groups <- function(object,
                                    group_by = NULL,
                                    features = NULL,
                                    npcs = 20,
                                    dims = 1:10,
                                    resolution = 0.2,
                                    max_fraction = 0.05,
                                    max_cells = 100) {
  groups <- if (!is.null(group_by)) {
    if (!group_by %in% colnames(object[[]])) {
      stop(glue("`rare_feature_group_by` column '{group_by}' was not found."), call. = FALSE)
    }
    as.character(object[[group_by, drop = TRUE]])
  } else {
    .sn_make_temporary_grouping(
      object = object,
      features = features,
      npcs = npcs,
      dims = dims,
      resolution = resolution
    )
  }

  group_sizes <- table(groups)
  rare_groups <- names(group_sizes)[
    group_sizes <= max_cells | (group_sizes / length(groups)) <= max_fraction
  ]

  list(
    groups = groups,
    rare_groups = rare_groups,
    group_sizes = group_sizes
  )
}

.sn_detect_rare_features_local_markers <- function(object,
                                                   groups,
                                                   rare_groups,
                                                   assay = "RNA",
                                                   layer = "data",
                                                   nfeatures = 200,
                                                   min_cells = 3) {
  if (length(rare_groups) == 0) {
    return(character(0))
  }

  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  local_markers <- lapply(rare_groups, function(current_group) {
    in_group <- groups == current_group
    if (sum(in_group) < min_cells || sum(!in_group) < min_cells) {
      return(character(0))
    }
    avg_in <- Matrix::rowMeans(expr[, in_group, drop = FALSE])
    avg_out <- Matrix::rowMeans(expr[, !in_group, drop = FALSE])
    pct_in <- Matrix::rowSums(expr[, in_group, drop = FALSE] > 0) / sum(in_group)
    pct_out <- Matrix::rowSums(expr[, !in_group, drop = FALSE] > 0) / sum(!in_group)
    score <- avg_in - avg_out
    keep <- is.finite(score) & pct_in > pct_out & pct_in > 0
    ranked <- names(sort(score[keep], decreasing = TRUE))
    utils::head(ranked, nfeatures)
  })

  utils::head(unique(unlist(local_markers, use.names = FALSE)), nfeatures)
}

.sn_resolve_rare_feature_control <- function(control = list()) {
  if (is.null(control)) {
    control <- list()
  }
  if (!is.list(control)) {
    stop("`rare_feature_control` must be a named list.", call. = FALSE)
  }

  defaults <- list(
    group_max_fraction = 0.05,
    group_max_cells = 100,
    gene_max_fraction = 0.1,
    min_cells = 3
  )
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0L) {
    stop(
      glue("Unknown `rare_feature_control` field(s): {paste(unknown, collapse = ', ')}."),
      call. = FALSE
    )
  }
  resolved <- utils::modifyList(defaults, control)

  resolved$group_max_fraction <- as.numeric(resolved$group_max_fraction)
  resolved$group_max_cells <- as.integer(resolved$group_max_cells)
  resolved$gene_max_fraction <- as.numeric(resolved$gene_max_fraction)
  resolved$min_cells <- as.integer(resolved$min_cells)

  if (!is.finite(resolved$group_max_fraction) || resolved$group_max_fraction <= 0 || resolved$group_max_fraction > 1) {
    stop("`rare_feature_control$group_max_fraction` must be in (0, 1].", call. = FALSE)
  }
  if (!is.finite(resolved$gene_max_fraction) || resolved$gene_max_fraction <= 0 || resolved$gene_max_fraction > 1) {
    stop("`rare_feature_control$gene_max_fraction` must be in (0, 1].", call. = FALSE)
  }
  if (is.na(resolved$group_max_cells) || resolved$group_max_cells < 1L) {
    stop("`rare_feature_control$group_max_cells` must be a positive integer.", call. = FALSE)
  }
  if (is.na(resolved$min_cells) || resolved$min_cells < 1L) {
    stop("`rare_feature_control$min_cells` must be a positive integer.", call. = FALSE)
  }

  resolved
}

.sn_select_rare_features <- function(object,
                                     base_features,
                                     method = "none",
                                     assay = "RNA",
                                     layer = "data",
                                     nfeatures = 200,
                                     group_by = NULL,
                                     control = list(),
                                     min_cells = 3,
                                     npcs = 20,
                                     dims = 1:10,
                                     resolution = 0.2,
                                     verbose = TRUE) {
  control <- .sn_resolve_rare_feature_control(control = control)
  method <- unique(match.arg(method, c("none", "gini", "local_markers"), several.ok = TRUE))
  method <- setdiff(method, "none")
  if (length(method) == 0) {
    return(list(
      features = character(0),
      metadata = data.frame(),
      groups = NULL
    ))
  }

  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  rare_features <- list()
  rare_meta <- list()
  rare_group_info <- NULL

  if ("local_markers" %in% method) {
    rare_group_info <- .sn_resolve_rare_groups(
      object = object,
      group_by = group_by,
      features = base_features,
      npcs = npcs,
      dims = dims,
      resolution = resolution,
      max_fraction = control$group_max_fraction,
      max_cells = control$group_max_cells
    )
  }

  if ("gini" %in% method) {
    gini_tbl <- .sn_detect_rare_features_gini(
      expr = expr,
      nfeatures = nfeatures,
      min_cells = control$min_cells %||% min_cells,
      max_fraction = control$gene_max_fraction
    )
    rare_features$gini <- gini_tbl$feature
    if (nrow(gini_tbl) > 0) {
      rare_meta[[length(rare_meta) + 1]] <- transform(gini_tbl, method = "gini")
    }
  }

  if ("local_markers" %in% method) {
    local_markers <- .sn_detect_rare_features_local_markers(
      object = object,
      groups = rare_group_info$groups,
      rare_groups = rare_group_info$rare_groups,
      assay = assay,
      layer = layer,
      nfeatures = nfeatures,
      min_cells = control$min_cells %||% min_cells
    )
    rare_features$local_markers <- local_markers
    if (length(local_markers) > 0) {
      rare_meta[[length(rare_meta) + 1]] <- data.frame(
        feature = local_markers,
        score = NA_real_,
        prevalence = NA_integer_,
        prevalence_fraction = NA_real_,
        mean_expression = NA_real_,
        method = "local_markers",
        stringsAsFactors = FALSE
      )
    }
  }

  combined <- unique(unlist(rare_features, use.names = FALSE))
  if (verbose) {
    .sn_log_info(
      "[rare_features] Added {length(combined)} rare-aware feature(s) from method(s): {paste(method, collapse = ', ')}"
    )
  }

  list(
    features = combined,
    metadata = .sn_bind_rows(rare_meta),
    groups = rare_group_info
  )
}

.sn_run_sccad <- function(expr,
                          cell_ids,
                          gene_ids,
                          python = NULL,
                          script = NULL,
                          normalization = FALSE,
                          seed = 2023,
                          rare_h = 0.01,
                          merge_h = 0.3,
                          overlap_h = 0.7,
                          save_full = FALSE) {
  python <- python %||% getOption("shennong.sccad_python", Sys.which("python"))
  if (!nzchar(python)) {
    stop("Could not find a Python executable for the scCAD backend.", call. = FALSE)
  }

  script <- script %||% getOption("shennong.sccad_script", Sys.getenv("SHENNONG_SCCAD_SCRIPT", unset = ""))
  if (!nzchar(script) || !file.exists(path.expand(script))) {
    stop(
      "Could not locate `scCAD.py`. Supply `sccad_script` or set `options(shennong.sccad_script = ...)`.",
      call. = FALSE
    )
  }
  script <- normalizePath(path.expand(script), winslash = "/", mustWork = TRUE)

  workdir <- tempfile("sccad_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  expr_path <- file.path(workdir, "expression.csv")
  output_path <- file.path(workdir, "sccad_result.json")
  runner_path <- file.path(workdir, "run_sccad.py")

  expr_df <- as.data.frame(t(as.matrix(expr)))
  colnames(expr_df) <- gene_ids
  rownames(expr_df) <- cell_ids
  utils::write.csv(expr_df, file = expr_path, quote = FALSE)

  runner_lines <- c(
    "import json",
    "import os",
    "import sys",
    "import numpy as np",
    "import pandas as pd",
    sprintf("sys.path.insert(0, %s)", shQuote(dirname(script))),
    sprintf("import %s as sccad_module", tools::file_path_sans_ext(basename(script))),
    sprintf("expr = pd.read_csv(%s, index_col=0)", shQuote(expr_path)),
    "data = expr.to_numpy(dtype=float)",
    "cell_names = expr.index.to_numpy()",
    "gene_names = expr.columns.to_numpy()",
    sprintf(
      paste0(
        "result, score, sub_clusters, degs_list = sccad_module.scCAD(",
        "data=data, dataName='Shennong', cellNames=cell_names, geneNames=gene_names, ",
        "normalization=%s, seed=%s, rare_h=%s, merge_h=%s, overlap_h=%s, save_full=%s, save_path=%s)"
      ),
      if (isTRUE(normalization)) "True" else "False",
      as.integer(seed),
      as.numeric(rare_h),
      as.numeric(merge_h),
      as.numeric(overlap_h),
      if (isTRUE(save_full)) "True" else "False",
      shQuote(workdir)
    ),
    "def _to_str_list(x):",
    "    out = []",
    "    for item in x:",
    "        if isinstance(item, bytes):",
    "            out.append(item.decode('utf-8'))",
    "        else:",
    "            out.append(str(item))",
    "    return out",
    "rare_sets = [_to_str_list(cluster) for cluster_by in result]",
    "sub_clusters = [str(x) for x in sub_clusters]",
    "score = [float(x) for x in score]",
    "payload = {",
    "    'rare_sets': rare_sets,",
    "    'scores': score,",
    "    'sub_clusters': sub_clusters,",
    "    'degs_list': [[str(g) for g in genes] for genes in degs_list]",
    "}",
    sprintf("with open(%s, 'w') as handle:", shQuote(output_path)),
    "    json.dump(payload, handle)"
  )
  writeLines(runner_lines, con = runner_path, useBytes = TRUE)

  status <- tryCatch(
    suppressWarnings(system2(
      command = python,
      args = c(shQuote(runner_path)),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) {
      stop(
        "SCA execution failed. Ensure the Python executable and the `shannonca` package are available. ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  exit_code <- attr(status, "status") %||% 0L
  if (!identical(exit_code, 0L) || !file.exists(output_path)) {
    stop(
      "scCAD execution failed. ",
      paste(status, collapse = "\n"),
      call. = FALSE
    )
  }

  jsonlite::read_json(output_path, simplifyVector = TRUE)
}

.sn_get_embedding_knn <- function(embeddings, k = 20, n_trees = 50) {
  k <- min(as.integer(k), nrow(embeddings) - 1L)
  if (k < 1L) {
    stop("At least two cells are required to score embedding rarity.", call. = FALSE)
  }

  if (rlang::is_installed("Seurat") && nrow(embeddings) > 1000L) {
    return(.sn_find_annoy_knn(
      embeddings = embeddings,
      k = k,
      n_trees = n_trees,
      include_distance = TRUE
    ))
  }

  .sn_exact_knn(
    embeddings = embeddings,
    k = k,
    include_distance = TRUE
  )
}

.sn_score_embedding_rarity <- function(embeddings, k = 20, n_trees = 50) {
  knn <- .sn_get_embedding_knn(embeddings = embeddings, k = k, n_trees = n_trees)
  rare_score <- rowMeans(knn$dist, na.rm = TRUE)
  rare_score[!is.finite(rare_score)] <- 0
  stats::setNames(as.numeric(rare_score), rownames(embeddings))
}

.sn_run_gapclust <- function(expr, k = 200) {
  check_installed_github("GapClust", "fabotao/GapClust", reason = "to detect rare cells with the GapClust backend.")

  result <- GapClust::GapClust(data = as.matrix(expr), k = as.integer(k))
  cell_ids <- colnames(expr)
  if (length(result) == 1 && is.na(result)) {
    return(data.frame(
      cell_id = cell_ids,
      rare_score = 0,
      rare_cell = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  rare_membership <- sort(unique(unlist(result$rare_cell_indices, use.names = FALSE)))
  rare_score <- apply(result$rare_score, 1, function(x) {
    current <- x[is.finite(x)]
    if (length(current) == 0) {
      return(0)
    }
    max(current)
  })

  data.frame(
    cell_id = cell_ids,
    rare_score = as.numeric(rare_score),
    rare_cell = seq_along(cell_ids) %in% rare_membership,
    stringsAsFactors = FALSE
  )
}

.sn_run_sca <- function(expr,
                        python = NULL,
                        n_comps = 20,
                        iters = 3,
                        nbhd_size = 15,
                        model = "wilcoxon",
                        seed = 717) {
  python <- python %||% getOption("shennong.sca_python", Sys.which("python"))
  if (!nzchar(python)) {
    stop("Could not find a Python executable for the SCA backend.", call. = FALSE)
  }

  workdir <- tempfile("sca_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  expr_path <- file.path(workdir, "expression.csv")
  output_path <- file.path(workdir, "sca_result.json")
  runner_path <- file.path(workdir, "run_sca.py")

  expr_df <- as.data.frame(t(as.matrix(expr)))
  colnames(expr_df) <- rownames(expr)
  rownames(expr_df) <- colnames(expr)
  utils::write.csv(expr_df, file = expr_path, quote = FALSE)

  runner_lines <- c(
    "import json",
    "import pandas as pd",
    "from shannonca.dimred import reduce",
    sprintf("expr = pd.read_csv(%s, index_col=0)", shQuote(expr_path)),
    "reduction = reduce(",
    "    expr.to_numpy(dtype=float),",
    sprintf("    n_comps=%s,", min(as.integer(n_comps), ncol(expr_df) - 1L)),
    sprintf("    iters=%s,", as.integer(iters)),
    sprintf("    nbhd_size=%s,", as.integer(nbhd_size)),
    sprintf("    model=%s,", shQuote(model)),
    sprintf("    seed=%s", as.integer(seed)),
    ")",
    "payload = {'reduction': reduction.tolist()}",
    sprintf("with open(%s, 'w') as handle:", shQuote(output_path)),
    "    json.dump(payload, handle)"
  )
  writeLines(runner_lines, con = runner_path, useBytes = TRUE)

  status <- tryCatch(
    suppressWarnings(system2(
      command = python,
      args = c(shQuote(runner_path)),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) {
      stop(
        "SCA execution failed. Ensure the Python executable and the `shannonca` package are available. ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  exit_code <- attr(status, "status") %||% 0L
  if (!identical(exit_code, 0L) || !file.exists(output_path)) {
    stop(
      "SCA execution failed. Ensure the Python package `shannonca` is installed. ",
      paste(status, collapse = "\n"),
      call. = FALSE
    )
  }

  result <- jsonlite::read_json(output_path, simplifyVector = TRUE)
  embedding <- as.matrix(result$reduction)
  rownames(embedding) <- colnames(expr)
  embedding
}

#' Detect rare cells with native or optional rare-cell backends
#'
#' @param object A \code{Seurat} object.
#' @param method Rare-cell method. Supported values are \code{"gini"},
#'   \code{"sccad"}, \code{"sca"}, \code{"gapclust"}, and
#'   \code{"challenging_groups"}.
#' @param group_by Optional metadata column used with
#'   \code{method = "challenging_groups"}.
#' @param reduction Reduction used by graph-based methods. Defaults to
#'   \code{"harmony"} when present, otherwise \code{"pca"}.
#' @param dims Optional embedding dimensions to use.
#' @param assay Assay used to extract expression values.
#' @param layer Layer used to extract expression values for score-based methods.
#' @param nfeatures Number of rare-aware genes to use for score construction.
#' @param min_cells Minimum number of cells a gene must be detected in before it
#'   is considered by score-based methods.
#' @param max_fraction Maximum expressing-cell fraction for score-based rare
#'   genes.
#' @param threshold Optional explicit threshold on the rare-cell score. When
#'   \code{NULL}, the function uses the upper-IQR rule.
#' @param k Number of neighbors for graph-based methods.
#' @param seed Random seed used by stochastic backends.
#' @param sccad_python Optional Python executable used by the scCAD backend.
#' @param sccad_script Optional path to the upstream \code{scCAD.py} script.
#' @param sccad_normalization Whether scCAD should normalize the provided
#'   matrix internally.
#' @param sccad_rare_h Rare threshold passed to scCAD.
#' @param sccad_merge_h Merge threshold passed to scCAD.
#' @param sccad_overlap_h Overlap threshold passed to scCAD.
#' @param gapclust_k Upper limit of the minor-cluster size used by GapClust.
#' @param sca_python Optional Python executable used by the SCA backend.
#' @param sca_n_comps Number of SCA components used before rarity scoring.
#' @param sca_iters Number of SCA iterations.
#' @param sca_nbhd_size Neighborhood size passed to SCA.
#' @param sca_model Scoring model passed to SCA.
#'
#' @return A data frame with one row per cell, including a \code{rare_score}
#'   column and a logical \code{rare_cell} flag.
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(
#'   Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#' ))
#' rare_tbl <- sn_detect_rare_cells(pbmc, method = "gini")
#' head(rare_tbl)
#' }
#'
#' @export
sn_detect_rare_cells <- function(object,
                                 method = c("gini", "sccad", "sca", "gapclust", "challenging_groups"),
                                 group_by = NULL,
                                 reduction = .sn_default_metric_reduction(object),
                                 dims = NULL,
                                 assay = "RNA",
                                 layer = "data",
                                 nfeatures = 200,
                                 min_cells = 3,
                                 max_fraction = 0.1,
                                 threshold = NULL,
                                 k = 20,
                                 seed = 717,
                                 sccad_python = NULL,
                                 sccad_script = NULL,
                                 sccad_normalization = FALSE,
                                 sccad_rare_h = 0.01,
                                 sccad_merge_h = 0.3,
                                 sccad_overlap_h = 0.7,
                                 gapclust_k = 200,
                                 sca_python = NULL,
                                 sca_n_comps = 20,
                                 sca_iters = 3,
                                 sca_nbhd_size = 15,
                                 sca_model = "wilcoxon") {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.", call. = FALSE)
  }

  method <- rlang::arg_match(method)
  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  cell_ids <- colnames(object)

  if (identical(method, "gini")) {
    gini_tbl <- .sn_detect_rare_features_gini(
      expr = expr,
      nfeatures = nfeatures,
      min_cells = min_cells,
      max_fraction = max_fraction
    )
    if (nrow(gini_tbl) == 0) {
      return(data.frame(
        cell_id = cell_ids,
        method = method,
        rare_score = 0,
        rare_cell = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    feature_mat <- expr[gini_tbl$feature, , drop = FALSE]
    feature_scale <- pmax(Matrix::rowMeans(feature_mat), 1e-8)
    normalized <- Matrix::Diagonal(x = 1 / feature_scale) %*% feature_mat
    rare_score <- Matrix::colMeans(normalized)
  } else if (identical(method, "sccad")) {
    sccad_result <- .sn_run_sccad(
      expr = expr,
      cell_ids = cell_ids,
      gene_ids = rownames(expr),
      python = sccad_python,
      script = sccad_script,
      normalization = sccad_normalization,
      seed = seed,
      rare_h = sccad_rare_h,
      merge_h = sccad_merge_h,
      overlap_h = sccad_overlap_h
    )
    rare_membership <- unique(unlist(sccad_result$rare_sets, use.names = FALSE))
    subcluster_score <- stats::setNames(
      as.numeric(sccad_result$scores),
      unique(sccad_result$sub_clusters)
    )
    rare_score <- unname(subcluster_score[as.character(sccad_result$sub_clusters)])
    rare_flag <- cell_ids %in% rare_membership
    return(data.frame(
      cell_id = cell_ids,
      method = method,
      rare_score = as.numeric(rare_score),
      rare_cell = rare_flag,
      subcluster = as.character(sccad_result$sub_clusters),
      stringsAsFactors = FALSE
    ))
  } else if (identical(method, "gapclust")) {
    return(transform(
      .sn_run_gapclust(expr = expr, k = gapclust_k),
      method = method
    ))
  } else if (identical(method, "sca")) {
    sca_embedding <- .sn_run_sca(
      expr = expr,
      python = sca_python,
      n_comps = sca_n_comps,
      iters = sca_iters,
      nbhd_size = sca_nbhd_size,
      model = sca_model,
      seed = seed
    )
    rare_score <- .sn_score_embedding_rarity(sca_embedding, k = k)
  } else {
    if (is.null(group_by)) {
      stop("`group_by` must be supplied when `method = \"challenging_groups\"`.", call. = FALSE)
    }
    group_tbl <- sn_identify_challenging_groups(
      x = object,
      group_by = group_by,
      reduction = reduction,
      dims = dims,
      k = k,
      neighbor_method = "auto",
      seed = seed
    )
    group_scores <- stats::setNames(group_tbl$challenge_score, group_tbl[[group_by]])
    rare_score <- unname(group_scores[as.character(object[[group_by, drop = TRUE]])])
  }

  score_threshold <- threshold %||% as.numeric(stats::quantile(rare_score, 0.75, na.rm = TRUE) + 1.5 * stats::IQR(rare_score, na.rm = TRUE))
  if (!is.finite(score_threshold)) {
    score_threshold <- Inf
  }

  data.frame(
    cell_id = cell_ids,
    method = method,
    rare_score = as.numeric(rare_score),
    rare_cell = as.numeric(rare_score) >= score_threshold,
    stringsAsFactors = FALSE
  )
}

