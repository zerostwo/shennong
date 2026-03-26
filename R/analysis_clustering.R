# Internal helper to normalize the clustering dimension arguments.
.sn_resolve_cluster_dims <- function(dims = NULL, npcs = 50) {
  dims <- dims %||% seq_len(min(20, npcs))
  dims <- as.integer(dims)
  dims[dims > 0]
}

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
  temp_object <- suppressWarnings(
    Seurat::ScaleData(
      object = object,
      features = features,
      verbose = FALSE
    )
  )
  temp_object <- suppressWarnings(
    Seurat::RunPCA(
      object = temp_object,
      features = features,
      npcs = npcs,
      verbose = FALSE,
      seed.use = 717
    )
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

.sn_detect_rare_features_local_hvg <- function(object,
                                               groups,
                                               rare_groups,
                                               nfeatures = 200) {
  if (length(rare_groups) == 0) {
    return(character(0))
  }

  local_features <- lapply(rare_groups, function(current_group) {
    current_cells <- colnames(object)[groups == current_group]
    if (length(current_cells) < 3) {
      return(character(0))
    }
    current_object <- subset(object, cells = current_cells)
    current_object <- Seurat::FindVariableFeatures(
      current_object,
      nfeatures = nfeatures,
      verbose = FALSE
    )
    Seurat::VariableFeatures(current_object)
  })

  unique(unlist(local_features, use.names = FALSE))
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

  unique(unlist(local_markers, use.names = FALSE))
}

.sn_detect_rare_features_ciara <- function(object,
                                           assay = "RNA",
                                           nfeatures = 200,
                                           min_cells = 3,
                                           max_cells = 100,
                                           k = 20) {
  check_installed("CIARA", reason = "to select rare-cell markers with the CIARA backend.")

  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = "data")
  temp_object <- Seurat::FindNeighbors(
    object = object,
    reduction = "pca",
    dims = seq_len(min(20, ncol(SeuratObject::Embeddings(object, "pca")))),
    k.param = k,
    verbose = FALSE
  )
  nn_graph <- methods::as(temp_object@graphs[[.sn_guess_graph_name(temp_object)]], "dgCMatrix")
  background <- CIARA::get_background_full(
    norm_matrix = as.matrix(expr),
    threshold = 1,
    n_cells_low = min_cells,
    n_cells_high = max_cells
  )
  result <- CIARA::CIARA(
    norm_matrix = as.matrix(expr),
    knn_matrix = as.matrix(nn_graph),
    background = background,
    cores_number = 1,
    p_value = 0.001,
    local_region = 1,
    approximation = FALSE
  )
  ranked <- rownames(result)[order(as.numeric(result[, 1]), rownames(result))]
  utils::head(ranked, nfeatures)
}

.sn_select_rare_features <- function(object,
                                     base_features,
                                     method = "none",
                                     assay = "RNA",
                                     layer = "data",
                                     nfeatures = 200,
                                     group_by = NULL,
                                     rare_group_max_fraction = 0.05,
                                     rare_group_max_cells = 100,
                                     rare_gene_max_fraction = 0.1,
                                     min_cells = 3,
                                     npcs = 20,
                                     dims = 1:10,
                                     resolution = 0.2,
                                     verbose = TRUE) {
  method <- unique(match.arg(method, c("none", "gini", "local_hvg", "local_markers", "ciara"), several.ok = TRUE))
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

  if (any(method %in% c("local_hvg", "local_markers"))) {
    rare_group_info <- .sn_resolve_rare_groups(
      object = object,
      group_by = group_by,
      features = base_features,
      npcs = npcs,
      dims = dims,
      resolution = resolution,
      max_fraction = rare_group_max_fraction,
      max_cells = rare_group_max_cells
    )
  }

  if ("gini" %in% method) {
    gini_tbl <- .sn_detect_rare_features_gini(
      expr = expr,
      nfeatures = nfeatures,
      min_cells = min_cells,
      max_fraction = rare_gene_max_fraction
    )
    rare_features$gini <- gini_tbl$feature
    if (nrow(gini_tbl) > 0) {
      rare_meta[[length(rare_meta) + 1]] <- transform(gini_tbl, method = "gini")
    }
  }

  if ("local_hvg" %in% method) {
    local_hvg <- .sn_detect_rare_features_local_hvg(
      object = object,
      groups = rare_group_info$groups,
      rare_groups = rare_group_info$rare_groups,
      nfeatures = nfeatures
    )
    rare_features$local_hvg <- local_hvg
    if (length(local_hvg) > 0) {
      rare_meta[[length(rare_meta) + 1]] <- data.frame(
        feature = local_hvg,
        score = NA_real_,
        prevalence = NA_integer_,
        prevalence_fraction = NA_real_,
        mean_expression = NA_real_,
        method = "local_hvg",
        stringsAsFactors = FALSE
      )
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
      min_cells = min_cells
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

  if ("ciara" %in% method) {
    ciara_features <- .sn_detect_rare_features_ciara(
      object = object,
      assay = assay,
      nfeatures = nfeatures,
      min_cells = min_cells,
      max_cells = rare_group_max_cells
    )
    rare_features$ciara <- ciara_features
    if (length(ciara_features) > 0) {
      rare_meta[[length(rare_meta) + 1]] <- data.frame(
        feature = ciara_features,
        score = NA_real_,
        prevalence = NA_integer_,
        prevalence_fraction = NA_real_,
        mean_expression = NA_real_,
        method = "ciara",
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
    "rare_sets = [_to_str_list(cluster) for cluster in result]",
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
    system2(
      command = python,
      args = c(shQuote(runner_path)),
      stdout = TRUE,
      stderr = TRUE
    ),
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
    system2(
      command = python,
      args = c(shQuote(runner_path)),
      stdout = TRUE,
      stderr = TRUE
    ),
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
#' @param group Optional metadata column used with
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
#' data("pbmc_small", package = "Shennong")
#' rare_tbl <- sn_detect_rare_cells(pbmc_small, method = "gini")
#' head(rare_tbl)
#' }
#'
#' @export
sn_detect_rare_cells <- function(object,
                                 method = c("gini", "sccad", "sca", "gapclust", "challenging_groups"),
                                 group = NULL,
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
    if (is.null(group)) {
      stop("`group` must be supplied when `method = \"challenging_groups\"`.", call. = FALSE)
    }
    group_tbl <- sn_identify_challenging_groups(
      x = object,
      group = group,
      reduction = reduction,
      dims = dims,
      k = k,
      neighbor_method = "auto",
      seed = seed
    )
    group_scores <- stats::setNames(group_tbl$challenge_score, group_tbl[[group]])
    rare_score <- unname(group_scores[as.character(object[[group, drop = TRUE]])])
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

.sn_select_variable_features <- function(object,
                                         nfeatures = 3000,
                                         split_by = NULL,
                                         verbose = TRUE) {
  if (is_null(split_by) || identical(split_by, "global")) {
    object <- Seurat::FindVariableFeatures(object, nfeatures = nfeatures, verbose = verbose)
    return(list(
      object = object,
      features = Seurat::VariableFeatures(object = object)
    ))
  }

  if (!split_by %in% colnames(object[[]])) {
    stop(glue("`hvg_group_by` must be NULL or a metadata column name. '{split_by}' was not found."))
  }

  split_values <- as.character(object[[split_by, drop = TRUE]])
  split_levels <- unique(split_values)
  feature_lists <- lapply(split_levels, function(current_group) {
    current_cells <- colnames(object)[split_values == current_group]
    current_object <- object[, current_cells]
    current_object <- Seurat::FindVariableFeatures(
      current_object,
      nfeatures = nfeatures,
      verbose = FALSE
    )
    Seurat::VariableFeatures(current_object)
  })

  feature_frequency <- sort(table(unlist(feature_lists, use.names = FALSE)), decreasing = TRUE)
  all_features <- names(feature_frequency)
  feature_ranks <- lapply(feature_lists, function(features) {
    stats::setNames(seq_along(features), features)
  })
  mean_rank <- vapply(all_features, function(feature) {
    mean(vapply(feature_ranks, function(rank_map) {
      if (feature %in% names(rank_map)) {
        rank_map[[feature]]
      } else {
        nfeatures + 1
      }
    }, numeric(1)))
  }, numeric(1))
  selected <- all_features[order(-as.numeric(feature_frequency), mean_rank, all_features)]
  selected <- utils::head(selected, nfeatures)

  Seurat::VariableFeatures(object = object) <- selected
  list(object = object, features = selected)
}

#' Run clustering for a single dataset or Harmony integration workflow
#'
#' This function is the main clustering entry point in `Shennong`.
#' When `batch = NULL`, it performs single-dataset clustering with either the
#' standard Seurat workflow or an SCTransform workflow. When `batch` is
#' supplied, it performs Harmony-based integration followed by clustering
#' and UMAP.
#'
#' @param object A \code{Seurat} object.
#' @param batch A column name in \code{object@meta.data} specifying batch info.
#'   If \code{NULL}, no integration is performed.
#' @param normalization_method One of \code{"seurat"}, \code{"scran"}, or
#'   \code{"sctransform"}. The SCTransform and scran workflows are currently
#'   only supported when \code{batch = NULL}.
#' @param nfeatures Number of variable features to select.
#' @param vars_to_regress Covariates to regress out in \code{ScaleData}.
#' @param resolution Resolution parameter for \code{FindClusters}.
#' @param hvg_group_by Optional metadata column used to compute highly variable
#'   genes within groups before merging and ranking them. Use \code{NULL} to
#'   compute HVGs on the full object.
#' @param rare_feature_method Optional rare-cell-aware feature methods appended
#'   to the base HVG set before PCA/clustering. Supported values are
#'   \code{"none"}, \code{"gini"}, \code{"local_hvg"},
#'   \code{"local_markers"}, and \code{"ciara"}.
#' @param rare_feature_group_by Optional metadata column used to define groups
#'   for rare-group-aware feature extraction. When \code{NULL}, Shennong builds
#'   a temporary coarse clustering from the base HVGs.
#' @param rare_feature_n Number of rare-aware features to add per selected
#'   method.
#' @param rare_group_max_fraction Maximum cluster fraction used when identifying
#'   rare groups for rare-aware feature selection.
#' @param rare_group_max_cells Maximum absolute cluster size used when
#'   identifying rare groups for rare-aware feature selection.
#' @param rare_gene_max_fraction Maximum expressing-cell fraction used by
#'   score-based rare-feature methods such as \code{"gini"}.
#' @param block_genes Either a character vector of predefined bundled signature
#'   categories (for example \code{c("ribo","mito")}) or a custom vector of
#'   gene symbols to exclude from HVGs.
#' @param theta The \code{theta} parameter for \code{harmony::RunHarmony}, controlling batch
#'   diversity preservation vs. correction.
#' @param group_by_vars Optional column name or character vector passed to
#'   \code{harmony::RunHarmony(group.by.vars = ...)}. Defaults to \code{batch}.
#' @param npcs Number of PCs to compute in \code{RunPCA}.
#' @param dims A numeric vector of PCs (dimensions) to use for neighbor search,
#'   clustering, and UMAP.
#' @param species Optional species label. Used when block genes must be resolved
#'   from built-in signatures.
#' @param assay Assay used for clustering. Defaults to \code{"RNA"}.
#' @param layer Layer used as the input count matrix. Defaults to \code{"counts"}.
#' @param return_cluster If \code{TRUE}, return only the cluster assignments.
#' @param verbose Whether to print/log progress messages.
#'
#' @return A \code{Seurat} object with clustering results and embeddings, or a
#'   cluster vector if \code{return_cluster = TRUE}.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- sn_run_cluster(
#'   object = seurat_obj,
#'   normalization_method = "seurat",
#'   resolution = 0.8
#' )
#'
#' seurat_obj <- sn_run_cluster(
#'   object = seurat_obj,
#'   batch = "sample_id",
#'   normalization_method = "seurat",
#'   hvg_group_by = "sample_id",
#'   nfeatures = 3000,
#'   resolution = 0.5,
#'   block_genes = c("ribo", "mito") # or a custom vector of gene symbols
#' )
#' }
#' @export
sn_run_cluster <- function(object,
                           batch = NULL,
                           normalization_method = c("seurat", "scran", "sctransform"),
                           nfeatures = 3000,
                           vars_to_regress = NULL,
                           resolution = 0.8,
                           hvg_group_by = NULL,
                           rare_feature_method = "none",
                           rare_feature_group_by = NULL,
                           rare_feature_n = 200,
                           rare_group_max_fraction = 0.05,
                           rare_group_max_cells = 100,
                           rare_gene_max_fraction = 0.1,
                           block_genes = c("heatshock", "ribo", "mito", "tcr", "immunoglobulins", "pseudogenes"),
                           theta = 2,
                           group_by_vars = NULL,
                           npcs = 50,
                           dims = NULL,
                           species = NULL,
                           assay = "RNA",
                           layer = "counts",
                           return_cluster = FALSE,
                           verbose = TRUE) {
  check_installed("Seurat")
  check_installed("HGNChelper")

  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  normalization_method <- match.arg(normalization_method)
  rare_feature_method <- unique(match.arg(
    rare_feature_method,
    c("none", "gini", "local_hvg", "local_markers", "ciara"),
    several.ok = TRUE
  ))
  if (!is_null(hvg_group_by) && !hvg_group_by %in% colnames(object[[]])) {
    stop(glue("`hvg_group_by` must be NULL or a metadata column name. '{hvg_group_by}' was not found."))
  }
  dims <- .sn_resolve_cluster_dims(dims = dims, npcs = npcs)

  prepared <- .sn_prepare_seurat_analysis_input(
    object = object,
    assay = assay,
    layer = layer
  )
  object <- prepared$object

  if (!is_null(x = batch)) {
    if (!(batch %in% colnames(object@meta.data))) {
      stop(glue("Batch variable '{batch}' not found in metadata."))
    }
    if (normalization_method != "seurat") {
      stop("`normalization_method = \"scran\"` and `\"sctransform\"` are currently only supported when `batch = NULL`.")
    }
    check_installed("harmony")
    if (verbose) .sn_log_info("[sn_run_cluster] Starting integration for batch = '{batch}'.")
  }

  if (verbose) {
    .sn_log_info("[sn_run_cluster] Normalization method = {normalization_method}; batch = {batch %||% 'none'}.")
  }

  predefined_genesets <- c(
    "tcr", "immunoglobulins", "ribo", "mito",
    "heatshock", "noncoding", "pseudogenes", "g1s", "g2m"
  )

  if (normalization_method == "sctransform" && !is_null(block_genes)) {
    .sn_log_warn("`block_genes` is only applied in log-normalization workflows; ignoring it for SCTransform.")
    block_genes <- NULL
  }
  if (normalization_method == "sctransform" && !identical(rare_feature_method, "none")) {
    .sn_log_warn("`rare_feature_method` is only applied in log-normalization workflows; ignoring it for SCTransform.")
    rare_feature_method <- "none"
  }

  if (!is_null(block_genes)) {
    species <- sn_get_species(object = object, species = species)
    if (verbose) .sn_log_info("Processing blocked genes.")

    if (is_character(block_genes) && all(block_genes %in% predefined_genesets)) {
      block_genes <- sn_get_signatures(species = species, category = block_genes)
      if (verbose) {
        .sn_log_info("Loaded {length(block_genes)} blocked genes from built-in sets.")
      }
    } else {
      checked <- HGNChelper::checkGeneSymbols(block_genes, species = species)
      valid_genes <- checked$Suggested.Symbol[!is_na(checked$Suggested.Symbol)]
      invalid <- setdiff(block_genes, valid_genes)

      if (length(invalid) > 0) {
        .sn_log_warn(
          "Removed {length(invalid)} invalid genes from block list (e.g. {paste(utils::head(invalid, 3), collapse=', ')})."
        )
      }
      block_genes <- unique(valid_genes)
      if (verbose) .sn_log_info("Using {length(block_genes)} custom blocked genes.")
    }
  }

  hvg_candidate_nfeatures <- nfeatures
  if (!is_null(block_genes)) {
    hvg_candidate_nfeatures <- min(
      nrow(object),
      nfeatures + length(intersect(block_genes, rownames(object)))
    )
  }

  if (normalization_method == "sctransform") {
    check_installed("glmGamPoi", reason = "for the SCTransform workflow.")

    if (verbose) .sn_log_info("[1/4] Running SCTransform.")
    object <- Seurat::SCTransform(
      object = object,
      variable.features.n = nfeatures,
      vars.to.regress = vars_to_regress,
      verbose = verbose,
      seed.use = 717
    )

    if (verbose) .sn_log_info("[2/4] Running PCA.")
    object <- Seurat::RunPCA(
      object,
      npcs = npcs,
      verbose = verbose,
      seed.use = 717
    )

    reduction <- "pca"
  } else {
    if (normalization_method == "scran") {
      object <- sn_normalize_data(
        object = object,
        method = "scran",
        assay = assay,
        layer = "counts"
      )
    } else {
      object <- Seurat::NormalizeData(object = object, verbose = verbose)
    }

    if (!is_null(species)) {
      if (verbose) .sn_log_info("[1/6] Scoring cell cycle.")
      object <- sn_score_cell_cycle(object = object, species = species)
    }

    if (verbose) {
      .sn_log_info("[2/6] Selecting highly variable features with hvg_group_by = {hvg_group_by %||% 'global'}.")
    }
    hvg_info <- .sn_select_variable_features(
      object = object,
      nfeatures = hvg_candidate_nfeatures,
      split_by = hvg_group_by,
      verbose = verbose
    )
    object <- hvg_info$object
    hvg <- hvg_info$features

    if (!is_null(block_genes)) {
      n_before <- length(hvg)
      hvg <- setdiff(hvg, block_genes)
      n_after_filter <- length(hvg)
      hvg <- utils::head(hvg, nfeatures)
      n_removed <- n_before - n_after_filter

      if (verbose) {
        .sn_log_info(
          "  Removed {n_removed} genes ({round(n_removed/n_before*100,1)}%) via block_genes"
        )
      }

      if (length(hvg) < nfeatures) {
        .sn_log_warn(
          "Only {length(hvg)} HVGs left (< requested {nfeatures}).\n",
          "Consider adjusting 'nfeatures' or 'block_genes'."
        )
      }
      hvg <- utils::head(hvg, nfeatures)
    }

    rare_feature_info <- .sn_select_rare_features(
      object = object,
      base_features = hvg,
      method = rare_feature_method,
      assay = assay,
      layer = "data",
      nfeatures = rare_feature_n,
      group_by = rare_feature_group_by,
      rare_group_max_fraction = rare_group_max_fraction,
      rare_group_max_cells = rare_group_max_cells,
      rare_gene_max_fraction = rare_gene_max_fraction,
      min_cells = 3,
      npcs = min(npcs, 20),
      dims = seq_len(min(10, npcs)),
      resolution = min(resolution, 0.4),
      verbose = verbose
    )
    selected_rare_features <- utils::head(rare_feature_info$features, rare_feature_n)
    if (length(selected_rare_features) > 0) {
      hvg <- unique(c(hvg, selected_rare_features))
    }
    Seurat::VariableFeatures(object = object) <- hvg
    object@misc$rare_feature_selection <- list(
      method = rare_feature_method,
      base_hvg_n = nfeatures,
      rare_feature_n = length(selected_rare_features),
      selected_features = hvg,
      rare_features = selected_rare_features,
      rare_feature_table = rare_feature_info$metadata,
      rare_groups = if (!is.null(rare_feature_info$groups)) rare_feature_info$groups$rare_groups else character(0)
    )

    #-- Scaling
    if (verbose) .sn_log_info("[3/6] Scaling data.")
    object <- Seurat::ScaleData(
      object = object,
      vars.to.regress = vars_to_regress,
      features = hvg,
      verbose = verbose
    )

    #-- PCA
    if (verbose) .sn_log_info("[4/6] Running PCA.")
    object <- Seurat::RunPCA(
      object,
      npcs = npcs,
      features = hvg,
      verbose = verbose,
      seed.use = 717
    )

    if (is_null(x = batch)) {
      reduction <- "pca"
    } else {
      if (verbose) .sn_log_info("[5/6] Running Harmony integration.")
      group_by_vars <- group_by_vars %||% batch
      object <- harmony::RunHarmony(
        object = object,
        group.by.vars = group_by_vars,
        theta = theta,
        reduction.use = "pca",
        verbose = verbose
      )
      reduction <- "harmony"
    }
  }

  if (verbose) .sn_log_info("[6/6] Clustering with integrated embeddings.")
  object <- Seurat::FindNeighbors(object, reduction = reduction, dims = dims, verbose = verbose)
  object <- Seurat::FindClusters(object, resolution = resolution, random.seed = 717, verbose = verbose)

  if (return_cluster) {
    object <- .sn_restore_seurat_analysis_input(object = object, context = prepared$context)
    if (verbose) .sn_log_info("Integration completed successfully.")
    return(object@meta.data[, "seurat_clusters"])
  } else {
    if (verbose) .sn_log_info("[7/7] Running UMAP.")
    object <- suppressWarnings(Seurat::RunUMAP(
      object,
      reduction = reduction,
      dims = dims,
      umap.method = "uwot",
      metric = "cosine",
      verbose = verbose,
      seed.use = 717
    ))
    object <- .sn_restore_seurat_analysis_input(object = object, context = prepared$context)
    if (verbose) .sn_log_info("Integration completed successfully.")
    return(.sn_log_seurat_command(object = object, assay = assay, name = "sn_run_cluster"))
  }
}


#' Run CellTypist for automated cell type annotation
#'
#' @param x A Seurat object or a path to a count matrix / AnnData file that CellTypist can consume.
#' @param celltypist Path to the `celltypist` binary. Defaults to "/opt/mambaforge/envs/scverse/bin/celltypist".
#' @param model Model used for predictions. Defaults to "Immune_All_Low.pkl".
#' @param outdir Directory to store the output files. If NULL, use a temporary directory.
#' @param prefix Prefix for the output files. By default, use the model name plus a dot.
#' @param mode Choose the cell type with the largest score/probability (`"best_match"`) or enable multi-label classification (`"prob_match"`).
#' @param p_thres Probability threshold for the multi-label classification. Ignored if `mode = "best_match"`.
#' @param majority_voting Logical. Whether to refine labels using majority voting after over-clustering.
#' @param over_clustering Input file or a string key specifying an existing metadata column in the AnnData object, or "auto".
#' @param min_prop For the dominant cell type within a subcluster, the minimum proportion of cells required to name the subcluster by this cell type.
#' @param transpose_input Logical. If `TRUE`, add the `--transpose-input` argument when calling `celltypist`.
#' @param gene_file If the provided input is in the `mtx` format, path to the file storing gene information. Otherwise ignored.
#' @param cell_file If the provided input is in the `mtx` format, path to the file storing cell information. Otherwise ignored.
#' @param assay Assay used when exporting Seurat counts to CellTypist. Defaults to \code{"RNA"}.
#' @param layer Layer used as the input count matrix. Defaults to \code{"counts"}.
#' @param xlsx Logical. If `TRUE`, merge output tables into a single Excel (.xlsx). Defaults to `FALSE`.
#' @param plot_results Logical. If `TRUE`, plot the prediction results. Defaults to `FALSE`.
#' @param quiet Logical. If `TRUE`, hide the banner and config info from `celltypist`. Defaults to `FALSE`.
#'
#' @return A Seurat object with three new columns in its metadata: `_predicted_labels`, `_over_clustering`, `_majority_voting`.
#'
#' @examples
#' \dontrun{
#' data("pbmc_small", package = "Shennong")
#' pbmc <- sn_run_cluster(pbmc, normalization_method = "seurat", verbose = FALSE)
#' pbmc <- sn_run_celltypist(pbmc, model = "Immune_All_Low.pkl")
#' head(colnames(pbmc[[]]))
#' }
#' @export
sn_run_celltypist <- function(x,
                              celltypist = NULL,
                              model = "Immune_All_Low.pkl",
                              outdir = NULL,
                              prefix = NULL,
                              mode = c("best_match", "prob_match"),
                              p_thres = 0.5,
                              majority_voting = TRUE,
                              over_clustering = "auto",
                              min_prop = 0,
                              transpose_input = TRUE,
                              gene_file = NULL,
                              cell_file = NULL,
                              assay = "RNA",
                              layer = "counts",
                              xlsx = FALSE,
                              plot_results = FALSE,
                              quiet = FALSE) {
  check_installed(c("Seurat", "data.table", "logger", "glue"),
    reason = "to run CellTypist analysis."
  )

  mode <- match.arg(mode)
  celltypist <- celltypist %||% getOption("shennong.celltypist_path", Sys.which("celltypist"))
  sn_check_file(celltypist)

  .sn_log_info("Starting CellTypist analysis with model = {model}.")
  tictoc::tic("Total CellTypist runtime")

  if (is_null(outdir)) {
    outdir <- tempfile("celltypist_")
    dir.create(outdir, recursive = TRUE)
    .sn_log_info("Using a temporary output directory: {outdir}.")
  } else {
    outdir <- sn_set_path(outdir)
    .sn_log_info("Using the user-specified output directory: {outdir}.")
  }

  on.exit(
    {
      if (dir.exists(outdir) && grepl("^celltypist_", basename(outdir))) {
        unlink(outdir, recursive = TRUE)
        log_debug("Cleaned temporary directory: {outdir}")
      }
    },
    add = TRUE
  )

  if (inherits(x, "Seurat")) {
    .sn_log_info("Converting the Seurat object to CellTypist input format.")
    input_data <- file.path(outdir, "counts.csv")

    counts <- tryCatch(
      .sn_get_seurat_layer_data(object = x, assay = assay, layer = layer),
      error = function(e) {
        .sn_log_error("Failed to extract the counts layer: {e$message}.")
        stop("Counts layer extraction failed")
      }
    )

    data.table::fwrite(
      x = as.data.frame(counts),
      file = input_data,
      quote = FALSE,
      row.names = TRUE,
      showProgress = FALSE
    )
    log_debug("Count matrix written to {input_data} ({file.size(input_data)} bytes)")
  } else {
    input_data <- x
    .sn_log_info("Using the precomputed input matrix: {input_data}.")
  }

  over_clustering_path <- NULL

  if (over_clustering != "auto" && over_clustering %in% colnames(x@meta.data)) {
    .sn_log_info("Using the existing clustering column: {over_clustering}.")
    over_clustering_path <- file.path(outdir, "over_clustering.txt")
    writeLines(
      as.character(x[[over_clustering, drop = TRUE]]),
      over_clustering_path
    )
  } else if (over_clustering == "auto") {
    .sn_log_warn("Automatic over-clustering is not implemented yet.")
  }

  model_name <- tools::file_path_sans_ext(basename(model))
  prefix <- prefix %||% glue("{model_name}.")
  predicted_labels_path <- file.path(outdir, glue("{prefix}predicted_labels.csv"))

  cmd_args <- c(
    "--indata", shQuote(input_data),
    "--model", shQuote(model),
    "--mode", mode,
    "--outdir", shQuote(outdir),
    "--prefix", prefix
  )

  add_arg <- function(args, flag, value, condition = TRUE) {
    if (condition && !is_null(value)) c(args, flag, shQuote(value)) else args
  }

  cmd_args <- add_arg(cmd_args, "--gene-file", gene_file, !is_null(gene_file))
  cmd_args <- add_arg(cmd_args, "--cell-file", cell_file, !is_null(cell_file))
  cmd_args <- add_arg(cmd_args, "--p-thres", p_thres, mode == "prob_match")
  cmd_args <- add_arg(cmd_args, "--over-clustering", over_clustering_path, !is_null(over_clustering_path))
  cmd_args <- add_arg(cmd_args, "--min-prop", min_prop, majority_voting)

  if (transpose_input) cmd_args <- c(cmd_args, "--transpose-input")
  if (xlsx) cmd_args <- c(cmd_args, "--xlsx")
  if (plot_results) cmd_args <- c(cmd_args, "--plot-results")
  if (quiet) cmd_args <- c(cmd_args, "--quiet")
  if (majority_voting) cmd_args <- c(cmd_args, "--majority-voting")
  .sn_log_info("Executing CellTypist with command:\n{celltypist} {paste(cmd_args, collapse = ' ')}")

  exit_code <- system2(
    command = celltypist,
    args = cmd_args,
    stdout = if (quiet) FALSE else "",
    stderr = if (quiet) FALSE else ""
  )

  if (exit_code != 0) {
    .sn_log_error("CellTypist failed with exit code {exit_code}.")
    stop("CellTypist execution failed. Check logs for details.")
  }
  log_debug("CellTypist output written to {outdir}")

  if (!file.exists(predicted_labels_path)) {
    .sn_log_error("Prediction output is missing: {predicted_labels_path}.")
    stop("CellTypist did not generate expected output files")
  }

  predicted_labels <- utils::read.csv(predicted_labels_path, row.names = 1)
  if (majority_voting) {
    colnames(predicted_labels) <- c(
      paste0(gsub("\\.pkl$", "", model_name), "_predicted_labels"),
      paste0(gsub("\\.pkl$", "", model_name), "_over_clustering"),
      paste0(gsub("\\.pkl$", "", model_name), "_majority_voting")
    )
  } else {
    colnames(predicted_labels) <- c(
      paste0(gsub("\\.pkl$", "", model_name), "_predicted_labels")
    )
  }

  .sn_log_info("Adding {ncol(predicted_labels)} metadata columns to the Seurat object.")
  x <- SeuratObject::AddMetaData(x, metadata = predicted_labels)

  tictoc::toc()
  .sn_log_info("CellTypist analysis completed successfully.")

  .sn_log_seurat_command(object = x, assay = assay, name = "sn_run_celltypist")
}
