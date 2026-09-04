#' Identify rare or difficult-to-separate groups
#'
#' This helper summarizes group size, local neighbor purity, graph connectivity,
#' and silhouette width for a grouping column. It is intended to surface rare
#' populations or poorly separated groups that may be obscured in standard UMAP
#' visual inspection.
#'
#' @param x A Seurat object.
#' @param group_by Metadata column defining the groups to inspect.
#' @param graph Optional graph name stored in \code{x@graphs}. If \code{NULL},
#'   the function tries to reuse an existing nearest-neighbor graph and falls
#'   back to a kNN graph built from the selected embedding.
#' @param reduction Reduction name used when a graph must be built from
#'   embeddings. Defaults to \code{"harmony"} when present, otherwise
#'   \code{"pca"}.
#' @param dims Optional integer vector of embedding dimensions to retain.
#' @param cells Optional character vector of cell names to include.
#' @param k Number of neighbors used when building a graph from embeddings.
#' @param neighbor_method Strategy used when a graph must be built. One of
#'   \code{"auto"}, \code{"graph"}, \code{"annoy"}, or \code{"exact"}.
#' @param max_cells Optional integer cap used to subsample cells before running
#'   the metric.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling. Defaults to \code{group_by}.
#' @param rare_fraction Fraction-of-cells threshold for flagging rare groups.
#' @param rare_n Absolute cell-count threshold for flagging rare groups.
#' @param challenge_threshold Threshold on the derived \code{challenge_score}
#'   used to flag difficult groups.
#' @param seed Random seed used when \code{max_cells} triggers subsampling.
#' @param n_trees Number of Annoy trees when \code{neighbor_method = "annoy"}.
#'
#' @return A data frame with one row per group and rarity/separation flags.
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(
#'   Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#' ))
#' pbmc <- sn_run_cluster(
#'   pbmc,
#'   batch = "sample",
#'   species = "human",
#'   verbose = FALSE
#' )
#' difficult_tbl <- sn_identify_challenging_groups(
#'   pbmc,
#'   group_by = "seurat_clusters",
#'   reduction = "harmony"
#' )
#' difficult_tbl
#' }
#'
#' @export
sn_identify_challenging_groups <- function(
  x,
  group_by = NULL,
  graph = NULL,
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  cells = NULL,
  k = 20,
  neighbor_method = c("auto", "graph", "annoy", "exact"),
  max_cells = 3000,
  stratify_by = NULL,
  rare_fraction = 0.02,
  rare_n = 50,
  challenge_threshold = 0.5,
  seed = 717,
  n_trees = 50
) {
  if (is.null(group_by)) {
    stop("`group_by` must be supplied.", call. = FALSE)
  }
  stratify_by <- stratify_by %||% group_by
  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = cells,
    max_cells = max_cells,
    stratify_by = stratify_by,
    seed = seed,
    required_cols = group_by
  )
  neighbor_method <- rlang::arg_match(neighbor_method)

  graph_info <- .sn_get_metric_graph(
    x = x,
    embeddings = metric_input$embeddings,
    cells = metric_input$cells,
    graph = graph,
    k = k,
    neighbor_method = neighbor_method,
    n_trees = n_trees
  )

  groups <- metric_input$metadata[[group_by]]
  purity_tbl <- .sn_calculate_neighbor_purity(
    adjacency = graph_info$adjacency,
    labels = groups,
    label_name = group_by,
    cell_ids = metric_input$cells
  )
  connectivity_tbl <- .sn_calculate_connectivity_table(
    adjacency = graph_info$adjacency,
    groups = groups,
    group_name = group_by
  )

  silhouette_tbl <- tryCatch(
    .sn_compute_silhouette_table(
      embeddings = metric_input$embeddings,
      labels = groups,
      cell_ids = metric_input$cells,
      label_name = group_by
    ),
    error = function(e) NULL
  )

  group_levels <- unique(as.character(groups))
  group_n <- table(as.character(groups))
  total_cells <- length(groups)

  difficult_tbl <- lapply(group_levels, function(current_group) {
    current_cells <- metric_input$cells[as.character(groups) == current_group]
    purity_values <- purity_tbl$neighbor_purity[purity_tbl$cell_id %in% current_cells]
    current_connectivity <- connectivity_tbl$connectivity_score[
      connectivity_tbl[[group_by]] == current_group
    ]
    current_silhouette <- NA_real_
    if (!is.null(silhouette_tbl)) {
      current_silhouette <- mean(
        silhouette_tbl$silhouette_width[
          silhouette_tbl$cell_id %in% current_cells
        ]
      )
    }

    n_cells <- as.integer(group_n[[current_group]])
    fraction_cells <- n_cells / total_cells
    purity_score <- stats::median(purity_values, na.rm = TRUE)
    silhouette_score <- .sn_scale_silhouette(current_silhouette)
    connectivity_score <- current_connectivity[[1]]
    resolution_score <- mean(
      c(purity_score, silhouette_score, connectivity_score),
      na.rm = TRUE
    )
    challenge_score <- 1 - resolution_score

    data.frame(
      group = current_group,
      n_cells = n_cells,
      fraction_cells = fraction_cells,
      median_neighbor_purity = purity_score,
      mean_neighbor_purity = mean(purity_values, na.rm = TRUE),
      graph_connectivity = connectivity_score,
      mean_silhouette = current_silhouette,
      separation_score = resolution_score,
      challenge_score = challenge_score,
      rare_group = n_cells < rare_n || fraction_cells < rare_fraction,
      challenging_group = is.finite(challenge_score) &&
        challenge_score >= challenge_threshold,
      stringsAsFactors = FALSE
    )
  })

  difficult_tbl <- .sn_bind_rows(difficult_tbl)
  names(difficult_tbl)[names(difficult_tbl) == "group"] <- group_by
  difficult_tbl <- difficult_tbl[order(difficult_tbl$challenge_score, decreasing = TRUE), ]
  rownames(difficult_tbl) <- NULL

  attr(difficult_tbl, "graph_source") <- graph_info$source
  difficult_tbl
}

#' Assess integration quality across multiple metrics
#'
#' This wrapper combines fast, broadly useful metrics for batch_by mixing,
#' biological conservation, and cluster-level diagnostics. It reuses existing
#' neighbor graphs when available and otherwise builds a kNN graph with Annoy or
#' an exact fallback.
#'
#' @param x A Seurat object.
#' @param batch_by Metadata column containing batch labels.
#' @param label_by Optional metadata column containing biological labels such as
#'   cell type annotations.
#' @param cluster_by Optional metadata column containing the current clustering.
#' @param reduction Reduction name used for the primary integrated embedding.
#'   Defaults to \code{"harmony"} when present, otherwise \code{"pca"}.
#' @param baseline Optional Seurat object used as the baseline reference for
#'   PCR batch_by scoring. If \code{NULL}, the baseline is taken from \code{x}.
#' @param baseline_reduction Optional reduction name used as the baseline.
#' @param graph Optional graph name stored in \code{x@graphs}.
#' @param dims Optional integer vector of embedding dimensions to retain.
#' @param k Number of neighbors used when a graph must be built from embeddings.
#' @param metrics Optional character vector selecting which metrics to compute.
#'   Supported values are \code{"batch_lisi"}, \code{"label_lisi"},
#'   \code{"batch_silhouette"}, \code{"label_silhouette"},
#'   \code{"graph_connectivity"}, \code{"clustering_agreement"},
#'   \code{"isolated_label_score"}, \code{"cluster_purity"},
#'   \code{"cluster_entropy"}, \code{"pcr_batch"}, and
#'   \code{"challenging_groups"}.
#' @param neighbor_method Strategy used when a graph must be built. One of
#'   \code{"auto"}, \code{"graph"}, \code{"annoy"}, or \code{"exact"}.
#' @param max_cells Optional integer cap used to subsample cells before running
#'   the metrics.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling. Defaults to \code{batch_by}.
#' @param rare_fraction Fraction-of-cells threshold for flagging rare groups.
#' @param rare_n Absolute cell-count threshold for flagging rare groups.
#' @param challenge_threshold Threshold on the derived \code{challenge_score}
#'   used to flag difficult groups.
#' @param seed Random seed used when \code{max_cells} triggers subsampling.
#' @param n_trees Number of Annoy trees when \code{neighbor_method = "annoy"}.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#'
#' @return A list with four top-level elements:
#'   \itemize{
#'     \item \code{summary}: one row per summary metric, plus aggregate scores.
#'     \item \code{per_cell}: per-cell outputs such as LISI and silhouette.
#'     \item \code{per_group}: per-group outputs such as connectivity,
#'           composition, and difficult-group diagnostics.
#'     \item \code{parameters}: the effective runtime configuration.
#'   }
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(
#'   Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#' ))
#' pbmc <- sn_run_cluster(
#'   pbmc,
#'   batch = "sample",
#'   species = "human",
#'   verbose = FALSE
#' )
#' metrics <- sn_assess_integration(
#'   pbmc,
#'   batch = "sample",
#'   cluster_by = "seurat_clusters",
#'   reduction = "harmony",
#'   baseline_reduction = "pca"
#' )
#' metrics$summary
#' }
#'
#' @export
sn_assess_integration <- function(
  x,
  batch_by = NULL,
  label_by = NULL,
  cluster_by = NULL,
  reduction = .sn_default_metric_reduction(x),
  baseline = NULL,
  baseline_reduction = .sn_default_baseline_reduction(x, reduction),
  graph = NULL,
  dims = NULL,
  k = 20,
  metrics = NULL,
  neighbor_method = c("auto", "graph", "annoy", "exact"),
  max_cells = 5000,
  stratify_by = NULL,
  rare_fraction = 0.02,
  rare_n = 50,
  challenge_threshold = 0.5,
  seed = 717,
  n_trees = 50,
  object = NULL
) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  cluster_by <- cluster_by %||% "seurat_clusters"
  if (is.null(batch_by)) {
    stop("`batch_by` must be supplied.", call. = FALSE)
  }
  stratify_by <- stratify_by %||% batch_by
  neighbor_method <- rlang::arg_match(neighbor_method)
  metrics <- metrics %||% c(
    "batch_lisi",
    "label_lisi",
    "batch_silhouette",
    "label_silhouette",
    "graph_connectivity",
    "clustering_agreement",
    "isolated_label_score",
    "cluster_purity",
    "cluster_entropy",
    "pcr_batch",
    "challenging_groups"
  )

  required_cols <- unique(stats::na.omit(c(batch_by, label_by, stratify_by)))
  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = NULL,
    max_cells = max_cells,
    stratify_by = stratify_by,
    seed = seed,
    required_cols = required_cols
  )

  summary_rows <- list()
  per_cell <- list()
  per_group <- list()
  notes <- character()
  cell_count <- nrow(metric_input$metadata)

  if ("batch_silhouette" %in% metrics) {
    batch_silhouette <- .sn_calculate_batch_silhouette_score(
      embeddings = metric_input$embeddings,
      batch = metric_input$metadata[[batch_by]],
      label = if (!is.null(label_by)) metric_input$metadata[[label_by]] else NULL,
      cell_ids = metric_input$cells,
      batch_name = batch_by,
      label_name = label_by
    )
    summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
      metric = "batch_silhouette",
      category = "batch_removal",
      score = batch_silhouette$raw_score,
      scaled_score = batch_silhouette$scaled_score,
      n_cells = cell_count,
      source = if (is.null(label_by)) reduction else paste0(reduction, " + ", label_by),
      note = batch_silhouette$note
    )
    per_cell$batch_silhouette <- batch_silhouette$per_cell
    if (!is.null(batch_silhouette$per_group)) {
      per_group$batch_silhouette <- batch_silhouette$per_group
    }
  }

  if ("label_silhouette" %in% metrics && !is.null(label_by)) {
    label_silhouette <- .sn_compute_silhouette_table(
      embeddings = metric_input$embeddings,
      labels = metric_input$metadata[[label_by]],
      cell_ids = metric_input$cells,
      label_name = label_by
    )
    summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
      metric = "label_silhouette",
      category = "biology_conservation",
      score = mean(label_silhouette$silhouette_width),
      scaled_score = .sn_scale_silhouette(mean(label_silhouette$silhouette_width)),
      n_cells = nrow(label_silhouette),
      source = reduction
    )
    per_cell$label_silhouette <- label_silhouette
  }

  if ("batch_lisi" %in% metrics) {
    if (rlang::is_installed("lisi")) {
      batch_lisi <- sn_calculate_lisi(
        x = x,
        reduction = reduction,
        label_by = batch_by,
        dims = dims,
        cells = metric_input$cells,
        max_cells = NULL,
        stratify_by = NULL,
        seed = seed
      )
      batch_levels <- length(unique(metric_input$metadata[[batch_by]]))
      batch_mean <- mean(batch_lisi[[batch_by]])
      summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
        metric = "batch_lisi",
        category = "batch_removal",
        score = batch_mean,
        scaled_score = .sn_scale_lisi(batch_mean, n_levels = batch_levels),
        n_cells = nrow(batch_lisi),
        source = reduction
      )
      per_cell$batch_lisi <- batch_lisi
    } else {
      notes <- c(
        notes,
        "Skipped batch_lisi because package 'lisi' is not installed."
      )
    }
  }

  if ("label_lisi" %in% metrics && !is.null(label_by)) {
    if (rlang::is_installed("lisi")) {
      label_lisi <- sn_calculate_lisi(
        x = x,
        reduction = reduction,
        label_by = label_by,
        dims = dims,
        cells = metric_input$cells,
        max_cells = NULL,
        stratify_by = NULL,
        seed = seed
      )
      label_levels <- length(unique(metric_input$metadata[[label_by]]))
      label_mean <- mean(label_lisi[[label_by]])
      summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
        metric = "label_lisi",
        category = "biology_conservation",
        score = label_mean,
        scaled_score = .sn_scale_lisi(
          label_mean,
          n_levels = label_levels,
          inverse = TRUE
        ),
        n_cells = nrow(label_lisi),
        source = reduction
      )
      per_cell$label_lisi <- label_lisi
    } else {
      notes <- c(
        notes,
        "Skipped label_lisi because package 'lisi' is not installed."
      )
    }
  }

  if ("graph_connectivity" %in% metrics) {
    connectivity_label <- if (!is.null(label_by) && label_by %in% colnames(metric_input$metadata)) {
      label_by
    } else if (!is.null(cluster_by) && cluster_by %in% colnames(metric_input$metadata)) {
      cluster_by
    } else {
      NULL
    }
    if (!is.null(connectivity_label)) {
      connectivity_tbl <- sn_calculate_graph_connectivity(
        x = x,
        label_by = connectivity_label,
        graph = graph,
        reduction = reduction,
        dims = dims,
        cells = metric_input$cells,
        k = k,
        neighbor_method = neighbor_method,
        max_cells = NULL,
        stratify_by = NULL,
        seed = seed,
        n_trees = n_trees
      )
      connectivity_category <- if (!is.null(label_by)) {
        "biology_conservation"
      } else {
        "structure"
      }
      summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
        metric = if (!is.null(label_by)) {
          "graph_connectivity"
        } else {
          "cluster_connectivity"
        },
        category = connectivity_category,
        score = mean(connectivity_tbl$connectivity_score),
        scaled_score = mean(connectivity_tbl$connectivity_score),
        n_cells = cell_count,
        source = attr(connectivity_tbl, "graph_source")
      )
      per_group$graph_connectivity <- connectivity_tbl
    }
  }

  if ("clustering_agreement" %in% metrics &&
    !is.null(label_by) &&
    !is.null(cluster_by) &&
    label_by %in% colnames(metric_input$metadata) &&
    cluster_by %in% colnames(metric_input$metadata)) {
    agreement_tbl <- sn_calculate_clustering_agreement(
      x = metric_input$metadata,
      cluster_by = cluster_by,
      label_by = label_by
    )
    summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
      metric = "ari",
      category = "biology_conservation",
      score = agreement_tbl$ari,
      scaled_score = agreement_tbl$ari,
      n_cells = agreement_tbl$n_cells,
      source = paste(cluster_by, "vs", label_by)
    )
    summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
      metric = "nmi",
      category = "biology_conservation",
      score = agreement_tbl$nmi,
      scaled_score = agreement_tbl$nmi,
      n_cells = agreement_tbl$n_cells,
      source = paste(cluster_by, "vs", label_by)
    )
    per_group$clustering_agreement <- agreement_tbl
  }

  if ("isolated_label_score" %in% metrics &&
    !is.null(label_by) &&
    label_by %in% colnames(metric_input$metadata)) {
    isolated_tbl <- sn_calculate_isolated_label_score(
      x = x,
      label_by = label_by,
      reduction = reduction,
      dims = dims,
      cells = metric_input$cells,
      max_cells = NULL,
      stratify_by = NULL,
      isolated_fraction = rare_fraction,
      isolated_n = rare_n,
      seed = seed
    )
    isolated_score <- attr(isolated_tbl, "overall_score")
    if (!is.na(isolated_score)) {
      summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
        metric = "isolated_label_score",
        category = "biology_conservation",
        score = isolated_score,
        scaled_score = isolated_score,
        n_cells = sum(isolated_tbl$n_cells[isolated_tbl$isolated_label]),
        source = label_by,
        note = paste(length(attr(isolated_tbl, "isolated_labels")), "isolated labels")
      )
    } else {
      notes <- c(
        notes,
        paste0("Skipped isolated_label_score because no labels met the isolated thresholds for '", label_by, "'.")
      )
    }
    per_group$isolated_label_score <- isolated_tbl
  }

  if ("cluster_purity" %in% metrics &&
    !is.null(cluster_by) &&
    !is.null(label_by) &&
    cluster_by %in% colnames(metric_input$metadata) &&
    label_by %in% colnames(metric_input$metadata)) {
    purity_tbl <- sn_calculate_cluster_purity(
      x = metric_input$metadata,
      cluster_by = cluster_by,
      label_by = label_by
    )
    summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
      metric = "cluster_label_purity",
      category = "biology_conservation",
      score = mean(purity_tbl$purity_score),
      scaled_score = mean(purity_tbl$purity_score),
      n_cells = sum(purity_tbl$n_cells),
      source = paste(cluster_by, "vs", label_by)
    )
    per_group$cluster_purity <- purity_tbl
  }

  if ("cluster_entropy" %in% metrics &&
    !is.null(cluster_by) &&
    cluster_by %in% colnames(metric_input$metadata) &&
    batch_by %in% colnames(metric_input$metadata)) {
    entropy_tbl <- sn_calculate_cluster_entropy(
      x = metric_input$metadata,
      cluster_by = cluster_by,
      label_by = batch_by
    )
    summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
      metric = "cluster_batch_entropy",
      category = "batch_removal",
      score = mean(entropy_tbl$normalized_entropy),
      scaled_score = mean(entropy_tbl$normalized_entropy),
      n_cells = sum(entropy_tbl$n_cells),
      source = paste(cluster_by, "vs", batch_by)
    )
    per_group$cluster_entropy <- entropy_tbl
  }

  if ("pcr_batch" %in% metrics) {
    pcr_tbl <- sn_calculate_pcr_batch(
      x = x,
      batch_by = batch_by,
      reduction = reduction,
      dims = dims,
      baseline = baseline,
      baseline_reduction = baseline_reduction,
      cells = metric_input$cells,
      max_cells = NULL,
      stratify_by = NULL,
      seed = seed
    )
    summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
      metric = "pcr_batch",
      category = "batch_removal",
      score = pcr_tbl$batch_variance,
      scaled_score = pcr_tbl$scaled_score,
      n_cells = pcr_tbl$n_cells,
      source = if (is.na(pcr_tbl$baseline_reduction)) {
        reduction
      } else {
        paste(reduction, "vs", pcr_tbl$baseline_reduction)
      }
    )
    per_group$pcr_batch <- pcr_tbl
  }

  if ("challenging_groups" %in% metrics) {
    difficulty_group <- if (!is.null(cluster_by) && cluster_by %in% colnames(metric_input$metadata)) {
      cluster_by
    } else if (!is.null(label_by) && label_by %in% colnames(metric_input$metadata)) {
      label_by
    } else {
      NULL
    }
    if (!is.null(difficulty_group)) {
      difficult_tbl <- sn_identify_challenging_groups(
        x = x,
        group_by = difficulty_group,
        graph = graph,
        reduction = reduction,
        dims = dims,
        cells = metric_input$cells,
        k = k,
        neighbor_method = neighbor_method,
        max_cells = NULL,
        stratify_by = NULL,
        rare_fraction = rare_fraction,
        rare_n = rare_n,
        challenge_threshold = challenge_threshold,
        seed = seed,
        n_trees = n_trees
      )
      summary_rows[[length(summary_rows) + 1]] <- .sn_metric_row(
        metric = "well_resolved_group_fraction",
        category = "structure",
        score = mean(!difficult_tbl$challenging_group),
        scaled_score = mean(!difficult_tbl$challenging_group),
        n_cells = cell_count,
        source = difficulty_group,
        note = paste0(sum(difficult_tbl$rare_group), " rare groups flagged")
      )
      per_group$challenging_groups <- difficult_tbl
    }
  }

  if (!is.null(cluster_by) &&
    cluster_by %in% colnames(metric_input$metadata) &&
    batch_by %in% colnames(metric_input$metadata)) {
    per_group$composition <- sn_calculate_composition(
      x = metric_input$metadata,
      group_by = batch_by,
      variable = cluster_by,
      min_cells = 1
    )
  }

  summary_tbl <- .sn_bind_rows(summary_rows)
  aggregate_rows <- .sn_calculate_assessment_aggregates(summary_tbl)
  summary_tbl <- .sn_bind_rows(list(summary_tbl, aggregate_rows))
  rownames(summary_tbl) <- NULL

  result <- list(
    summary = summary_tbl,
    per_cell = per_cell,
    per_group = per_group,
    parameters = list(
      reduction = reduction,
      baseline_reduction = baseline_reduction,
      batch_by = batch_by,
      label_by = label_by,
      cluster_by = cluster_by,
      dims = dims,
      k = k,
      metrics = metrics,
      neighbor_method = neighbor_method,
      max_cells = max_cells,
      stratify_by = stratify_by,
      rare_fraction = rare_fraction,
      rare_n = rare_n,
      challenge_threshold = challenge_threshold,
      seed = seed,
      notes = unique(notes)
    )
  )
  class(result) <- c("sn_integration_assessment", "list")
  result
}

#' Sweep clustering resolutions and summarize cluster-quality diagnostics
#'
#' This helper reruns graph clustering across a grid of Seurat resolution
#' values and evaluates each partition with one or more quality diagnostics.
#' It is intended to support empirical resolution selection rather than relying
#' on a fixed default cluster_by count.
#'
#' @param x A Seurat object containing the reduction used for clustering.
#' @param resolutions Numeric vector of Seurat resolution values to evaluate.
#' @param reduction Reduction name used for neighbor construction and optional
#'   silhouette scoring. Defaults to \code{"harmony"} when present, otherwise
#'   \code{"pca"}.
#' @param dims Optional integer vector of embedding dimensions used by
#'   \code{FindNeighbors()}, \code{FindClusters()}, and silhouette scoring.
#' @param cluster_name Metadata column used to store the temporary cluster
#'   assignments. Defaults to \code{"seurat_clusters"}.
#' @param metrics Metrics to calculate for each resolution. Supported values are
#'   \code{"n_clusters"}, \code{"mean_silhouette"},
#'   \code{"graph_connectivity"}, \code{"cluster_purity"},
#'   \code{"clustering_agreement"}, and \code{"rogue"}.
#' @param label_by Optional metadata column containing reference labels. Required
#'   for \code{"cluster_purity"} and \code{"clustering_agreement"}.
#' @param max_cells Optional integer cap used by expensive embedding-based
#'   metrics such as silhouette.
#' @param rogue_max_cells Optional integer cap used by \code{"rogue"}.
#' @param assay Assay used when \code{"rogue"} is requested.
#' @param layer Layer used when \code{"rogue"} is requested.
#' @param neighbor_method Strategy used when graph connectivity must rebuild a
#'   graph from embeddings.
#' @param k Number of neighbors used for graph connectivity when rebuilding a
#'   graph from embeddings.
#' @param seed Random seed forwarded to Seurat clustering and subsampling.
#' @param n_trees Number of Annoy trees used by graph connectivity when
#'   \code{neighbor_method = "annoy"}.
#'
#' @return A list with a per-resolution summary table and the recommended
#'   resolution based on the mean of the available scaled quality metrics.
#'
#' @examples
#' \dontrun{
#' sweep <- sn_sweep_cluster_resolution(
#'   pbmc,
#'   resolutions = seq(0.2, 1.0, by = 0.2),
#'   reduction = "pca",
#'   dims = 1:20
#' )
#' sweep$summary
#' sweep$recommended_resolution
#' }
#'
#' @export
sn_sweep_cluster_resolution <- function(
  x,
  resolutions = seq(0.2, 1.2, by = 0.2),
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  cluster_name = "seurat_clusters",
  metrics = c("n_clusters", "mean_silhouette", "graph_connectivity"),
  label_by = NULL,
  max_cells = 3000,
  rogue_max_cells = max_cells,
  assay = "RNA",
  layer = "counts",
  neighbor_method = c("auto", "graph", "annoy", "exact"),
  k = 20,
  seed = 717,
  n_trees = 50
) {
  check_installed("Seurat")
  if (!inherits(x, "Seurat")) {
    stop("Input `x` must be a Seurat object.")
  }
  neighbor_method <- rlang::arg_match(neighbor_method)
  metrics <- unique(match.arg(
    metrics,
    c("n_clusters", "mean_silhouette", "graph_connectivity", "cluster_purity", "clustering_agreement", "rogue"),
    several.ok = TRUE
  ))
  resolutions <- sort(unique(as.numeric(resolutions)))
  resolutions <- resolutions[is.finite(resolutions) & resolutions > 0]
  if (length(resolutions) == 0) {
    stop("`resolutions` must contain at least one positive numeric value.")
  }
  if ((any(metrics %in% c("cluster_purity", "clustering_agreement"))) && is.null(label_by)) {
    stop("`label_by` must be supplied when requesting `cluster_purity` or `clustering_agreement`.")
  }

  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = NULL,
    max_cells = NULL,
    stratify_by = NULL,
    seed = seed,
    required_cols = label_by
  )
  dims_use <- dims %||% seq_len(ncol(metric_input$embeddings))

  base_object <- x
  base_object <- Seurat::FindNeighbors(
    object = base_object,
    reduction = reduction,
    dims = dims_use,
    verbose = FALSE
  )

  result_rows <- lapply(resolutions, function(current_resolution) {
    current_object <- Seurat::FindClusters(
      object = base_object,
      resolution = current_resolution,
      random.seed = seed,
      verbose = FALSE
    )
    clusters <- as.character(current_object[[cluster_name, drop = TRUE]])
    row <- data.frame(
      resolution = current_resolution,
      n_clusters = length(unique(clusters)),
      stringsAsFactors = FALSE
    )

    if ("mean_silhouette" %in% metrics) {
      silhouette_tbl <- sn_calculate_silhouette(
        current_object,
        label_by = cluster_name,
        reduction = reduction,
        dims = dims_use,
        max_cells = max_cells,
        stratify_by = cluster_name,
        seed = seed
      )
      row$mean_silhouette <- mean(silhouette_tbl$silhouette_width, na.rm = TRUE)
      row$scaled_silhouette <- .sn_scale_silhouette(row$mean_silhouette)
    }

    if ("graph_connectivity" %in% metrics) {
      connectivity_tbl <- sn_calculate_graph_connectivity(
        current_object,
        label_by = cluster_name,
        reduction = reduction,
        dims = dims_use,
        k = k,
        neighbor_method = neighbor_method,
        max_cells = max_cells,
        stratify_by = cluster_name,
        seed = seed,
        n_trees = n_trees
      )
      row$graph_connectivity <- mean(connectivity_tbl$connectivity_score, na.rm = TRUE)
      row$scaled_graph_connectivity <- row$graph_connectivity
    }

    if ("cluster_purity" %in% metrics) {
      purity_tbl <- sn_calculate_cluster_purity(current_object, cluster_by = cluster_name, label_by = label_by)
      row$cluster_purity <- mean(purity_tbl$purity_score, na.rm = TRUE)
      row$scaled_cluster_purity <- row$cluster_purity
    }

    if ("clustering_agreement" %in% metrics) {
      agreement_tbl <- sn_calculate_clustering_agreement(current_object, cluster_by = cluster_name, label_by = label_by)
      row$ari <- agreement_tbl$ari
      row$nmi <- agreement_tbl$nmi
      row$scaled_ari <- .sn_scale_silhouette(row$ari)
      row$scaled_nmi <- min(max(row$nmi, 0), 1)
    }

    if ("rogue" %in% metrics) {
      rogue_tbl <- sn_calculate_rogue(
        current_object,
        cluster_by = cluster_name,
        assay = assay,
        layer = layer,
        max_cells = rogue_max_cells,
        stratify_by = cluster_name,
        seed = seed
      )
      row$mean_rogue <- mean(rogue_tbl$rogue, na.rm = TRUE)
      row$scaled_rogue <- min(max(row$mean_rogue, 0), 1)
    }

    score_cols <- intersect(
      c(
        "scaled_silhouette",
        "scaled_graph_connectivity",
        "scaled_cluster_purity",
        "scaled_ari",
        "scaled_nmi",
        "scaled_rogue"
      ),
      colnames(row)
    )
    row$composite_score <- if (length(score_cols) > 0) mean(unlist(row[score_cols]), na.rm = TRUE) else NA_real_
    row
  })

  summary_tbl <- .sn_bind_rows(result_rows)
  best_idx <- which(summary_tbl$composite_score == max(summary_tbl$composite_score, na.rm = TRUE))
  best_idx <- best_idx[[1]]

  list(
    summary = summary_tbl,
    recommended_resolution = summary_tbl$resolution[[best_idx]],
    recommended_n_clusters = summary_tbl$n_clusters[[best_idx]]
  ) |>
    structure(class = c("sn_cluster_resolution_sweep", "list"))
}

.sn_validate_constant_within_sample <- function(metadata, sample_col, group_col) {
  summary_tbl <- metadata |>
    dplyr::filter(!is.na(.data[[sample_col]])) |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      n_group = length(unique(.data[[group_col]][!is.na(.data[[group_col]])])),
      .groups = "drop"
    )

  if (any(summary_tbl$n_group > 1L)) {
    stop(
      "Column `", group_col, "` is not constant within samples defined by `",
      sample_col, "`. Composition comparisons require one group label per sample.",
      call. = FALSE
    )
  }
}

.sn_prepare_sample_design <- function(metadata,
                                      sample_col,
                                      group_col,
                                      covariates = NULL,
                                      contrast = NULL) {
  design_cols <- unique(c(sample_col, group_col, covariates))
  missing_cols <- setdiff(design_cols, colnames(metadata))
  if (length(missing_cols) > 0L) {
    stop("Missing required design columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  for (col in unique(c(group_col, covariates))) {
    .sn_validate_constant_within_sample(
      metadata = metadata,
      sample_col = sample_col,
      group_col = col
    )
  }

  design_df <- metadata |>
    dplyr::filter(!is.na(.data[[sample_col]]), !is.na(.data[[group_col]])) |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(unique(c(group_col, covariates))), dplyr::first),
      .groups = "drop"
    ) |>
    as.data.frame(stringsAsFactors = FALSE)

  if (!is.null(contrast)) {
    design_df <- design_df[design_df[[group_col]] %in% contrast, , drop = FALSE]
    if (length(unique(as.character(design_df[[group_col]]))) != 2L) {
      stop("The requested `contrast` must retain exactly two groups in the design data.", call. = FALSE)
    }
    design_df[[group_col]] <- factor(
      as.character(design_df[[group_col]]),
      levels = c(contrast[[2]], contrast[[1]])
    )
  } else {
    group_levels <- unique(as.character(design_df[[group_col]]))
    group_levels <- group_levels[!is.na(group_levels)]
    if (length(group_levels) != 2L) {
      stop(
        "`sn_run_milo()` currently supports two-group comparisons. Supply `contrast = c(case, control)` or subset your data to two groups.",
        call. = FALSE
      )
    }
    design_df[[group_col]] <- factor(as.character(design_df[[group_col]]), levels = group_levels)
  }

  rownames(design_df) <- design_df[[sample_col]]
  design_df
}

#' Run neighborhood differential abundance testing with miloR
#'
#' This wrapper runs the standard miloR neighborhood differential abundance
#' workflow on a Seurat object using one embedding and one sample-level group
#' contrast. Cells are grouped into neighborhoods, counted per sample, and then
#' tested for differential abundance between the two sample groups.
#'
#' @param x A Seurat object.
#' @param sample_by Metadata column defining biological samples.
#' @param group_by Metadata column defining the sample-level comparison group.
#' @param contrast Optional character vector of length 2 giving the comparison
#'   as \code{c(case, control)}. When omitted, \code{group_by} must contain
#'   exactly two levels in the selected cells.
#' @param reduction Reduction name used to build neighborhoods. Defaults to
#'   \code{"pca"}.
#' @param dims Optional integer vector of embedding dimensions to retain.
#' @param cells Optional character vector of cells to include.
#' @param max_cells Optional integer cap used to subsample cells before
#'   neighborhood construction.
#' @param stratify_by Optional metadata column used to preserve representation
#'   during subsampling. Defaults to \code{sample_by}.
#' @param k Number of neighbors for graph and neighborhood construction.
#' @param d Number of embedding dimensions passed to miloR. Defaults to the
#'   selected embedding dimensionality.
#' @param prop Proportion of cells sampled as neighborhood indices.
#' @param refined Logical; if \code{TRUE}, use miloR's refined neighborhood
#'   sampling.
#' @param refinement_scheme Refinement scheme passed to
#'   \code{miloR::makeNhoods()}.
#' @param covariates Optional sample-level covariates added to the DA design
#'   formula alongside \code{group_by}.
#' @param annotation_by Optional cell-level metadata column used to annotate
#'   neighborhoods with \code{miloR::annotateNhoods()}.
#' @param fdr_weighting FDR weighting strategy passed to
#'   \code{miloR::testNhoods()}.
#' @param min_mean Minimum mean count threshold passed to
#'   \code{miloR::testNhoods()}.
#' @param norm_method Normalization method passed to
#'   \code{miloR::testNhoods()}.
#' @param result_id Stable identifier for the stored Milo result.
#' @param result_id Optional explicit result identifier. Overrides
#'   \code{result_id} when supplied.
#'   When supplied, the milo result is stored on the Seurat object.
#' @param return_object Logical; when \code{TRUE} and \code{result_id} is
#'   supplied, return the updated Seurat object.
#' @param return_intermediate Logical; if \code{TRUE}, return a list with the
#'   DA table, design data, and milo object.
#' @param verbose Logical; if \code{TRUE}, emit progress logs.
#' @param seed Optional random seed recorded in the stored result provenance.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return By default, a data frame of neighborhood-level DA statistics. When
#'   \code{result_id} is supplied, return the unified stored-result list, or
#'   the updated Seurat object when \code{return_object = TRUE}. When
#'   \code{return_intermediate = TRUE}, return a list with \code{table},
#'   \code{design_df}, and \code{milo}.
#'
#' @examples
#' \dontrun{
#' da_tbl <- sn_run_milo(
#'   seu,
#'   sample_by = "sample",
#'   group_by = "condition",
#'   contrast = c("treated", "control"),
#'   reduction = "pca"
#' )
#' }
#'
#' @export
sn_run_milo <- function(x,
                        sample_by = NULL,
                        group_by = NULL,
                        contrast = NULL,
                        reduction = "pca",
                        dims = NULL,
                        cells = NULL,
                        max_cells = NULL,
                        stratify_by = NULL,
                        k = 20,
                        d = NULL,
                        prop = 0.1,
                        refined = TRUE,
                        refinement_scheme = "reduced_dim",
                        covariates = NULL,
                        annotation_by = NULL,
                        fdr_weighting = c("k-distance", "neighbour-distance", "max", "graph-overlap", "none"),
                        min_mean = 0,
                        norm_method = c("TMM", "RLE", "logMS"),
                        result_id = NULL,
                        return_object = FALSE,
                        return_intermediate = FALSE,
                        verbose = TRUE,
                        seed = NULL,
                        object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  if (!is_null(result_id)) {
    result_id <- .sn_validate_result_id(result_id)
  }
  check_installed(c("Seurat", "SingleCellExperiment", "miloR"))
  stopifnot(is.character(sample_by), length(sample_by) == 1L)
  stopifnot(is.character(group_by), length(group_by) == 1L)
  stopifnot(is.null(contrast) || (is.character(contrast) && length(contrast) == 2L))
  stopifnot(is.null(covariates) || is.character(covariates))
  stopifnot(is.null(annotation_by) || (is.character(annotation_by) && length(annotation_by) == 1L))
  stopifnot(is.null(result_id) || (is.character(result_id) && length(result_id) == 1L))
  stopifnot(is.logical(return_object), length(return_object) == 1L)
  stopifnot(is.logical(return_intermediate), length(return_intermediate) == 1L)
  stopifnot(is.logical(verbose), length(verbose) == 1L)
  if (isTRUE(return_intermediate) && !is.null(result_id)) {
    stop("`return_intermediate = TRUE` cannot be combined with `result_id`.", call. = FALSE)
  }

  fdr_weighting <- match.arg(fdr_weighting)
  norm_method <- match.arg(norm_method)

  required_cols <- unique(c(sample_by, group_by, covariates, annotation_by))
  metric_input <- .sn_prepare_metric_input(
    x = x,
    reduction = reduction,
    dims = dims,
    cells = cells,
    max_cells = max_cells,
    stratify_by = stratify_by %||% sample_by,
    required_cols = required_cols
  )

  metadata <- metric_input$metadata
  if (!is.null(contrast)) {
    metadata <- metadata[metadata[[group_by]] %in% contrast, , drop = FALSE]
    metric_input$embeddings <- metric_input$embeddings[rownames(metadata), , drop = FALSE]
    metric_input$cells <- rownames(metadata)
  }

  if (nrow(metadata) < 2L) {
    stop("Too few cells remain after filtering to run miloR.", call. = FALSE)
  }

  design_df <- .sn_prepare_sample_design(
    metadata = metadata,
    sample_col = sample_by,
    group_col = group_by,
    covariates = covariates,
    contrast = contrast
  )

  if (nrow(design_df) < 2L) {
    stop("At least two samples are required to run neighborhood differential abundance testing.", call. = FALSE)
  }

  if (is.null(d)) {
    d <- ncol(metric_input$embeddings)
  }
  d <- min(as.integer(d), ncol(metric_input$embeddings))
  if (d < 1L) {
    stop("`d` must be at least 1.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    .sn_log_info("Converting the selected cells to a SingleCellExperiment for miloR.")
  }
  seu_subset <- x[, metric_input$cells]
  sce <- Seurat::as.SingleCellExperiment(seu_subset)
  SingleCellExperiment::reducedDim(sce, "shennong_milo") <- as.matrix(metric_input$embeddings)

  milo_obj <- miloR::Milo(sce)

  if (isTRUE(verbose)) {
    .sn_log_info("Building the miloR neighborhood graph with k = {k} and d = {d}.")
  }
  milo_obj <- miloR::buildGraph(
    milo_obj,
    k = k,
    d = d,
    reduced.dim = "shennong_milo"
  )

  if (isTRUE(verbose)) {
    .sn_log_info("Sampling and refining neighborhoods with prop = {prop}.")
  }
  milo_obj <- miloR::makeNhoods(
    milo_obj,
    prop = prop,
    k = k,
    d = d,
    refined = refined,
    reduced_dims = "shennong_milo",
    refinement_scheme = refinement_scheme
  )

  if (isTRUE(verbose)) {
    .sn_log_info("Counting cells per sample across miloR neighborhoods.")
  }
  milo_obj <- miloR::countCells(
    milo_obj,
    samples = sample_by,
    meta.data = as.data.frame(SummarizedExperiment::colData(milo_obj))
  )

  design_terms <- c(covariates, group_by)
  design_formula <- stats::reformulate(termlabels = design_terms)

  if (isTRUE(verbose)) {
    .sn_log_info("Testing neighborhoods for differential abundance across `{group_by}`.")
  }
  da_table <- miloR::testNhoods(
    milo_obj,
    design = design_formula,
    design.df = design_df,
    fdr.weighting = fdr_weighting,
    min.mean = min_mean,
    norm.method = norm_method
  )

  if (!is.null(annotation_by)) {
    if (isTRUE(verbose)) {
      .sn_log_info("Annotating neighborhoods using `{annotation_by}`.")
    }
    da_table <- miloR::annotateNhoods(
      milo_obj,
      da.res = da_table,
      coldata_col = annotation_by
    )
  }

  fdr_rank <- if ("SpatialFDR" %in% colnames(da_table)) da_table$SpatialFDR else da_table$FDR
  da_table <- da_table[order(fdr_rank, da_table$PValue), , drop = FALSE]
  da_table$comparison <- paste(levels(design_df[[group_by]])[[2]], "vs", levels(design_df[[group_by]])[[1]])
  da_table$sample_col <- sample_by
  da_table$group_col <- group_by
  da_table$reduction <- reduction

  if (!is.null(result_id)) {
    stored_result <- sn_store_milo(
      object = x,
      result = da_table,
      result_id = result_id,
      sample_by = sample_by,
      group_by = group_by,
      comparison = da_table$comparison[[1]],
      reduction = reduction,
      dims = colnames(metric_input$embeddings),
      annotation_by = annotation_by,
      random_seed = seed,
      return_object = FALSE
    )

    if (isTRUE(return_object)) {
      object <- sn_store_result(
        object = x,
        type = "milo",
        result_id = result_id,
        result = stored_result
      )
      return(.sn_log_seurat_command(object = object, name = "sn_run_milo"))
    }

    return(stored_result)
  }

  if (isTRUE(return_intermediate)) {
    return(list(
      table = da_table,
      design_df = design_df,
      milo = milo_obj
    ))
  }

  da_table
}

#' Store a miloR differential-abundance result on a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param result A neighborhood-level differential-abundance table.
#' @param result_id Stable identifier for the stored Milo result.
#' @param result_id Optional explicit result identifier. Overrides
#'   \code{result_id} when supplied.
#' @param sample_by Sample column used for the design.
#' @param group_by Group column used for the design.
#' @param comparison Human-readable comparison label.
#' @param reduction Reduction used to build the milo neighborhoods.
#' @param dims Optional embedding dimension names or indices used for milo.
#' @param annotation_by Optional neighborhood annotation column.
#' @param random_seed Optional random seed recorded in the result provenance.
#' @param return_object If \code{TRUE}, return the updated object.
#'
#' @return A \code{Seurat} object or stored-result list.
#' @export
sn_store_milo <- function(object,
                          result,
                          result_id = "default",
                          sample_by = NULL,
                          group_by = NULL,
                          comparison = NULL,
                          reduction = "pca",
                          dims = NULL,
                          annotation_by = NULL,
                          random_seed = NULL,
                          return_object = TRUE) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)
  stopifnot(is.character(sample_by), length(sample_by) == 1L)
  stopifnot(is.character(group_by), length(group_by) == 1L)

  stored_result <- list(
    schema_version = .sn_analysis_result_schema_version(),
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = tibble::as_tibble(result),
    analysis = "milo",
    method = "miloR",
    sample_col = sample_by,
    group_col = group_by,
    comparison = comparison,
    reduction = reduction,
    dims = dims,
    annotation_by = annotation_by,
    provenance = .sn_analysis_provenance(random_seed = random_seed %||% NA_integer_)
  )

  object <- sn_store_result(
    object = object,
    type = "milo",
    result_id = result_id,
    result = stored_result
  )

  if (isTRUE(return_object)) {
    return(.sn_log_seurat_command(object = object, name = "sn_store_milo"))
  }

  sn_get_result(
    object = object,
    type = "milo",
    result_id = result_id
  )
}

#' Retrieve a stored miloR result from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param result_id Name of the stored milo result.
#' @param annotation Optional subset of annotation labels to keep.
#' @param spatial_fdr Optional maximum \code{SpatialFDR} threshold.
#' @param with_metadata If \code{TRUE}, return the full stored-result list.
#'
#' @return A tibble or stored-result list.
#' @export
sn_get_milo_result <- function(object,
                               result_id = "default",
                               annotation = NULL,
                               spatial_fdr = NULL,
                               with_metadata = FALSE) {
  .sn_validate_seurat_object(object)

  stored <- sn_get_result(
    object = object,
    type = "milo",
    result_id = result_id
  )
  if (isTRUE(with_metadata)) {
    return(stored)
  }

  table <- tibble::as_tibble(stored$tables$primary)
  annotation_by <- stored$annotation_col %||% NULL
  if (!is.null(annotation) && !is.null(annotation_by) && annotation_by %in% colnames(table)) {
    table <- dplyr::filter(table, .data[[annotation_by]] %in% annotation)
  }
  if (!is.null(spatial_fdr) && "SpatialFDR" %in% colnames(table)) {
    table <- dplyr::filter(table, .data$SpatialFDR <= spatial_fdr)
  }
  table
}

.sn_qc_assessment_sample_col <- function(object, sample_col = NULL) {
  metadata <- object@meta.data
  if (!is_null(sample_col)) {
    if (!sample_col %in% colnames(metadata)) {
      stop("Sample column '", sample_col, "' was not found in the Seurat metadata.", call. = FALSE)
    }
    return(sample_col)
  }

  candidates <- c("sample", "orig.ident")
  available <- candidates[candidates %in% colnames(metadata)]
  if (length(available) > 0) {
    return(available[[1]])
  }

  NULL
}

.sn_qc_assessment_count_col <- function(metadata) {
  candidates <- c("nCount_RNA_corrected", "nCount_RNA")
  available <- candidates[candidates %in% colnames(metadata)]
  if (length(available) > 0) available[[1]] else NULL
}

.sn_qc_assessment_feature_col <- function(metadata) {
  candidates <- c("nFeature_RNA_corrected", "nFeature_RNA")
  available <- candidates[candidates %in% colnames(metadata)]
  if (length(available) > 0) available[[1]] else NULL
}

.sn_qc_assessment_doublet_col <- function(metadata) {
  candidates <- c("scDblFinder.class_corrected", "scDblFinder.class")
  available <- candidates[candidates %in% colnames(metadata)]
  if (length(available) > 0) available[[1]] else NULL
}

.sn_qc_assessment_zero_col <- function(metadata) {
  zero_cols <- grep("_zero_count$", colnames(metadata), value = TRUE)
  if (length(zero_cols) == 0) {
    return(NULL)
  }
  corrected_first <- grep("^decontaminated_counts_zero_count$", zero_cols, value = TRUE)
  if (length(corrected_first) > 0) {
    return(corrected_first[[1]])
  }
  zero_cols[[1]]
}

.sn_qc_assessment_flag_cols <- function(metadata) {
  grep("_qc$", colnames(metadata), value = TRUE)
}

.sn_qc_quality_label <- function(score) {
  dplyr::case_when(
    is.na(score) ~ "unknown",
    score >= 85 ~ "high",
    score >= 70 ~ "good",
    score >= 50 ~ "mixed",
    TRUE ~ "poor"
  )
}

.sn_qc_format_percent <- function(x, digits = 1) {
  if (is.na(x)) {
    return("NA")
  }
  paste0(formatC(100 * x, format = "f", digits = digits), "%")
}

.sn_qc_scale_score <- function(x, good, minimum = 0) {
  if (is.na(x)) {
    return(NA_real_)
  }
  100 * max(min((x - minimum) / max(good - minimum, .Machine$double.eps), 1), 0)
}

.sn_qc_inverse_score <- function(x, warning, failure) {
  if (is.na(x)) {
    return(NA_real_)
  }
  if (x <= warning) {
    return(100)
  }
  if (x >= failure) {
    return(0)
  }
  100 * (1 - (x - warning) / (failure - warning))
}

.sn_qc_current_summary <- function(metadata, sample_col = NULL) {
  metadata <- as.data.frame(metadata)
  count_col <- .sn_qc_assessment_count_col(metadata)
  feature_col <- .sn_qc_assessment_feature_col(metadata)
  doublet_col <- .sn_qc_assessment_doublet_col(metadata)
  zero_col <- .sn_qc_assessment_zero_col(metadata)
  qc_cols <- .sn_qc_assessment_flag_cols(metadata)

  if (is_null(sample_col)) {
    metadata$.sn_qc_sample <- "all_cells"
    sample_col <- ".sn_qc_sample"
  }

  split_meta <- split(metadata, metadata[[sample_col]], drop = TRUE)
  summary_tbl <- lapply(names(split_meta), function(sample_name) {
    current_meta <- split_meta[[sample_name]]
    failed_qc_fraction <- if (length(qc_cols) > 0) {
      mean(rowSums(current_meta[, qc_cols, drop = FALSE] == "Failed", na.rm = TRUE) > 0)
    } else {
      NA_real_
    }
    doublet_fraction <- if (!is_null(doublet_col)) {
      mean(tolower(as.character(current_meta[[doublet_col]])) == "doublet", na.rm = TRUE)
    } else {
      NA_real_
    }
    zero_fraction <- if (!is_null(zero_col)) {
      mean(as.logical(current_meta[[zero_col]]), na.rm = TRUE)
    } else {
      NA_real_
    }

    library_score <- mean(c(
      if (!is_null(feature_col)) .sn_qc_scale_score(stats::median(current_meta[[feature_col]], na.rm = TRUE), good = 1500, minimum = 200) else NA_real_,
      if (!is_null(count_col)) .sn_qc_scale_score(stats::median(current_meta[[count_col]], na.rm = TRUE), good = 5000, minimum = 500) else NA_real_
    ), na.rm = TRUE)
    if (is.nan(library_score)) {
      library_score <- NA_real_
    }

    stress_score <- mean(c(
      if ("percent.mt" %in% colnames(current_meta)) .sn_qc_inverse_score(stats::median(current_meta$percent.mt, na.rm = TRUE), warning = 10, failure = 25) else NA_real_,
      if (!is.na(failed_qc_fraction)) .sn_qc_inverse_score(failed_qc_fraction, warning = 0.1, failure = 0.4) else NA_real_,
      if (!is.na(doublet_fraction)) .sn_qc_inverse_score(doublet_fraction, warning = 0.08, failure = 0.2) else NA_real_,
      if (!is.na(zero_fraction)) .sn_qc_inverse_score(zero_fraction, warning = 0.01, failure = 0.1) else NA_real_
    ), na.rm = TRUE)
    if (is.nan(stress_score)) {
      stress_score <- NA_real_
    }

    score <- mean(c(library_score, stress_score), na.rm = TRUE)
    if (is.nan(score)) {
      score <- NA_real_
    }

    data.frame(
      sample = as.character(sample_name),
      n_cells = nrow(current_meta),
      median_nCount = if (!is_null(count_col)) stats::median(current_meta[[count_col]], na.rm = TRUE) else NA_real_,
      median_nFeature = if (!is_null(feature_col)) stats::median(current_meta[[feature_col]], na.rm = TRUE) else NA_real_,
      median_percent_mt = if ("percent.mt" %in% colnames(current_meta)) stats::median(current_meta$percent.mt, na.rm = TRUE) else NA_real_,
      failed_qc_fraction = failed_qc_fraction,
      doublet_fraction = doublet_fraction,
      zero_count_fraction = zero_fraction,
      library_score = library_score,
      stress_score = stress_score,
      qc_score = score,
      qc_label = .sn_qc_quality_label(score),
      stringsAsFactors = FALSE
    )
  })

  .sn_bind_rows(summary_tbl)
}

.sn_qc_compare_to_reference <- function(object, reference, sample_col = NULL) {
  current_meta <- object@meta.data
  reference_meta <- reference@meta.data
  current_cells <- rownames(current_meta)
  qc_cols <- .sn_qc_assessment_flag_cols(reference_meta)
  doublet_col <- .sn_qc_assessment_doublet_col(reference_meta)

  if (is_null(sample_col)) {
    reference_meta$.sn_qc_sample <- "all_cells"
    current_meta$.sn_qc_sample <- "all_cells"
    sample_col <- ".sn_qc_sample"
  }

  split_reference <- split(reference_meta, reference_meta[[sample_col]], drop = TRUE)
  compare_tbl <- lapply(names(split_reference), function(sample_name) {
    ref_meta <- split_reference[[sample_name]]
    keep <- rownames(ref_meta) %in% current_cells

    low_quality <- if (length(qc_cols) > 0) {
      rowSums(ref_meta[, qc_cols, drop = FALSE] == "Failed", na.rm = TRUE) > 0
    } else {
      rep(FALSE, nrow(ref_meta))
    }
    doublet <- if (!is_null(doublet_col)) {
      tolower(as.character(ref_meta[[doublet_col]])) == "doublet"
    } else {
      rep(FALSE, nrow(ref_meta))
    }
    clean <- !low_quality & !doublet

    data.frame(
      sample = as.character(sample_name),
      n_cells_reference = nrow(ref_meta),
      retention_fraction = mean(keep),
      low_quality_removed_fraction = if (any(low_quality)) mean(!keep[low_quality]) else NA_real_,
      doublet_removed_fraction = if (any(doublet)) mean(!keep[doublet]) else NA_real_,
      clean_retained_fraction = if (any(clean)) mean(keep[clean]) else NA_real_,
      comparison_score = mean(c(
        if (any(low_quality)) 100 * mean(!keep[low_quality]) else NA_real_,
        if (any(doublet)) 100 * mean(!keep[doublet]) else NA_real_,
        if (any(clean)) 100 * mean(keep[clean]) else NA_real_
      ), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  .sn_bind_rows(compare_tbl)
}

.sn_qc_assessment_messages <- function(overall, by_sample, has_reference = FALSE) {
  lines <- c(
    glue::glue(
      "QC assessment: {overall$qc_label[[1]]} quality (score {round(overall$qc_score[[1]], 1)}/100) across {overall$n_samples[[1]]} sample(s) and {overall$n_cells[[1]]} cell(s)."
    )
  )

  if (has_reference) {
    lines <- c(
      lines,
      glue::glue(
        "QC comparison: retained {.sn_qc_format_percent(overall$retention_fraction[[1]])} of reference cells, removed {.sn_qc_format_percent(overall$low_quality_removed_fraction[[1]])} of flagged low-quality cells, and removed {.sn_qc_format_percent(overall$doublet_removed_fraction[[1]])} of called doublets."
      )
    )
  }

  sample_lines <- apply(by_sample, 1, function(row) {
    glue::glue(
      "Sample {row[['sample']]}: {row[['qc_label']]} quality, {row[['n_cells']]} cells, median features {round(as.numeric(row[['median_nFeature']]), 1)}, median counts {round(as.numeric(row[['median_nCount']]), 1)}, median percent.mt {round(as.numeric(row[['median_percent_mt']]), 1)}."
    )
  })

  c(lines, sample_lines)
}

#' Assess overall QC status and before/after filtering outcomes
#'
#' This function summarizes the QC status of a Seurat object, optionally
#' compares it against a reference object captured before filtering, and
#' produces per-sample and overall quality scores. It is designed to answer
#' whether low-quality cells and called doublets were removed while retaining
#' clean cells, and whether ambient-RNA correction introduced problematic
#' zero-count cells.
#'
#' @param object A Seurat object to assess.
#' @param reference Optional Seurat object representing the pre-filter state.
#'   When supplied, the function compares the retained cells in \code{object}
#'   against the QC flags and doublet calls recorded in \code{reference}.
#' @param sample_by Optional metadata column defining samples. When
#'   \code{NULL}, the function uses \code{sample} or \code{orig.ident} when
#'   available and otherwise treats the object as one sample.
#' @param result_id Name used when storing the assessment under
#'   \code{object@misc$qc_assessments}.
#' @param result_id Optional explicit result identifier. Overrides
#'   \code{result_id} when supplied.
#' @param return_object Logical; when \code{TRUE}, store the assessment in the
#'   Seurat object and return the updated object.
#' @param verbose Logical; when \code{TRUE}, print a concise QC summary.
#'
#' @return A list with \code{overall}, \code{by_sample}, \code{comparison},
#'   and \code{messages} when \code{return_object = FALSE}; otherwise the
#'   updated Seurat object with the stored assessment.
#'
#' @examples
#' \dontrun{
#'   pbmc <- qs2::qs_read(file.path(
#'     Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
#'   ))
#'   qc_report <- sn_assess_qc(pbmc, verbose = FALSE)
#'   qc_report$tables$overall
#' }
#'
#' @export
sn_assess_qc <- function(object,
                         reference = NULL,
                         sample_by = NULL,
                         result_id = "default",
                         return_object = FALSE,
                         verbose = TRUE) {
  result_id <- .sn_validate_result_id(result_id)
  check_installed("SeuratObject")
  .sn_validate_seurat_object(object)
  if (!is_null(reference) && !inherits(reference, "Seurat")) {
    stop("`reference` must be NULL or a Seurat object.", call. = FALSE)
  }

  sample_by <- .sn_qc_assessment_sample_col(object, sample_col = sample_by)
  by_sample <- .sn_qc_current_summary(object@meta.data, sample_col = sample_by)
  comparison <- if (!is_null(reference)) {
    reference_sample_col <- .sn_qc_assessment_sample_col(reference, sample_col = sample_by)
    .sn_qc_compare_to_reference(object = object, reference = reference, sample_col = reference_sample_col)
  } else {
    NULL
  }

  if (!is_null(comparison) && nrow(comparison) > 0) {
    by_sample <- dplyr::left_join(by_sample, comparison, by = "sample")
  }

  overall_score <- mean(c(by_sample$qc_score, by_sample$comparison_score), na.rm = TRUE)
  if (is.nan(overall_score)) {
    overall_score <- mean(by_sample$qc_score, na.rm = TRUE)
  }
  overall <- data.frame(
    n_samples = nrow(by_sample),
    n_cells = sum(by_sample$n_cells),
    qc_score = overall_score,
    qc_label = .sn_qc_quality_label(overall_score),
    retention_fraction = if ("retention_fraction" %in% colnames(by_sample)) stats::weighted.mean(by_sample$retention_fraction, by_sample$n_cells_reference, na.rm = TRUE) else NA_real_,
    low_quality_removed_fraction = if ("low_quality_removed_fraction" %in% colnames(by_sample)) stats::weighted.mean(by_sample$low_quality_removed_fraction, by_sample$n_cells_reference, na.rm = TRUE) else NA_real_,
    doublet_removed_fraction = if ("doublet_removed_fraction" %in% colnames(by_sample)) stats::weighted.mean(by_sample$doublet_removed_fraction, by_sample$n_cells_reference, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )

  report <- list(
    schema_version = .sn_analysis_result_schema_version(),
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    overall = overall,
    by_sample = by_sample,
    comparison = comparison,
    sample_col = sample_by %||% NA_character_,
    reference_provided = !is_null(reference)
  )
  report$messages <- .sn_qc_assessment_messages(
    overall = overall,
    by_sample = by_sample,
    has_reference = !is_null(reference)
  )

  if (isTRUE(verbose)) {
    for (line in report$messages) {
      .sn_log_info("{line}")
    }
  }

  if (isTRUE(return_object)) {
    object <- sn_store_result(
      object = object,
      type = "qc_assessment",
      result_id = result_id,
      result = report
    )
    return(object)
  }

  .sn_prepare_result(
    type = "qc_assessment",
    result_id = result_id,
    result = report
  )
}

.sn_default_metric_reduction <- function(x) {
  reductions <- names(x@reductions)
  candidates <- c("harmony", "integrated", "integrated.rpca", "pca", "umap")
  available <- candidates[candidates %in% reductions]
  if (length(available) > 0) {
    return(available[[1]])
  }
  if (length(reductions) == 0) {
    stop("No dimensional reductions were found in the Seurat object.")
  }
  reductions[[1]]
}

.sn_default_baseline_reduction <- function(x, reduction) {
  reductions <- names(x@reductions)
  if (!identical(reduction, "pca") && "pca" %in% reductions) {
    return("pca")
  }
  NULL
}

.sn_check_metric_columns <- function(metadata, columns) {
  columns <- unique(stats::na.omit(columns))
  missing_cols <- setdiff(columns, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

.sn_prepare_metric_input <- function(
  x,
  reduction,
  dims = NULL,
  cells = NULL,
  max_cells = NULL,
  stratify_by = NULL,
  seed = 717,
  required_cols = NULL
) {
  check_installed("SeuratObject")
  if (!inherits(x, "Seurat")) {
    stop("Input `x` must be a Seurat object.")
  }

  metadata <- x@meta.data
  .sn_check_metric_columns(metadata, required_cols)

  reduction_names <- names(x@reductions)
  if (!reduction %in% reduction_names) {
    stop("Reduction '", reduction, "' was not found in the Seurat object.")
  }

  embeddings <- SeuratObject::Embeddings(object = x, reduction = reduction)
  if (!is.null(dims)) {
    dims <- as.integer(dims)
    if (any(dims < 1 | dims > ncol(embeddings))) {
      stop("`dims` must refer to existing embedding dimensions.")
    }
    embeddings <- embeddings[, dims, drop = FALSE]
  }

  available_cells <- rownames(embeddings)
  cells_use <- cells %||% available_cells
  missing_cells <- setdiff(cells_use, available_cells)
  if (length(missing_cells) > 0) {
    stop("Unknown cells requested: ", paste(utils::head(missing_cells, 5), collapse = ", "))
  }

  cells_use <- intersect(available_cells, cells_use)
  metadata <- metadata[cells_use, , drop = FALSE]
  embeddings <- embeddings[cells_use, , drop = FALSE]

  if (!is.null(required_cols) && length(required_cols) > 0) {
    keep <- stats::complete.cases(metadata[, required_cols, drop = FALSE])
    if (!all(keep)) {
      warning(
        "Dropping ", sum(!keep), " cells with missing values in required metric columns.",
        call. = FALSE
      )
      metadata <- metadata[keep, , drop = FALSE]
      embeddings <- embeddings[keep, , drop = FALSE]
      cells_use <- rownames(metadata)
    }
  }

  if (!is.null(max_cells) && nrow(embeddings) > max_cells) {
    cells_use <- .sn_subsample_metric_cells(
      cells = rownames(embeddings),
      metadata = metadata,
      max_cells = max_cells,
      stratify_by = stratify_by,
      seed = seed
    )
    metadata <- metadata[cells_use, , drop = FALSE]
    embeddings <- embeddings[cells_use, , drop = FALSE]
  }

  if (nrow(embeddings) < 2) {
    stop("At least two cells are required to calculate this metric.")
  }

  list(
    embeddings = embeddings,
    metadata = metadata,
    cells = rownames(embeddings)
  )
}

.sn_subsample_metric_cells <- function(cells, metadata, max_cells, stratify_by = NULL, seed = 717) {
  if (length(cells) <= max_cells) {
    return(cells)
  }

  set.seed(seed)
  if (!is.null(stratify_by) && stratify_by %in% colnames(metadata)) {
    strata <- as.character(metadata[[stratify_by]])
    strata[is.na(strata)] <- "__NA__"
    strata_count <- table(strata)
    weights <- 1 / as.numeric(strata_count[strata])
    return(sort(sample(cells, size = max_cells, replace = FALSE, prob = weights)))
  }

  sort(sample(cells, size = max_cells, replace = FALSE))
}

.sn_compute_silhouette_table <- function(embeddings, labels, cell_ids, label_name) {
  labels <- as.factor(labels)
  if (length(unique(labels)) < 2) {
    stop("At least two distinct groups are required to calculate silhouette widths.")
  }
  if (nrow(embeddings) < 3) {
    stop("At least three cells are required to calculate silhouette widths.")
  }

  sil <- cluster::silhouette(
    x = as.integer(labels),
    dist = stats::dist(embeddings)
  )

  data.frame(
    cell_id = cell_ids,
    label = as.character(labels),
    silhouette_width = sil[, "sil_width"],
    stringsAsFactors = FALSE
  ) |>
    stats::setNames(c("cell_id", label_name, "silhouette_width"))
}

.sn_calculate_batch_silhouette_score <- function(
  embeddings,
  batch,
  label = NULL,
  cell_ids,
  batch_name,
  label_name = NULL
) {
  if (is.null(label)) {
    batch_tbl <- .sn_compute_silhouette_table(
      embeddings = embeddings,
      labels = batch,
      cell_ids = cell_ids,
      label_name = batch_name
    )
    raw_score <- mean(abs(batch_tbl$silhouette_width))
    return(list(
      raw_score = raw_score,
      scaled_score = max(0, 1 - raw_score),
      per_cell = batch_tbl,
      per_group = NULL,
      note = "Global batch_by silhouette"
    ))
  }

  split_idx <- split(seq_along(label), as.character(label))
  per_group <- list()
  per_cell <- list()
  notes <- character()

  for (current_group in names(split_idx)) {
    idx <- split_idx[[current_group]]
    batch_levels <- unique(batch[idx])
    if (length(idx) < 3 || length(batch_levels) < 2) {
      notes <- c(
        notes,
        paste0("Skipped batch_by silhouette inside label_by '", current_group, "'.")
      )
      next
    }

    current_tbl <- .sn_compute_silhouette_table(
      embeddings = embeddings[idx, , drop = FALSE],
      labels = batch[idx],
      cell_ids = cell_ids[idx],
      label_name = batch_name
    )
    current_tbl[[label_name]] <- current_group
    per_cell[[length(per_cell) + 1]] <- current_tbl
    per_group[[length(per_group) + 1]] <- data.frame(
      label = current_group,
      n_cells = length(idx),
      mean_abs_batch_silhouette = mean(abs(current_tbl$silhouette_width)),
      scaled_score = max(0, 1 - mean(abs(current_tbl$silhouette_width))),
      stringsAsFactors = FALSE
    )
  }

  per_group <- .sn_bind_rows(per_group)
  if (nrow(per_group) == 0) {
    stop("No label_by groups contained enough cells and batch_by diversity for batch_by silhouette scoring.")
  }

  raw_score <- stats::weighted.mean(
    x = per_group$mean_abs_batch_silhouette,
    w = per_group$n_cells
  )
  scaled_score <- stats::weighted.mean(
    x = per_group$scaled_score,
    w = per_group$n_cells
  )

  list(
    raw_score = raw_score,
    scaled_score = scaled_score,
    per_cell = .sn_bind_rows(per_cell),
    per_group = per_group,
    note = paste(unique(notes), collapse = " ")
  )
}

.sn_guess_graph_name <- function(x, graph = NULL) {
  graph_names <- names(x@graphs)
  if (length(graph_names) == 0) {
    return(NULL)
  }
  if (!is.null(graph)) {
    if (!graph %in% graph_names) {
      stop("Graph '", graph, "' was not found in the Seurat object.")
    }
    return(graph)
  }

  nn_graph <- graph_names[grepl("_nn$", graph_names)]
  if (length(nn_graph) > 0) {
    return(nn_graph[[1]])
  }
  snn_graph <- graph_names[grepl("_snn$", graph_names)]
  if (length(snn_graph) > 0) {
    return(snn_graph[[1]])
  }
  graph_names[[1]]
}

.sn_get_metric_graph <- function(
  x,
  embeddings,
  cells,
  graph = NULL,
  k = 20,
  neighbor_method = "auto",
  n_trees = 50
) {
  inferred_graph <- .sn_guess_graph_name(x, graph)

  if (neighbor_method %in% c("auto", "graph") && !is.null(inferred_graph)) {
    adjacency <- .sn_subset_graph_adjacency(
      graph = x@graphs[[inferred_graph]],
      cells = cells
    )
    return(list(adjacency = adjacency, source = inferred_graph))
  }

  if (identical(neighbor_method, "graph")) {
    stop("No graph was available for `neighbor_method = \"graph\"`.")
  }

  if ((identical(neighbor_method, "auto") || identical(neighbor_method, "annoy")) &&
    rlang::is_installed("Seurat")) {
    adjacency <- .sn_build_annoy_adjacency(
      embeddings = embeddings,
      k = k,
      n_trees = n_trees
    )
    return(list(adjacency = adjacency, source = "annoy"))
  }

  if (identical(neighbor_method, "annoy")) {
    stop("Package 'Seurat' is required for `neighbor_method = \"annoy\"`.")
  }

  adjacency <- .sn_build_exact_adjacency(embeddings = embeddings, k = k)
  list(adjacency = adjacency, source = "exact")
}

.sn_subset_graph_adjacency <- function(graph, cells) {
  adjacency <- methods::as(graph, "dgCMatrix")
  adjacency <- adjacency[cells, cells, drop = FALSE]
  adjacency <- .sn_binarize_adjacency(adjacency)
  adjacency
}

.sn_build_annoy_adjacency <- function(embeddings, k = 20, n_trees = 50) {
  k <- min(as.integer(k), nrow(embeddings) - 1L)
  if (k < 1) {
    stop("At least two cells are required to build a neighbor graph.")
  }

  nn_idx <- .sn_find_annoy_knn(
    embeddings = embeddings,
    k = k,
    n_trees = n_trees,
    include_distance = FALSE
  )
  .sn_neighbor_index_to_adjacency(nn_idx, cell_names = rownames(embeddings))
}

.sn_build_exact_adjacency <- function(embeddings, k = 20) {
  k <- min(as.integer(k), nrow(embeddings) - 1L)
  if (k < 1) {
    stop("At least two cells are required to build a neighbor graph.")
  }

  if (nrow(embeddings) > 5000) {
    warning(
      "Exact kNN search on more than 5000 cells can be slow. ",
      "Prefer an existing graph or `neighbor_method = \"annoy\"`.",
      call. = FALSE
    )
  }

  nn_idx <- .sn_exact_knn(
    embeddings = embeddings,
    k = k,
    include_distance = FALSE
  )$idx
  .sn_neighbor_index_to_adjacency(nn_idx, cell_names = rownames(embeddings))
}

.sn_neighbor_index_to_adjacency <- function(nn_idx, cell_names) {
  if (is.null(dim(nn_idx))) {
    nn_idx <- matrix(nn_idx, ncol = 1)
  }

  row_index <- rep(seq_len(nrow(nn_idx)), times = ncol(nn_idx))
  col_index <- as.vector(nn_idx)
  keep <- !is.na(col_index) & row_index != col_index
  adjacency <- Matrix::sparseMatrix(
    i = row_index[keep],
    j = col_index[keep],
    x = 1,
    dims = c(length(cell_names), length(cell_names)),
    dimnames = list(cell_names, cell_names)
  )
  .sn_binarize_adjacency(adjacency)
}

.sn_binarize_adjacency <- function(adjacency) {
  adjacency <- methods::as(adjacency, "dgCMatrix")
  adjacency <- adjacency + Matrix::t(adjacency)
  if (length(adjacency@x) > 0) {
    adjacency@x[] <- 1
  }
  if (nrow(adjacency) > 0) {
    adjacency[cbind(seq_len(nrow(adjacency)), seq_len(nrow(adjacency)))] <- 0
  }
  Matrix::drop0(adjacency)
}

.sn_calculate_connectivity_table <- function(adjacency, groups, group_name) {
  split_idx <- split(seq_along(groups), as.character(groups))
  connectivity_tbl <- lapply(names(split_idx), function(current_group) {
    idx <- split_idx[[current_group]]
    current_adjacency <- adjacency[idx, idx, drop = FALSE]
    largest_component <- .sn_largest_component_size(current_adjacency)
    data.frame(
      group = current_group,
      n_cells = length(idx),
      largest_component = largest_component,
      connectivity_score = largest_component / length(idx),
      stringsAsFactors = FALSE
    )
  })

  connectivity_tbl <- .sn_bind_rows(connectivity_tbl)
  names(connectivity_tbl)[names(connectivity_tbl) == "group"] <- group_name
  connectivity_tbl
}

.sn_largest_component_size <- function(adjacency) {
  adjacency <- methods::as(adjacency, "dgCMatrix")
  n_cells <- nrow(adjacency)
  if (n_cells == 0) {
    return(0L)
  }
  if (n_cells == 1) {
    return(1L)
  }

  col_ptr <- adjacency@p
  row_idx <- adjacency@i
  visited <- rep(FALSE, n_cells)
  max_size <- 0L

  for (start in seq_len(n_cells)) {
    if (visited[[start]]) {
      next
    }

    queue <- integer(n_cells)
    head <- 1L
    tail <- 1L
    queue[[tail]] <- start
    visited[[start]] <- TRUE
    component_size <- 0L

    while (head <= tail) {
      node <- queue[[head]]
      head <- head + 1L
      component_size <- component_size + 1L

      if (col_ptr[[node]] < col_ptr[[node + 1L]]) {
        neighbors <- row_idx[(col_ptr[[node]] + 1L):col_ptr[[node + 1L]]] + 1L
        neighbors <- neighbors[!visited[neighbors]]
        if (length(neighbors) > 0) {
          visited[neighbors] <- TRUE
          queue[(tail + 1L):(tail + length(neighbors))] <- neighbors
          tail <- tail + length(neighbors)
        }
      }
    }

    max_size <- max(max_size, component_size)
  }

  max_size
}

.sn_calculate_neighbor_purity <- function(adjacency, labels, label_name, cell_ids) {
  adjacency <- methods::as(adjacency, "dgCMatrix")
  edge_tbl <- Matrix::summary(adjacency)
  degree <- as.numeric(Matrix::rowSums(adjacency))
  same_group <- labels[edge_tbl$i] == labels[edge_tbl$j]
  same_count <- tabulate(edge_tbl$i[same_group], nbins = nrow(adjacency))
  purity <- ifelse(degree > 0, same_count / degree, NA_real_)

  data.frame(
    cell_id = cell_ids,
    label = as.character(labels),
    neighbor_purity = purity,
    stringsAsFactors = FALSE
  ) |>
    stats::setNames(c("cell_id", label_name, "neighbor_purity"))
}

.sn_scale_lisi <- function(score, n_levels, inverse = FALSE) {
  if (n_levels <= 1) {
    return(1)
  }

  scaled <- (score - 1) / (n_levels - 1)
  scaled <- min(max(scaled, 0), 1)
  if (inverse) {
    scaled <- 1 - scaled
  }
  scaled
}

.sn_scale_silhouette <- function(score) {
  if (is.na(score)) {
    return(NA_real_)
  }
  min(max((score + 1) / 2, 0), 1)
}

.sn_metric_row <- function(metric, category, score, scaled_score, n_cells, source, note = "") {
  data.frame(
    metric = metric,
    category = category,
    score = score,
    scaled_score = scaled_score,
    n_cells = n_cells,
    source = source,
    note = note,
    stringsAsFactors = FALSE
  )
}

.sn_calculate_assessment_aggregates <- function(summary_tbl) {
  if (nrow(summary_tbl) == 0) {
    return(summary_tbl)
  }

  rows <- list()
  batch_rows <- summary_tbl$category == "batch_removal" & !is.na(summary_tbl$scaled_score)
  biology_rows <- summary_tbl$category == "biology_conservation" & !is.na(summary_tbl$scaled_score)
  structure_rows <- summary_tbl$category == "structure" & !is.na(summary_tbl$scaled_score)

  batch_score <- if (any(batch_rows)) {
    mean(summary_tbl$scaled_score[batch_rows])
  } else {
    NA_real_
  }
  biology_score <- if (any(biology_rows)) {
    mean(summary_tbl$scaled_score[biology_rows])
  } else {
    NA_real_
  }
  structure_score <- if (any(structure_rows)) {
    mean(summary_tbl$scaled_score[structure_rows])
  } else {
    NA_real_
  }

  if (!is.na(batch_score)) {
    rows[[length(rows) + 1]] <- .sn_metric_row(
      metric = "batch_mixing_score",
      category = "aggregate",
      score = batch_score,
      scaled_score = batch_score,
      n_cells = max(summary_tbl$n_cells, na.rm = TRUE),
      source = "aggregate"
    )
  }
  if (!is.na(biology_score)) {
    rows[[length(rows) + 1]] <- .sn_metric_row(
      metric = "biology_conservation_score",
      category = "aggregate",
      score = biology_score,
      scaled_score = biology_score,
      n_cells = max(summary_tbl$n_cells, na.rm = TRUE),
      source = "aggregate"
    )
  }
  if (!is.na(structure_score)) {
    rows[[length(rows) + 1]] <- .sn_metric_row(
      metric = "structure_score",
      category = "aggregate",
      score = structure_score,
      scaled_score = structure_score,
      n_cells = max(summary_tbl$n_cells, na.rm = TRUE),
      source = "aggregate"
    )
  }
  if (!is.na(batch_score) && !is.na(biology_score)) {
    rows[[length(rows) + 1]] <- .sn_metric_row(
      metric = "overall_integration_score",
      category = "aggregate",
      score = 0.4 * batch_score + 0.6 * biology_score,
      scaled_score = 0.4 * batch_score + 0.6 * biology_score,
      n_cells = max(summary_tbl$n_cells, na.rm = TRUE),
      source = "0.4 batch_by + 0.6 biology"
    )
  }

  .sn_bind_rows(rows)
}

.sn_bind_rows <- function(rows) {
  rows <- rows[lengths(rows) > 0]
  if (length(rows) == 0) {
    return(data.frame())
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
