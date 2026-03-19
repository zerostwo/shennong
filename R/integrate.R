# Internal helper to normalize the clustering dimension arguments.
.sn_resolve_cluster_dims <- function(dims = NULL, npcs = 50) {
  dims <- dims %||% seq_len(min(20, npcs))
  dims <- as.integer(dims)
  dims[dims > 0]
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

  split_objects <- Seurat::SplitObject(object, split.by = split_by)
  feature_lists <- lapply(split_objects, function(current_object) {
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
#' @param block_genes Either a character vector of predefined sets (e.g. \code{c("ribo","mito")})
#'   if \code{SignatuR} is installed, or a custom vector of gene symbols to exclude from HVGs.
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
    if (verbose) log_info(glue("[sn_run_cluster] Starting integration for batch='{batch}'..."))
  }

  if (verbose) {
    log_info(glue("[sn_run_cluster] Normalization method = {normalization_method}; batch = {batch %||% 'none'}"))
  }

  predefined_genesets <- c(
    "tcr", "immunoglobulins", "ribo", "mito",
    "heatshock", "noncoding", "pseudogenes", "g1s", "g2m"
  )

  if (normalization_method == "sctransform" && !is_null(block_genes)) {
    log_warn("`block_genes` is only applied in the log-normalization workflows; ignoring it for SCTransform.")
    block_genes <- NULL
  }

  if (!is_null(block_genes)) {
    species <- sn_get_species(object = object, species = species)
    if (verbose) log_info("Processing block genes...")

    if (is_character(block_genes) && all(block_genes %in% predefined_genesets)) {
      block_genes <- sn_get_signatures(species = species, category = block_genes)
      if (verbose) {
        log_info(glue("  Loaded {length(block_genes)} blocked genes from built-in sets."))
      }
    } else {
      checked <- HGNChelper::checkGeneSymbols(block_genes, species = species)
      valid_genes <- checked$Suggested.Symbol[!is_na(checked$Suggested.Symbol)]
      invalid <- setdiff(block_genes, valid_genes)

      if (length(invalid) > 0) {
        log_warn(glue(
          "Removed {length(invalid)} invalid genes from block list (e.g. {paste(utils::head(invalid, 3), collapse=', ')})."
        ))
      }
      block_genes <- unique(valid_genes)
      if (verbose) log_info(glue("  Using {length(block_genes)} custom block genes."))
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

    if (verbose) log_info("[1/4] Running SCTransform...")
    object <- Seurat::SCTransform(
      object = object,
      variable.features.n = nfeatures,
      vars.to.regress = vars_to_regress,
      verbose = verbose,
      seed.use = 717
    )

    if (verbose) log_info("[2/4] Running PCA...")
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
      if (verbose) log_info("[1/6] Cell cycle scoring...")
      object <- sn_score_cell_cycle(object = object, species = species)
    }

    if (verbose) {
      log_info(glue("[2/6] Selecting highly variable features with hvg_group_by = {hvg_group_by %||% 'global'}..."))
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
        log_info(glue(
          "  Removed {n_removed} genes ({round(n_removed/n_before*100,1)}%) via block_genes"
        ))
      }

      if (length(hvg) < nfeatures) {
        log_warn(glue(
          "Only {length(hvg)} HVGs left (< requested {nfeatures}).\n",
          "Consider adjusting 'nfeatures' or 'block_genes'."
        ))
      }
      hvg <- utils::head(hvg, nfeatures)
    }
    Seurat::VariableFeatures(object = object) <- hvg

    #-- Scaling
    if (verbose) log_info("[3/6] Scaling data...")
    object <- Seurat::ScaleData(
      object = object,
      vars.to.regress = vars_to_regress,
      features = hvg,
      verbose = verbose
    )

    #-- PCA
    if (verbose) log_info("[4/6] Running PCA...")
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
      if (verbose) log_info("[5/6] Running Harmony integration...")
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

  if (verbose) log_info("[6/6] Clustering with integrated embeddings...")
  object <- Seurat::FindNeighbors(object, reduction = reduction, dims = dims, verbose = verbose)
  object <- Seurat::FindClusters(object, resolution = resolution, random.seed = 717, verbose = verbose)

  if (return_cluster) {
    object <- .sn_restore_seurat_analysis_input(object = object, context = prepared$context)
    if (verbose) log_info("Integration completed successfully!")
    return(object@meta.data[, "seurat_clusters"])
  } else {
    if (verbose) log_info("[7/7] Running UMAP...")
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
    if (verbose) log_info("Integration completed successfully!")
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

  log_info("Starting CellTypist analysis with model: {model}")
  tictoc::tic("Total CellTypist runtime")

  if (is_null(outdir)) {
    outdir <- tempfile("celltypist_")
    dir.create(outdir, recursive = TRUE)
    log_info("Using temporary output directory: {outdir}")
  } else {
    outdir <- sn_set_path(outdir)
    log_info("Using user-specified output directory: {outdir}")
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
    log_info("Converting Seurat object to CellTypist input format")
    input_data <- file.path(outdir, "counts.csv")

    counts <- tryCatch(
      .sn_get_seurat_layer_data(object = x, assay = assay, layer = layer),
      error = function(e) {
        log_error("Failed to extract counts layer: {e$message}")
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
    log_info("Using precomputed input matrix: {input_data}")
  }

  over_clustering_path <- NULL

  if (over_clustering != "auto" && over_clustering %in% colnames(x@meta.data)) {
    log_info("Using existing clustering column: {over_clustering}")
    over_clustering_path <- file.path(outdir, "over_clustering.txt")
    writeLines(
      as.character(x[[over_clustering, drop = TRUE]]),
      over_clustering_path
    )
  } else if (over_clustering == "auto") {
    log_warn("Automatic over-clustering not yet implemented")
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
  log_info("Executing CellTypist with command:\n{celltypist} {paste(cmd_args, collapse=' ')}")

  exit_code <- system2(
    command = celltypist,
    args = cmd_args,
    stdout = if (quiet) FALSE else "",
    stderr = if (quiet) FALSE else ""
  )

  if (exit_code != 0) {
    log_error("CellTypist failed with exit code {exit_code}")
    stop("CellTypist execution failed. Check logs for details.")
  }
  log_debug("CellTypist output written to {outdir}")

  if (!file.exists(predicted_labels_path)) {
    log_error("Prediction output missing: {predicted_labels_path}")
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

  log_info("Adding {ncol(predicted_labels)} metadata columns to Seurat object")
  x <- SeuratObject::AddMetaData(x, metadata = predicted_labels)

  tictoc::toc()
  log_info("CellTypist analysis completed successfully")

  .sn_log_seurat_command(object = x, assay = assay, name = "sn_run_celltypist")
}
