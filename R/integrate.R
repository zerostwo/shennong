#' Quickly run a clustering pipeline (standard or SCTransform-based)
#'
#' This function provides a simplified clustering workflow for a Seurat object,
#' including data normalization, variable feature selection, scaling, PCA,
#' neighbor finding, clustering, and optional UMAP embedding.
#'
#' @param object A \code{Seurat} object.
#' @param pipeline One of \code{"standard"} or \code{"sctransform"} (alias:
#'   \code{"sct"}).
#' @param nfeatures Number of variable features to find.
#' @param vars_to_regress Covariates for regression in \code{ScaleData} or
#'   \code{SCTransform}.
#' @param resolution Resolution parameter for \code{FindClusters}.
#' @param dims Dimensions to use for PCA/UMAP/FindNeighbors.
#' @param return_cluster If \code{TRUE}, returns just the cluster vector instead
#'   of the entire Seurat object.
#' @param verbose Whether to print/log progress messages.
#'
#' @return A \code{Seurat} object with clustering and UMAP performed, or a
#'   cluster vector if \code{return_cluster=TRUE}.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- sn_quick_cluster(seurat_obj, pipeline = "standard", resolution = 0.8)
#' }
#' @export
sn_quick_cluster <- function(object,
                             pipeline = "standard",
                             nfeatures = 3000,
                             vars_to_regress = NULL,
                             resolution = 0.8,
                             dims = 1:30,
                             return_cluster = FALSE,
                             verbose = TRUE,
                             species = NULL) {
  pipeline <- match.arg(pipeline, c("standard", "sctransform", "sct"))

  if (verbose) logger::log_info(glue::glue("[sn_quick_cluster] Pipeline = {pipeline}"))

  #-- Standard pipeline
  if (pipeline == "standard") {
    if (verbose) logger::log_info("Running standard normalization + variable feature selection + scaling...")
    object <- object |>
      Seurat::NormalizeData(verbose = verbose) |>
      Seurat::FindVariableFeatures(nfeatures = nfeatures, verbose = verbose) |>
      Seurat::ScaleData(vars.to.regress = vars_to_regress, verbose = verbose)
  } else {
    # sctransform pipeline
    rlang::check_installed("glmGamPoi", reason = "for SCTransform workflow (suggested in Seurat).")

    if (verbose) logger::log_info("Running SCTransform pipeline...")
    object <- Seurat::SCTransform(
      object              = object,
      variable.features.n = nfeatures,
      vars.to.regress     = vars_to_regress,
      verbose             = verbose,
      seed.use            = 717
    )
  }

  #-- PCA / Neighbors / Clusters
  if (verbose) logger::log_info("Running PCA, FindNeighbors, and FindClusters...")
  object <- object |>
    Seurat::RunPCA(seed.use = 717, verbose = verbose) |>
    Seurat::FindNeighbors(dims = dims, verbose = verbose) |>
    Seurat::FindClusters(resolution = resolution, random.seed = 717, verbose = verbose)

  if (return_cluster) {
    return(object@meta.data[, "seurat_clusters"])
  } else {
    if (verbose) logger::log_info("Running UMAP for visualization...")
    object <- object |>
      Seurat::RunUMAP(dims = dims, seed.use = 717, verbose = verbose)

    return(object)
  }
}

#' Integrate multiple batches using Harmony + optional blocklist genes
#'
#' This function demonstrates a workflow to integrate multiple batches
#' using \code{Harmony}, with optional blocklist-based HVG filtering.
#' It also does cell cycle scoring, HVG detection, scaling, PCA, clustering,
#' and UMAP. Compatible with Seurat v5 functions for integration.
#'
#' @param object A \code{Seurat} object.
#' @param batch A column name in \code{object@meta.data} specifying batch info.
#' @param nfeatures Number of variable features to select.
#' @param vars_to_regress Covariates to regress out in \code{ScaleData}.
#' @param resolution Resolution parameter for \code{FindClusters}.
#' @param block_genes Either a character vector of predefined sets (e.g. \code{c("ribo","mito")})
#'   if \code{SignatuR} is installed, or a custom vector of gene symbols to exclude from HVGs.
#' @param theta The \code{theta} parameter for \code{harmony::RunHarmony}, controlling batch
#'   diversity preservation vs. correction.
#' @param npcs Number of PCs to compute in \code{RunPCA}.
#' @param dims_use A numeric vector of PCs (dimensions) to use for \code{Harmony} and subsequent steps.
#' @param verbose Whether to print/log progress messages.
#'
#' @return A \code{Seurat} object with integrated data (Harmony reduction), cell cycle scores,
#'   clusters, and UMAP embeddings.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- sn_run_cluster(
#'   object = seurat_obj,
#'   batch = "sample_id",
#'   nfeatures = 3000,
#'   resolution = 0.5,
#'   block_genes = c("ribo", "mito") # or a custom vector of gene symbols
#' )
#' }
#' @export
sn_run_cluster <- function(object,
                           batch = NULL,
                           nfeatures = 3000,
                           vars_to_regress = NULL,
                           resolution = 0.8,
                           block_genes = c("ribo", "mito", "tcr", "immunoglobulins", "noncoding", "pseudogenes"),
                           theta = 2,
                           npcs = 50,
                           dims_use = 1:20,
                           species = NULL,
                           return_cluster = FALSE,
                           verbose = TRUE) {
  rlang::check_installed("Seurat")
  rlang::check_installed("HGNChelper")
  rlang::check_installed("harmony")

  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (!is.null(x = batch)) {
    if (!(batch %in% colnames(object@meta.data))) {
      stop(glue::glue("Batch variable '{batch}' not found in metadata."))
    }
    if (verbose) logger::log_info(glue::glue("[sn_quick_integration] Starting integration for batch='{batch}'..."))
  }


  predefined_genesets <- c(
    "tcr", "immunoglobulins", "ribo", "mito",
    "heatshock", "noncoding", "pseudogenes", "g1s", "g2m"
  )

  if (!is.null(block_genes)) {
    species <- sn_get_species(object = object, species = species)
    if (verbose) logger::log_info("Processing block genes...")

    if (is.character(block_genes) && all(block_genes %in% predefined_genesets)) {
      block_genes <- sn_get_signatures(species = species, category = block_genes)
      if (verbose) {
        logger::log_info(glue::glue("  Loaded {length(block_genes)} blocked genes from built-in sets."))
      }
    } else {
      checked <- HGNChelper::checkGeneSymbols(block_genes, species = species)
      valid_genes <- checked$Suggested.Symbol[!is.na(checked$Suggested.Symbol)]
      invalid <- setdiff(block_genes, valid_genes)

      if (length(invalid) > 0) {
        logger::log_warn(glue::glue(
          "Removed {length(invalid)} invalid genes from block list (e.g. {paste(head(invalid, 3), collapse=', ')})."
        ))
      }
      block_genes <- unique(valid_genes)
      if (verbose) logger::log_info(glue::glue("  Using {length(block_genes)} custom block genes."))
    }
  }
  object <- Seurat::NormalizeData(object = object, verbose = verbose)

  if (!is.null(species)) {
    if (verbose) logger::log_info("[1/8] Cell cycle scoring...")
    object <- sn_score_cell_cycle(object = object, species = species)
  }
  #-- （Seurat v5）Split by batch & Join
  if (!is.null(x = batch)) {
    if (verbose) logger::log_info("[2/8] Splitting object by batch and preparing for integration (Seurat v5 style)...")
    object[["RNA"]] <- split(object[["RNA"]], f = object[[batch, drop = TRUE]])
  }
  if (verbose) logger::log_info("[3/8] Selecting highly variable features...")
  object <- Seurat::FindVariableFeatures(object, nfeatures = nfeatures, verbose = verbose)

  if (is.null(x = batch)) {
    hvg <- Seurat::VariableFeatures(object = object)
  } else {
    hvg <- Seurat::SelectIntegrationFeatures5(object, nfeatures = nfeatures + 1000, verbose = verbose)
    object <- SeuratObject::JoinLayers(object)
  }

  if (!is.null(block_genes)) {
    n_before <- length(hvg)
    hvg <- setdiff(hvg, block_genes)
    n_removed <- n_before - length(hvg)

    if (verbose) {
      logger::log_info(glue::glue(
        "  Removed {n_removed} genes ({round(n_removed/n_before*100,1)}%) via block_genes"
      ))
    }

    if (length(hvg) < nfeatures) {
      logger::log_warn(glue::glue(
        "Only {length(hvg)} HVGs left (< requested {nfeatures}).\n",
        "Consider adjusting 'nfeatures' or 'block_genes'."
      ))
    }
    hvg <- head(hvg, nfeatures)
  }
  Seurat::VariableFeatures(object = object) <- hvg

  #-- Scaling
  if (verbose) logger::log_info("[4/8] Scaling data...")
  object <- Seurat::ScaleData(
    object          = object,
    vars.to.regress = vars_to_regress,
    features        = hvg,
    verbose         = verbose
  )

  #-- PCA
  if (verbose) logger::log_info("[5/8] Running PCA...")
  object <- Seurat::RunPCA(
    object,
    npcs     = npcs,
    features = hvg,
    verbose  = verbose
  )

  if (is.null(x = batch)) {
    reduction <- "pca"
  } else {
    #-- Harmony
    if (verbose) logger::log_info("[6/8] Running Harmony integration...")
    object <- harmony::RunHarmony(
      object        = object,
      group.by.vars = batch,
      theta         = theta,
      reduction.use = "pca",
      verbose       = verbose
    )
    reduction <- "harmony"
  }

  if (verbose) logger::log_info("[7/8] Clustering with integrated embeddings...")
  object <- Seurat::FindNeighbors(object, reduction = reduction, dims = dims_use, verbose = verbose)
  object <- Seurat::FindClusters(object, resolution = resolution, verbose = verbose)

  if (return_cluster) {
    if (verbose) logger::log_info("Integration completed successfully!")
    return(object@meta.data[, "seurat_clusters"])
  } else {
    #-- UMAP
    if (verbose) logger::log_info("[8/8] Running UMAP...")
    object <- Seurat::RunUMAP(object, reduction = reduction, dims = dims_use, verbose = verbose)
    if (verbose) logger::log_info("Integration completed successfully!")
    return(object)
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
                              xlsx = FALSE,
                              plot_results = FALSE,
                              quiet = FALSE) {
  rlang::check_installed(c("Seurat", "data.table", "logger", "glue"),
    reason = "to run CellTypist analysis."
  )

  mode <- match.arg(mode)
  celltypist <- celltypist %||% default_celltypist_path
  sn_check_file(celltypist)

  logger::log_info("Starting CellTypist analysis with model: {model}")
  tictoc::tic("Total CellTypist runtime")

  if (is.null(outdir)) {
    outdir <- tempfile("celltypist_")
    dir.create(outdir, recursive = TRUE)
    logger::log_info("Using temporary output directory: {outdir}")
  } else {
    outdir <- sn_set_path(outdir)
    logger::log_info("Using user-specified output directory: {outdir}")
  }

  on.exit(
    {
      if (dir.exists(outdir) && grepl("^celltypist_", basename(outdir))) {
        unlink(outdir, recursive = TRUE)
        logger::log_debug("Cleaned temporary directory: {outdir}")
      }
    },
    add = TRUE
  )

  if (inherits(x, "Seurat")) {
    logger::log_info("Converting Seurat object to CellTypist input format")
    input_data <- file.path(outdir, "counts.csv")

    counts <- tryCatch(
      SeuratObject::LayerData(object = x, layer = "counts"),
      error = function(e) {
        logger::log_error("Failed to extract counts layer: {e$message}")
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
    logger::log_debug("Count matrix written to {input_data} ({file.size(input_data)} bytes)")
  } else {
    input_data <- x
    logger::log_info("Using precomputed input matrix: {input_data}")
  }

  over_clustering_path <- NULL

  if (over_clustering != "auto" && over_clustering %in% colnames(x@meta.data)) {
    logger::log_info("Using existing clustering column: {over_clustering}")
    over_clustering_path <- file.path(outdir, "over_clustering.txt")
    writeLines(
      as.character(x[[over_clustering, drop = TRUE]]),
      over_clustering_path
    )
  } else if (over_clustering == "auto") {
    logger::log_warn("Automatic over-clustering not yet implemented")
  }

  model_name <- tools::file_path_sans_ext(basename(model))
  prefix <- prefix %||% glue::glue("{model_name}.")
  predicted_labels_path <- file.path(outdir, glue::glue("{prefix}predicted_labels.csv"))

  cmd_args <- c(
    "--indata", shQuote(input_data),
    "--model", shQuote(model),
    "--mode", mode,
    "--outdir", shQuote(outdir),
    "--prefix", prefix
  )

  add_arg <- function(args, flag, value, condition = TRUE) {
    if (condition && !is.null(value)) c(args, flag, shQuote(value)) else args
  }

  cmd_args <- add_arg(cmd_args, "--gene-file", gene_file, !is.null(gene_file))
  cmd_args <- add_arg(cmd_args, "--cell-file", cell_file, !is.null(cell_file))
  cmd_args <- add_arg(cmd_args, "--p-thres", p_thres, mode == "prob_match")
  cmd_args <- add_arg(cmd_args, "--over-clustering", over_clustering_path, !is.null(over_clustering_path))
  cmd_args <- add_arg(cmd_args, "--min-prop", min_prop, majority_voting)

  if (transpose_input) cmd_args <- c(cmd_args, "--transpose-input")
  if (xlsx) cmd_args <- c(cmd_args, "--xlsx")
  if (plot_results) cmd_args <- c(cmd_args, "--plot-results")
  if (quiet) cmd_args <- c(cmd_args, "--quiet")
  if (majority_voting) cmd_args <- c(cmd_args, "--majority-voting")
  logger::log_info("Executing CellTypist with command:\n{celltypist} {paste(cmd_args, collapse=' ')}")

  exit_code <- system2(
    command = celltypist,
    args = cmd_args,
    stdout = if (quiet) FALSE else "",
    stderr = if (quiet) FALSE else ""
  )

  if (exit_code != 0) {
    logger::log_error("CellTypist failed with exit code {exit_code}")
    stop("CellTypist execution failed. Check logs for details.")
  }
  logger::log_debug("CellTypist output written to {outdir}")

  if (!file.exists(predicted_labels_path)) {
    logger::log_error("Prediction output missing: {predicted_labels_path}")
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

  logger::log_info("Adding {ncol(predicted_labels)} metadata columns to Seurat object")
  x <- SeuratObject::AddMetaData(x, metadata = predicted_labels)

  tictoc::toc()
  logger::log_info("CellTypist analysis completed successfully")

  return(x)
}
