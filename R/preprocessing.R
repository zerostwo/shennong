#' Score Cell Cycle Phases
#'
#' This function scores the cell cycle phase of single-cell RNA-seq data using S and G2M phase marker genes.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param species (Optional) A character string indicating the species (e.g., "human" or "mouse").
#'                If NULL, the function will attempt to retrieve species information from `Seurat::Misc(object)`.
#'
#' @return A Seurat object with cell cycle scores added, including:
#'         - `S.Score`: Score for the S phase
#'         - `G2M.Score`: Score for the G2M phase
#'         - `Phase`: Assigned cell cycle phase
#'         - `CC.Difference`: Difference between S and G2M scores
#'
#' @export
sn_score_cell_cycle <- function(object, species = NULL) {
  # Attempt to retrieve species if not provided
  species <- sn_get_species(object = object, species = species)

  # Retrieve S-phase and G2M-phase markers
  s_features <- sn_get_signatures(species = species, category = "g1s")
  g2m_features <- sn_get_signatures(species = species, category = "g2m")

  # Perform cell cycle scoring
  object <- Seurat::CellCycleScoring(object, s.features = s_features, g2m.features = g2m_features)

  # Compute cell cycle difference
  object$CC.Difference <- object$S.Score - object$G2M.Score

  return(object)
}

#' Initialize a Seurat object with optional QC metrics
#'
#' This function creates a Seurat object from counts (and optional metadata),
#' then calculates common QC metrics such as mitochondrial and ribosomal gene
#' percentages if `qc = TRUE`. Currently supports human and mouse patterns
#' for these gene sets.
#'
#' @param counts A matrix, data.frame, or sparse matrix of counts.
#' @param metadata Optional metadata (data.frame or similar) to add to the Seurat object.
#' @param names_field Passed to \code{SeuratObject::CreateSeuratObject}, indicating how to parse cell names.
#' @param names_delim Passed to \code{SeuratObject::CreateSeuratObject}, indicating the delimiter for cell names.
#' @param min_cells Filter out genes expressed in fewer than \code{min_cells} cells.
#' @param min_features Filter out cells with fewer than \code{min_features} genes.
#' @param project A project name for the Seurat object.
#' @param species Either "human" or "mouse" (case-sensitive). Affects QC metric patterns.
#' @param qc Whether to calculate mitochondrial, ribosomal, and hemoglobin gene percentages.
#'
#' @return A \code{Seurat} object with optional QC metadata in its meta.data slot.
#' @examples
#' \dontrun{
#' # Minimal example:
#' seurat_obj <- sn_initialize_seurat_object(x = my_counts, project = "ExampleProject")
#' }
#' @export
sn_initialize_seurat_object <- function(
    x,
    metadata = NULL,
    names_field = 1L,
    names_delim = "_",
    project = "Shennong",
    min_cells = 0,
    min_features = 0,
    species = NULL,
    standardize_gene_symbols = FALSE,
    is_gene_id = FALSE, ...) {
  # -- Logging
  log_info(glue("Initializing Seurat object for project: {project}"))

  if (inherits(x, what = c("Matrix", "data.frame", "MatrixDir"))) {
    counts <- x
  } else {
    counts <- sn_read(path = x)
  }

  # -- Create Seurat Object
  seurat_obj <- SeuratObject::CreateSeuratObject(
    counts       = counts,
    meta.data    = metadata,
    project      = project,
    names.field  = names_field,
    names.delim  = names_delim,
    min.cells    = min_cells,
    min.features = min_features,
    ...
  )
  # -- Add metadata fields
  seurat_obj$study <- project

  # -- Calculate QC metrics if required
  if (!is_null(species)) {
    Seurat::Misc(object = seurat_obj, slot = "species") <- species

    if (standardize_gene_symbols) {
      seurat_obj <- seurat_obj |>
        sn_standardize_gene_symbols(species = species, is_gene_id = is_gene_id)
    }

    log_info(glue("Running QC metrics for {species} ..."))
    mt_featuers <- sn_get_signatures(species = species, category = "mito")
    mt_featuers <- mt_featuers[mt_featuers %in% rownames(x = seurat_obj)]
    ribo_features <- sn_get_signatures(species = species, category = "ribo")
    ribo_features <- ribo_features[ribo_features %in% rownames(x = seurat_obj)]
    if (species == "human") {
      seurat_obj <- seurat_obj |>
        Seurat::PercentageFeatureSet(features = mt_featuers, col.name = "percent.mt") |>
        Seurat::PercentageFeatureSet(features = ribo_features, col.name = "percent.ribo") |>
        Seurat::PercentageFeatureSet(pattern = "^HB[^(P)]", col.name = "percent.hb")
    } else if (species == "mouse") {
      seurat_obj <- seurat_obj |>
        Seurat::PercentageFeatureSet(features = mt_featuers, col.name = "percent.mt") |>
        Seurat::PercentageFeatureSet(features = ribo_features, col.name = "percent.ribo") |>
        Seurat::PercentageFeatureSet(pattern = "^Hb[^(p)]", col.name = "percent.hb")
    } else {
      log_warn("Unsupported species for QC metrics. Skipping QC calculation.")
    }
  }

  log_info("Seurat object initialization complete.")
  return(seurat_obj)
}

#' Normalize data in a Seurat object using scran
#'
#' This function converts a Seurat object to a SingleCellExperiment,
#' computes size factors with \code{scran}, and then applies log-normalization.
#'
#' @param object A \code{Seurat} object.
#' @param method Currently only "scran" is supported.
#' @param clusters Optional cluster assignments for \code{scran::quickCluster}.
#'
#' @return A \code{Seurat} object with normalized data in the \code{"data"} layer.
#' @examples
#' \dontrun{
#' seurat_obj <- sn_normalize_data(seurat_obj, method = "scran")
#' }
#' @export
sn_normalize_data <- function(
    object,
    method = "scran",
    clusters = NULL) {
  if (method == "scran") {
    check_installed("scran", reason = "to perform scran normalization.")
    check_installed("SingleCellExperiment")

    log_info("Converting Seurat object to SingleCellExperiment for scran normalization...")
    sce <- Seurat::as.SingleCellExperiment(x = object)

    if (is_null(clusters)) {
      log_info("Using scran::quickCluster to assign clusters.")
      clusters <- scran::quickCluster(
        x         = sce,
        use.ranks = FALSE,
        min.size  = 50
      )
    }

    log_info("Computing size factors via scran::computeSumFactors...")
    sce <- scran::computeSumFactors(sce, clusters = clusters, min.mean = 0.1)
    size_factors <- SingleCellExperiment::sizeFactors(object = sce)

    # log_info(glue("Size factor summary: {summary(size_factors)}"))
    log_info(glue("Size factor summary: {paste(round(summary(size_factors), 3), collapse = ', ')}"))
    object$size.factor <- size_factors

    # -- Apply normalization
    log_info("Applying log-normalization...")
    counts <- SeuratObject::LayerData(object = object, layer = "counts")
    normalized_counts <- sweep(counts, 2, size_factors, "/")
    object$total_counts_normalized <- Matrix::colSums(x = normalized_counts)

    # -- Store normalized data
    SeuratObject::LayerData(
      object = object,
      layer  = "data"
    ) <- as(log(normalized_counts + 1), "CsparseMatrix")

    log_info("scran normalization complete.")
    return(object)
  } else {
    stop("Currently only 'scran' method is supported in sn_normalize_data().")
  }
}

#' Standardize gene symbols in a count matrix or Seurat object
#'
#' This function helps unify gene symbols to a standard format. It can:
#' \enumerate{
#'   \item Convert gene IDs to gene symbols (if \code{is_gene_id = TRUE}).
#'   \item Check and correct gene symbols using \code{HGNChelper::checkGeneSymbols}.
#'   \item Aggregate duplicated gene symbols by summing their counts.
#' }
#'
#' @param x A count matrix or a \code{Seurat} object.
#' @param species The species for gene symbol checking, passed to \code{HGNChelper}.
#' @param is_gene_id If \code{TRUE}, \code{x} is assumed to contain gene IDs (e.g., ENSEMBL IDs) as rownames.
#'
#' @return If \code{x} is a matrix, returns a matrix with standardized gene symbols.
#'         If \code{x} is a Seurat object, returns the modified Seurat object.
#' @examples
#' \dontrun{
#' # For a Seurat object:
#' seurat_obj <- sn_standardize_gene_symbols(seurat_obj, species = "human", is_gene_id = FALSE)
#'
#' # For a matrix:
#' counts_mat <- sn_standardize_gene_symbols(counts_mat, species = "human", is_gene_id = TRUE)
#' }
#' @export
sn_standardize_gene_symbols <- function(
    x,
    species = NULL,
    is_gene_id = FALSE) {
  # -- Check required packages
  check_installed("HGNChelper", reason = "to check or correct gene symbols.")
  check_installed("dplyr")
  check_installed("readr")

  log_info("Starting standardization of gene symbols...")

  # -- Extract counts
  if (inherits(x, "Seurat")) {
    log_info("Detected Seurat object. Extracting counts from RNA assay...")
    species <- sn_get_species(object = x, species = species)
    counts <- SeuratObject::LayerData(object = x, layer = "counts")
  } else {
    counts <- x
  }

  # -- Convert to matrix if needed
  if (inherits(counts, "data.frame") || inherits(counts, "DelayedMatrix")) {
    counts <- as.matrix(counts)
  }

  # -- Convert gene IDs to symbols if requested
  if (is_gene_id) {
    gene_ids <- rownames(counts)
    log_info("Converting gene IDs to symbols...")

    genes <- readr::read_csv("/mnt/reference_genomes/gencode/human/47/genes.csv")
    df <- dplyr::tibble(gene_id = gene_ids) |>
      dplyr::left_join(genes, by = "gene_id")

    keep_genes <- !is_na(df$gene_name)
    if (sum(keep_genes) < length(keep_genes)) {
      log_warn("Some gene IDs did not match. They will be removed.")
    }

    counts <- counts[keep_genes, , drop = FALSE]
    df <- df[keep_genes, ]
    rownames(counts) <- df$gene_name
  }

  # -- Check and correct gene symbols
  log_info(glue("Using HGNChelper to check gene symbols for species: {species}"))
  check_gene_symbols <- HGNChelper::checkGeneSymbols(
    x       = rownames(counts),
    species = species
  )
  stopifnot(identical(check_gene_symbols$x, rownames(counts)))

  check_gene_symbols <- check_gene_symbols |>
    dplyr::mutate(Suggested.Symbol = dplyr::if_else(
      stringr::str_detect(Suggested.Symbol, "///"),
      NA_character_,
      Suggested.Symbol
    ))

  # -- Filter out NA suggested symbols
  valid_idx <- !is_na(check_gene_symbols$Suggested.Symbol)
  counts <- counts[valid_idx, , drop = FALSE]
  genes <- check_gene_symbols$Suggested.Symbol[valid_idx]
  rownames(counts) <- genes

  # -- Aggregate duplicates
  duplicated_genes <- rownames(counts)[duplicated(rownames(counts)) | duplicated(rownames(counts), fromLast = TRUE)]
  if (length(duplicated_genes) > 0) {
    log_info(glue("Found {length(unique(duplicated_genes))} duplicated gene symbol(s). Aggregating..."))

    duplicate_rows <- counts[duplicated_genes, , drop = FALSE]
    counts <- counts[!(rownames(counts) %in% duplicated_genes), , drop = FALSE]

    row_names <- rownames(duplicate_rows)
    aggregated_rows <- apply(
      X = duplicate_rows,
      MARGIN = 2,
      FUN = function(col) tapply(col, row_names, sum, na.rm = TRUE)
    )

    counts <- rbind(counts, aggregated_rows)
  }

  log_info("Gene symbol standardization complete.")

  # -- Return updated Seurat or matrix
  if (inherits(x, "Seurat")) {
    x[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = counts)
    return(x)
  } else {
    return(counts)
  }
}
