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
#' seurat_obj <- sn_initialize_seurat_object(counts = my_counts, project = "ExampleProject")
#' }
#' @export
sn_initialize_seurat_object <- function(
    counts,
    metadata = NULL,
    names_field = 1L,
    names_delim = "_",
    min_cells = 0,
    min_features = 0,
    project = "SeuratProject",
    species = "human",
    qc = TRUE) {
  # -- Logging
  logger::log_info(glue::glue("Initializing Seurat object for project: {project}"))

  # -- Create Seurat Object
  seurat_obj <- SeuratObject::CreateSeuratObject(
    counts       = counts,
    meta.data    = metadata,
    project      = project,
    names.field  = names_field,
    names.delim  = names_delim,
    min.cells    = min_cells,
    min.features = min_features
  )

  # -- Calculate QC metrics if required
  if (qc) {
    logger::log_info(glue::glue("Running QC metrics for {species} ..."))

    if (species == "human") {
      seurat_obj <- seurat_obj |>
        Seurat::PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") |>
        Seurat::PercentageFeatureSet(pattern = "^RPS|^RPL", col.name = "percent.ribo") |>
        Seurat::PercentageFeatureSet(pattern = "^HB[^(P)]", col.name = "percent.hb")
    } else if (species == "mouse") {
      seurat_obj <- seurat_obj |>
        Seurat::PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt") |>
        Seurat::PercentageFeatureSet(pattern = "^Rps|^Rpl", col.name = "percent.ribo") |>
        Seurat::PercentageFeatureSet(pattern = "^Hb[^(p)]", col.name = "percent.hb")
    } else {
      logger::log_warn("Unsupported species for QC metrics. Skipping QC calculation.")
    }
  }

  # -- Add metadata fields
  seurat_obj$study <- project
  seurat_obj@misc$species <- species

  logger::log_info("Seurat object initialization complete.")
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
    rlang::check_installed("scran", reason = "to perform scran normalization.")
    rlang::check_installed("SingleCellExperiment")

    logger::log_info("Converting Seurat object to SingleCellExperiment for scran normalization...")
    sce <- Seurat::as.SingleCellExperiment(x = object)

    if (rlang::is_null(clusters)) {
      logger::log_info("Using scran::quickCluster to assign clusters.")
      clusters <- scran::quickCluster(
        x         = sce,
        use.ranks = FALSE,
        min.size  = 50
      )
    }

    logger::log_info("Computing size factors via scran::computeSumFactors...")
    sce <- scran::computeSumFactors(sce, clusters = clusters, min.mean = 0.1)
    size_factors <- SingleCellExperiment::sizeFactors(object = sce)

    # logger::log_info(glue::glue("Size factor summary: {summary(size_factors)}"))
    logger::log_info(glue::glue("Size factor summary: {paste(round(summary(size_factors), 3), collapse = ', ')}"))
    object$size.factor <- size_factors

    # -- Apply normalization
    logger::log_info("Applying log-normalization...")
    counts <- SeuratObject::LayerData(object = object, layer = "counts")
    normalized_counts <- sweep(counts, 2, size_factors, "/")
    object$total_counts_normalized <- Matrix::colSums(x = normalized_counts)

    # -- Store normalized data
    SeuratObject::LayerData(
      object = object,
      layer  = "data"
    ) <- as(log(normalized_counts + 1), "CsparseMatrix")

    logger::log_info("scran normalization complete.")
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
    species = "human",
    is_gene_id = FALSE) {
  # -- Check required packages
  rlang::check_installed("HGNChelper", reason = "to check or correct gene symbols.")
  rlang::check_installed("dplyr")
  rlang::check_installed("readr")

  logger::log_info("Starting standardization of gene symbols...")

  # -- Extract counts
  if (inherits(x, "Seurat")) {
    logger::log_info("Detected Seurat object. Extracting counts from RNA assay...")
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
    logger::log_info("Converting gene IDs to symbols...")

    genes <- readr::read_csv("/mnt/reference_genomes/gencode/human/47/genes.csv")
    df <- dplyr::tibble(gene_id = gene_ids) |>
      dplyr::left_join(genes, by = "gene_id")

    keep_genes <- !is.na(df$gene_name)
    if (sum(keep_genes) < length(keep_genes)) {
      logger::log_warn("Some gene IDs did not match. They will be removed.")
    }

    counts <- counts[keep_genes, , drop = FALSE]
    df <- df[keep_genes, ]
    rownames(counts) <- df$gene_name
  }

  # -- Check and correct gene symbols
  logger::log_info(glue::glue("Using HGNChelper to check gene symbols for species: {species}"))
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
  valid_idx <- !is.na(check_gene_symbols$Suggested.Symbol)
  counts <- counts[valid_idx, , drop = FALSE]
  genes <- check_gene_symbols$Suggested.Symbol[valid_idx]
  rownames(counts) <- genes

  # -- Aggregate duplicates
  duplicated_genes <- rownames(counts)[duplicated(rownames(counts)) | duplicated(rownames(counts), fromLast = TRUE)]
  if (length(duplicated_genes) > 0) {
    logger::log_info(glue::glue("Found {length(unique(duplicated_genes))} duplicated gene symbol(s). Aggregating..."))

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

  logger::log_info("Gene symbol standardization complete.")

  # -- Return updated Seurat or matrix
  if (inherits(x, "Seurat")) {
    x[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = counts)
    return(x)
  } else {
    return(counts)
  }
}
