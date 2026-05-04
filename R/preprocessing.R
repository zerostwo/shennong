.sn_is_seurat_object <- function(object) {
  methods::is(object, "Seurat")
}

.sn_extract_feature_names <- function(object) {
  if (.sn_is_seurat_object(object)) {
    assay <- SeuratObject::DefaultAssay(object)
    return(rownames(object[[assay]]))
  }

  if (inherits(object, c("matrix", "Matrix", "data.frame"))) {
    return(rownames(object))
  }

  if (is.character(object)) {
    return(object)
  }

  NULL
}

.sn_infer_species_from_features <- function(features) {
  features <- unique(stats::na.omit(as.character(features)))
  features <- features[nzchar(features)]
  if (length(features) == 0) {
    return(NULL)
  }

  homology_env <- new.env(parent = emptyenv())
  utils::data("hom_genes", package = "Shennong", envir = homology_env)
  homology_table <- get("hom_genes", envir = homology_env, inherits = FALSE)

  human_matches <- sum(features %in% homology_table$human)
  mouse_matches <- sum(features %in% homology_table$mouse)

  human_mt <- sum(grepl("^MT-", features))
  mouse_mt <- sum(grepl("^mt-", features))
  human_like <- human_matches + human_mt
  mouse_like <- mouse_matches + mouse_mt

  min_confident_matches <- max(3L, ceiling(length(features) * 0.01))
  if (max(human_like, mouse_like) < min_confident_matches) {
    return(NULL)
  }

  if (human_like == mouse_like) {
    return(NULL)
  }

  if (human_like > mouse_like) "human" else "mouse"
}

#' Retrieve or infer species information
#'
#' This helper returns the explicit \code{species} argument when supplied,
#' otherwise it tries the species stored in a Seurat object and then falls back
#' to feature-name based inference using the packaged \code{hom_genes} mapping
#' plus common mitochondrial naming patterns.
#'
#' @param object A Seurat object, matrix-like object, or character vector of
#'   gene symbols.
#' @param species Optional explicit species label. If provided, it is returned
#'   directly.
#'
#' @return A character string indicating the inferred or explicit species.
#'
#' @examples
#' sn_get_species(c("CD3D", "LTB", "MS4A1"))
#'
#' m <- matrix(0, nrow = 3, ncol = 2)
#' rownames(m) <- c("Cd3d", "Ltb", "Ms4a1")
#' sn_get_species(m)
#' @export
sn_get_species <- function(object, species = NULL) {
  if (!is_null(species)) {
    return(arg_match(species, c("human", "mouse")))
  }

  if (.sn_is_seurat_object(object)) {
    species <- Seurat::Misc(object, slot = "species")
    if (!is_null(species)) {
      return(species)
    }
  }

  inferred_species <- .sn_infer_species_from_features(.sn_extract_feature_names(object))
  if (is_null(inferred_species)) {
    stop(
      paste(
        "Species information is required but could not be inferred from feature names.",
        "Please provide `species = \"human\"` or `species = \"mouse\"`."
      ),
      call. = FALSE
    )
  }

  if (.sn_is_seurat_object(object)) {
    Seurat::Misc(object = object, slot = "species") <- inferred_species
  }

  inferred_species
}

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
#' @examples
#' \dontrun{
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(400 * 20, lambda = 3), nrow = 400, ncol = 20)
#'   rownames(counts) <- c(
#'     "MCM5", "PCNA", "TYMS", "FEN1", "MCM2",
#'     "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
#'     "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS",
#'     "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN",
#'     "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3",
#'     "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45",
#'     "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
#'     "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
#'     "BRIP1", "E2F8",
#'     "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5",
#'     "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
#'     "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3",
#'     "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
#'     "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
#'     "KIF20B", "HJURP", "CDCA3", "CDC20", "TTK",
#'     "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5",
#'     "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
#'     "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
#'     "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3",
#'     "CBX5", "CENPA",
#'     paste0("GENE", 1:316)
#'   )
#'   colnames(counts) <- paste0("cell", 1:20)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_score_cell_cycle(obj, species = "human")
#'   head(obj[[]][, c("S.Score", "G2M.Score", "Phase")])
#' }
#' }
#'
#' @export
sn_score_cell_cycle <- function(object, species = NULL) {
  # Attempt to retrieve species if not provided
  species <- sn_get_species(object = object, species = species)

  # Retrieve S-phase and G2M-phase markers
  s_features <- sn_get_signatures(species = species, category = "Programs/cellCycle.G1S")
  g2m_features <- sn_get_signatures(species = species, category = "Programs/cellCycle.G2M")
  feature_names <- rownames(object)
  s_features <- intersect(s_features, feature_names)
  g2m_features <- intersect(g2m_features, feature_names)

  if (length(s_features) == 0 || length(g2m_features) == 0) {
    .sn_log_warn(
      "Skipping cell cycle scoring because the selected assay has insufficient overlap ",
      "with {species} cell-cycle markers (S: {length(s_features)}, G2M: {length(g2m_features)})."
    )
    return(object)
  }

  # Perform cell cycle scoring
  object <- Seurat::CellCycleScoring(object, s.features = s_features, g2m.features = g2m_features)

  # Compute cell cycle difference
  object$CC.Difference <- object$S.Score - object$G2M.Score
  .sn_log_seurat_command(object = object, name = "sn_score_cell_cycle")
}

#' Initialize a Seurat object with optional QC metrics
#'
#' This function creates a Seurat object from counts (and optional metadata),
#' then calculates common QC metrics such as mitochondrial and ribosomal gene
#' percentages when `species` is supplied. Currently supports human and mouse
#' patterns for these gene sets. When `x` points to a 10x Genomics `outs`
#' directory, the function automatically reads the filtered matrix and stores
#' discovered source metadata such as the raw matrix path and
#' `metrics_summary.csv` contents in `Seurat::Misc(object, "input_source")`.
#' When `x` is a character vector of multiple detected 10x paths, such as the
#' output of [sn_list_10x_paths()], the function returns a named list of Seurat
#' objects and imports each path in one call.
#'
#' @param x A matrix, data.frame, sparse matrix, path to counts data, or a
#'   character vector of multiple 10x paths.
#' @param metadata Optional metadata (data.frame or similar) to add to the Seurat object.
#' @param names_field Passed to \code{SeuratObject::CreateSeuratObject}, indicating how to parse cell names.
#' @param names_delim Passed to \code{SeuratObject::CreateSeuratObject}, indicating the delimiter for cell names.
#' @param min_cells Filter out genes expressed in fewer than \code{min_cells} cells.
#' @param min_features Filter out cells with fewer than \code{min_features} genes.
#' @param project A project name for the Seurat object.
#' @param sample_name Optional sample name to store in \code{meta.data$sample}.
#' @param study Optional study name to store in \code{meta.data$study}.
#' @param species Either "human" or "mouse" (case-sensitive). Affects QC metric patterns.
#' @param standardize_gene_symbols Logical; standardize gene symbols after object creation.
#' @param is_gene_id Logical; treat row names as gene IDs and convert them to symbols when standardizing.
#' @param ... Additional arguments passed to \code{SeuratObject::CreateSeuratObject()}.
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
  sample_name = NULL,
  study = NULL,
  species = NULL,
  standardize_gene_symbols = FALSE,
  is_gene_id = FALSE, ...
) {
  if (is.character(x) && length(x) > 1L) {
    if (!is_null(metadata) && (!is.list(metadata) || is.data.frame(metadata))) {
      stop(
        "When `x` contains multiple paths, `metadata` must be NULL or a list parallel to `x`.",
        call. = FALSE
      )
    }

    .resolve_parallel_arg <- function(value, i, default = NULL) {
      if (is.null(value)) {
        return(default)
      }
      if (length(value) == 1L) {
        return(value[[1]])
      }
      value[[i]]
    }

    object_list <- lapply(seq_along(x), function(i) {
      sample_i <- .resolve_parallel_arg(sample_name, i, default = names(x)[[i]] %||% NULL)
      metadata_i <- if (is.list(metadata)) metadata[[i]] else NULL
      sn_initialize_seurat_object(
        x = x[[i]],
        metadata = metadata_i,
        names_field = names_field,
        names_delim = names_delim,
        project = .resolve_parallel_arg(project, i, default = project),
        min_cells = .resolve_parallel_arg(min_cells, i, default = min_cells),
        min_features = .resolve_parallel_arg(min_features, i, default = min_features),
        sample_name = sample_i,
        study = .resolve_parallel_arg(study, i, default = study),
        species = .resolve_parallel_arg(species, i, default = species),
        standardize_gene_symbols = .resolve_parallel_arg(
          standardize_gene_symbols,
          i,
          default = standardize_gene_symbols
        ),
        is_gene_id = .resolve_parallel_arg(is_gene_id, i, default = is_gene_id),
        ...
      )
    })

    list_names <- names(x)
    if (is.null(list_names) || any(!nzchar(list_names))) {
      fallback_names <- vapply(
        object_list,
        FUN = function(object) unique(as.character(object$sample))[[1]] %||% "",
        FUN.VALUE = character(1)
      )
      fallback_names[!nzchar(fallback_names)] <- paste0("sample_", seq_along(object_list))[!nzchar(fallback_names)]
      list_names <- fallback_names
    }
    names(object_list) <- list_names
    return(object_list)
  }

  # -- Logging
  .sn_log_info("Initializing Seurat object for project: {project}.")

  source_info <- NULL
  inherited_sample_name <- NULL
  if (inherits(x, what = c("Matrix", "matrix", "data.frame", "MatrixDir"))) {
    counts <- x
  } else {
    if (is.character(x) && length(x) == 1L) {
      inherited_sample_name <- names(x)[[1]] %||% NULL
      if (!is.null(inherited_sample_name) && !nzchar(inherited_sample_name)) {
        inherited_sample_name <- NULL
      }
    }
    source_info <- .sn_detect_10x_outs_source(x)
    read_path <- source_info$filtered_path %||% x
    counts <- sn_read(path = read_path)
  }

  if (is_null(sample_name)) {
    sample_name <- inherited_sample_name %||% source_info$sample_name %||% NULL
  }

  counts <- .sn_as_sparse_matrix(counts)

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
  if (!is_null(sample_name)) {
    seurat_obj$sample <- sample_name
  }
  if (!is_null(study)) {
    seurat_obj$study <- study
  }

  if (is_null(species)) {
    species <- tryCatch(
      sn_get_species(object = seurat_obj, species = NULL),
      error = function(e) NULL
    )
    if (is_null(species)) {
      .sn_log_warn("Could not infer species during initialization; skipping species-specific QC metrics.")
    }
  }

  # -- Calculate QC metrics if required
  if (!is_null(species)) {
    Seurat::Misc(object = seurat_obj, slot = "species") <- species

    if (standardize_gene_symbols) {
      seurat_obj <- seurat_obj |>
        sn_standardize_gene_symbols(species = species, is_gene_id = is_gene_id)
    }

    .sn_log_info("Running QC metrics for {species}.")
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
      .sn_log_warn("Unsupported species for QC metrics; skipping QC calculation.")
    }
  }

  if (!is_null(source_info)) {
    Seurat::Misc(object = seurat_obj, slot = "input_source") <- source_info
  }

  .sn_log_info("Seurat object initialization complete.")
  return(.sn_log_seurat_command(object = seurat_obj, name = "sn_initialize_seurat_object"))
}

#' Normalize data in a Seurat object
#'
#' This function provides a unified normalization entry point for Seurat-style
#' log-normalization, scran normalization, and SCTransform.
#'
#' @param object A \code{Seurat} object.
#' @param method One of \code{"seurat"}, \code{"scran"}, or
#'   \code{"sctransform"} (alias \code{"sct"}).
#' @param clusters Optional cluster assignments for \code{scran::quickCluster}.
#' @param assay Assay used for normalization. Defaults to \code{"RNA"}.
#' @param layer Layer used as the input count matrix. Defaults to \code{"counts"}.
#' @param ... Additional method-specific arguments passed to
#'   \code{Seurat::NormalizeData()}, \code{scran::computeSumFactors()}, or
#'   \code{Seurat::SCTransform()}.
#'
#' @return A \code{Seurat} object with normalized data stored according to the
#'   chosen method.
#' @examples
#' \dontrun{
#' seurat_obj <- sn_normalize_data(seurat_obj, method = "scran")
#' }
#' @export
sn_normalize_data <- function(
  object,
  method = c("seurat", "scran", "sctransform", "sct"),
  clusters = NULL,
  assay = "RNA",
  layer = "counts",
  ...
) {
  method <- match.arg(method)
  if (method == "sct") {
    method <- "sctransform"
  }

  if (method == "seurat") {
    prepared <- .sn_prepare_seurat_analysis_input(
      object = object,
      assay = assay,
      layer = layer
    )
    object <- Seurat::NormalizeData(object = prepared$object, ...)
    object <- .sn_restore_seurat_analysis_input(object = object, context = prepared$context)
    return(.sn_log_seurat_command(object = object, assay = assay, name = "sn_normalize_data"))
  }

  if (method == "scran") {
    check_installed("scran", reason = "to perform scran normalization.")
    check_installed("SingleCellExperiment")
    counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)

    .sn_log_info("Converting the Seurat object to SingleCellExperiment for scran normalization.")
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = counts),
      colData = object[[]]
    )

    if (is_null(clusters)) {
      .sn_log_info("Running `scran::quickCluster()` to assign clusters.")
      clusters <- scran::quickCluster(
        x         = sce,
        use.ranks = FALSE,
        min.size  = 50
      )
    }

    .sn_log_info("Computing size factors with `scran::computeSumFactors()`.")
    sce <- scran::computeSumFactors(sce, clusters = clusters, min.mean = 0.1, ...)
    size_factors <- SingleCellExperiment::sizeFactors(object = sce)

    .sn_log_info("Size factor summary: {paste(round(summary(size_factors), 3), collapse = ', ')}.")
    object$size.factor <- size_factors

    # -- Apply normalization
    .sn_log_info("Applying log-normalization.")
    normalized_counts <- sweep(counts, 2, size_factors, "/")
    object$total_counts_normalized <- Matrix::colSums(x = normalized_counts)

    # -- Store normalized data
    SeuratObject::LayerData(
      object = object,
      assay  = assay,
      layer  = "data"
    ) <- methods::as(log(normalized_counts + 1), "CsparseMatrix")

    .sn_log_info("scran normalization complete.")
    return(.sn_log_seurat_command(object = object, assay = assay, name = "sn_normalize_data"))
  }

  check_installed("glmGamPoi", reason = "for the SCTransform workflow.")
  prepared <- .sn_prepare_seurat_analysis_input(
    object = object,
    assay = assay,
    layer = layer
  )
  object <- Seurat::SCTransform(
    object = prepared$object,
    verbose = TRUE,
    seed.use = 717,
    ...
  )

  if (isTRUE(prepared$context$needs_temp_counts)) {
    if (isTRUE(prepared$context$had_exact_counts)) {
      SeuratObject::LayerData(
        object = object,
        assay = prepared$context$analysis_assay,
        layer = "counts"
      ) <- prepared$context$original_counts
    } else {
      SeuratObject::LayerData(
        object = object,
        assay = prepared$context$analysis_assay,
        layer = "counts"
      ) <- NULL
    }
  }

  .sn_log_seurat_command(object = object, assay = assay, name = "sn_normalize_data")
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
  is_gene_id = FALSE
) {
  # -- Check required packages
  check_installed("HGNChelper", reason = "to check or correct gene symbols.")
  check_installed("dplyr")

  .sn_log_info("Starting gene-symbol standardization.")
  # -- Extract counts
  if (inherits(x, "Seurat")) {
    .sn_log_info("Detected a Seurat object; extracting counts from the RNA assay.")
    species <- sn_get_species(object = x, species = species)
    counts <- SeuratObject::LayerData(object = x, layer = "counts")
  } else {
    counts <- x
    species <- sn_get_species(object = counts, species = species)
  }

  # -- Convert to matrix if needed
  if (inherits(counts, "data.frame") || inherits(counts, "DelayedMatrix")) {
    counts <- as.matrix(counts)
  }

  # -- Convert gene IDs to symbols if requested
  if (is_gene_id) {
    gene_ids <- rownames(counts)
    .sn_log_info("Converting gene IDs to gene symbols.")

    genes <- .sn_get_gene_annotation_table(species = species)
    df <- dplyr::tibble(gene_id = gene_ids) |>
      dplyr::left_join(
        dplyr::select(genes, .data$gene_id, .data$gene_id_base, .data$gene_name),
        by = "gene_id"
      ) |>
      dplyr::mutate(
        gene_id_base = sub("\\..*$", "", .data$gene_id),
        gene_name = dplyr::coalesce(
          .data$gene_name,
          genes$gene_name[match(.data$gene_id_base, genes$gene_id_base)]
        )
      )

    keep_genes <- !is_na(df$gene_name)
    if (sum(keep_genes) < length(keep_genes)) {
      .sn_log_warn("Some gene IDs did not match and will be removed.")
    }

    counts <- counts[keep_genes, , drop = FALSE]
    df <- df[keep_genes, ]
    rownames(counts) <- df$gene_name
  }

  # -- Check and correct gene symbols
  .sn_log_info("Checking gene symbols for species = {species} with `HGNChelper`.")
  check_gene_symbols <- HGNChelper::checkGeneSymbols(
    x       = rownames(counts),
    species = species
  )
  stopifnot(identical(check_gene_symbols$x, rownames(counts)))

  check_gene_symbols <- check_gene_symbols |>
    dplyr::mutate(Suggested.Symbol = dplyr::if_else(
      stringr::str_detect(.data$Suggested.Symbol, "///"),
      NA_character_,
      .data$Suggested.Symbol
    ))

  suggested_symbols <- as.character(check_gene_symbols$Suggested.Symbol)
  original_symbols <- as.character(check_gene_symbols$x)
  invalid_suggested <- is.na(suggested_symbols) | !nzchar(suggested_symbols)
  valid_original <- !is.na(original_symbols) & nzchar(original_symbols)
  suggested_symbols[invalid_suggested & valid_original] <- original_symbols[invalid_suggested & valid_original]

  # -- Filter out truly missing symbols, but preserve valid original names when
  # HGNChelper has no unambiguous suggestion.
  valid_idx <- !is.na(suggested_symbols) & nzchar(suggested_symbols)
  counts <- counts[valid_idx, , drop = FALSE]
  genes <- suggested_symbols[valid_idx]
  rownames(counts) <- genes

  # -- Aggregate duplicates
  duplicated_genes <- rownames(counts)[duplicated(rownames(counts)) | duplicated(rownames(counts), fromLast = TRUE)]
  if (length(duplicated_genes) > 0) {
    .sn_log_info("Found {length(unique(duplicated_genes))} duplicated gene symbol(s); aggregating counts.")
    counts <- .sn_aggregate_rows_by_group(
      x = counts,
      groups = rownames(counts)
    )
  }

  .sn_log_info("Gene-symbol standardization complete.")

  # -- Return updated Seurat or matrix
  if (inherits(x, "Seurat")) {
    x[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = counts)
    return(.sn_log_seurat_command(object = x, name = "sn_standardize_gene_symbols"))
  } else {
    return(counts)
  }
}

.sn_get_gene_annotation_table <- function(species = NULL) {
  annotation_env <- new.env(parent = emptyenv())
  utils::data("shennong_gene_annotations", package = "Shennong", envir = annotation_env)
  annotations <- get("shennong_gene_annotations", envir = annotation_env, inherits = FALSE)

  if (!is_null(species)) {
    species <- arg_match(species, c("human", "mouse"))
    annotations <- annotations[annotations$species == species, , drop = FALSE]
  }

  annotations
}

.sn_match_gene_annotations <- function(features, species) {
  annotations <- .sn_get_gene_annotation_table(species = species)
  features <- as.character(features)

  by_name <- annotations[!is.na(annotations$gene_name) & nzchar(annotations$gene_name), , drop = FALSE]
  by_name <- by_name[!duplicated(by_name$gene_name), , drop = FALSE]
  rownames(by_name) <- by_name$gene_name

  by_id <- annotations[!is.na(annotations$gene_id) & nzchar(annotations$gene_id), , drop = FALSE]
  by_id <- by_id[!duplicated(by_id$gene_id), , drop = FALSE]
  rownames(by_id) <- by_id$gene_id

  by_id_base <- annotations[!is.na(annotations$gene_id_base) & nzchar(annotations$gene_id_base), , drop = FALSE]
  by_id_base <- by_id_base[!duplicated(by_id_base$gene_id_base), , drop = FALSE]
  rownames(by_id_base) <- by_id_base$gene_id_base

  feature_id_base <- sub("\\..*$", "", features)
  matched_index <- rep(NA_character_, length(features))
  matched_by <- rep(NA_character_, length(features))

  name_hits <- features %in% rownames(by_name)
  matched_index[name_hits] <- features[name_hits]
  matched_by[name_hits] <- "gene_name"

  id_hits <- is.na(matched_index) & features %in% rownames(by_id)
  matched_index[id_hits] <- features[id_hits]
  matched_by[id_hits] <- "gene_id"

  id_base_hits <- is.na(matched_index) & feature_id_base %in% rownames(by_id_base)
  matched_index[id_base_hits] <- feature_id_base[id_base_hits]
  matched_by[id_base_hits] <- "gene_id_base"

  out <- data.frame(
    feature = features,
    matched = !is.na(matched_index),
    matched_by = matched_by,
    gene_id = NA_character_,
    gene_id_base = NA_character_,
    gene_name = NA_character_,
    gene_type = NA_character_,
    gene_class = NA_character_,
    stringsAsFactors = FALSE
  )

  if (any(out$matched)) {
    out$gene_id[name_hits] <- by_name[matched_index[name_hits], "gene_id", drop = TRUE]
    out$gene_id_base[name_hits] <- by_name[matched_index[name_hits], "gene_id_base", drop = TRUE]
    out$gene_name[name_hits] <- by_name[matched_index[name_hits], "gene_name", drop = TRUE]
    out$gene_type[name_hits] <- by_name[matched_index[name_hits], "gene_type", drop = TRUE]
    out$gene_class[name_hits] <- by_name[matched_index[name_hits], "gene_class", drop = TRUE]

    out$gene_id[id_hits] <- by_id[matched_index[id_hits], "gene_id", drop = TRUE]
    out$gene_id_base[id_hits] <- by_id[matched_index[id_hits], "gene_id_base", drop = TRUE]
    out$gene_name[id_hits] <- by_id[matched_index[id_hits], "gene_name", drop = TRUE]
    out$gene_type[id_hits] <- by_id[matched_index[id_hits], "gene_type", drop = TRUE]
    out$gene_class[id_hits] <- by_id[matched_index[id_hits], "gene_class", drop = TRUE]

    out$gene_id[id_base_hits] <- by_id_base[matched_index[id_base_hits], "gene_id", drop = TRUE]
    out$gene_id_base[id_base_hits] <- by_id_base[matched_index[id_base_hits], "gene_id_base", drop = TRUE]
    out$gene_name[id_base_hits] <- by_id_base[matched_index[id_base_hits], "gene_name", drop = TRUE]
    out$gene_type[id_base_hits] <- by_id_base[matched_index[id_base_hits], "gene_type", drop = TRUE]
    out$gene_class[id_base_hits] <- by_id_base[matched_index[id_base_hits], "gene_class", drop = TRUE]
  }

  out
}

#' Filter genes based on the number of expressing cells
#'
#' This function filters genes in a Seurat object based on the number of cells
#' in which they are expressed. Optionally, it can visualize the effect of
#' different filtering thresholds and retain only genes matching bundled
#' GENCODE-based annotation classes.
#'
#' @param x A Seurat object.
#' @param min_cells An integer specifying the minimum number of cells in which a
#'   gene must be expressed to be retained. Default is 3.
#' @param plot Logical; if TRUE, a bar plot is generated showing the number of
#'   remaining genes at different filtering thresholds. Default is TRUE.
#' @param filter Logical; if TRUE, returns the filtered Seurat object. If
#'   FALSE, returns the original object. Default is TRUE.
#' @param assay Assay used when extracting expression values. Defaults to
#'   \code{"RNA"}.
#' @param layer Layer used for gene filtering. Defaults to \code{"counts"}.
#' @param species Optional species label used for annotation-aware filtering.
#'   Required only when it cannot be inferred from the object.
#' @param gene_class Optional coarse annotation class to retain. One of
#'   \code{"coding"} or \code{"noncoding"}.
#' @param gene_type Optional character vector of exact GENCODE \code{gene_type}
#'   values to retain.
#'
#' @return A filtered Seurat object if \code{filter = TRUE}, otherwise the
#'   original object.
#'
#' @details
#' The function computes the number of cells expressing each gene in the Seurat
#' object and filters out genes expressed in fewer than \code{min_cells} cells.
#' If \code{plot = TRUE}, it visualizes the effect of filtering using
#' \code{ggplot2}. When \code{gene_class} or \code{gene_type} is supplied, the
#' function also matches the feature set against the bundled
#' \code{shennong_gene_annotations} data and keeps only the requested
#' annotation subset.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   data("pbmc_small", package = "Shennong")
#'   pbmc_filtered <- sn_filter_genes(
#'     pbmc_small,
#'     min_cells = 5,
#'     plot = FALSE,
#'     filter = TRUE
#'   )
#'   pbmc_coding <- sn_filter_genes(
#'     pbmc_small,
#'     min_cells = 1,
#'     plot = FALSE,
#'     filter = TRUE,
#'     species = "human",
#'     gene_class = "coding"
#'   )
#' }
#'
#' @export
sn_filter_genes <- function(x,
                            min_cells = 3,
                            plot = TRUE,
                            filter = TRUE,
                            assay = "RNA",
                            layer = "counts",
                            species = NULL,
                            gene_class = NULL,
                            gene_type = NULL) {
  if (!inherits(x, "Seurat")) {
    stop("The input object x is not a Seurat object.")
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1 || is.na(min_cells) ||
    min_cells < 0 || min_cells != as.integer(min_cells)) {
    stop("`min_cells` must be a single non-negative integer.")
  }
  min_cells <- as.integer(min_cells)

  counts <- .sn_get_seurat_layer_data(object = x, assay = assay, layer = layer)
  if (!is_null(gene_class)) {
    gene_class <- arg_match(gene_class, c("coding", "noncoding"))
  }
  if (!is_null(gene_type)) {
    gene_type <- unique(as.character(gene_type))
    gene_type <- gene_type[!is.na(gene_type) & nzchar(gene_type)]
    if (length(gene_type) == 0) {
      gene_type <- NULL
    }
  }

  gene_counts <- Matrix::rowSums(x = counts > 0)
  gene_distribution <- table(factor(x = gene_counts, levels = 0:max(gene_counts)))
  cumulative_genes <- rev(x = cumsum(x = rev(x = gene_distribution)))

  min_cells_seq <- 0:min_cells
  remaining_genes <- numeric(length(min_cells_seq))
  valid_thresholds <- min_cells_seq <= max(gene_counts)
  remaining_genes[valid_thresholds] <- cumulative_genes[min_cells_seq[valid_thresholds] + 1]

  if (plot) {
    plot_data <- data.frame(
      min_cell = factor(x = min_cells_seq, levels = min_cells_seq),
      remaining_gene = remaining_genes
    )

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = .data$min_cell, x = .data$remaining_gene)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::geom_text(ggplot2::aes(label = .data$remaining_gene), hjust = -0.2) +
      ggplot2::labs(
        y = "Filtering threshold (min_cells)",
        x = "Remaining genes",
        title = "Remaining genes at different filtering thresholds"
      ) +
      ggplot2::theme_minimal()
    print(x = p)
  }

  if (filter) {
    keep_genes <- gene_counts >= min_cells

    if (!is_null(gene_class) || !is_null(gene_type)) {
      species <- sn_get_species(x, species = species)
      annotations <- .sn_match_gene_annotations(rownames(x), species = species)

      if (any(!annotations$matched)) {
        unmatched_features <- utils::head(rownames(x)[!annotations$matched], 10)
        unmatched_msg <- if (length(unmatched_features) > 0) {
          paste0(" Examples: ", paste(unmatched_features, collapse = ", "), ".")
        } else {
          ""
        }
        .sn_log_warn(
          "Annotation-based gene filtering could not match {sum(!annotations$matched)} features for species '{species}'. ",
          "Those unmatched features will be dropped.{unmatched_msg}"
        )
      }

      annotation_keep <- annotations$matched
      if (!is_null(gene_class)) {
        annotation_keep <- annotation_keep & annotations$gene_class %in% gene_class
      }
      if (!is_null(gene_type)) {
        annotation_keep <- annotation_keep & annotations$gene_type %in% gene_type
      }

      keep_genes <- keep_genes & annotation_keep
    }

    x <- x[keep_genes, ]
  }

  .sn_log_seurat_command(object = x, assay = assay, name = "sn_filter_genes")
}

#' @title Filter cells in a Seurat object based on QC metrics
#'
#' @description This function filters cells in a Seurat object using Median Absolute Deviation (MAD)
#' to identify outliers across specified metadata features. Supports grouped analysis and provides
#' visual diagnostics.
#'
#' @param x A Seurat object
#' @param features Character vector of metadata column names to use for filtering
#' @param group_by (Optional) Metadata column name to group by for group-wise calculations
#' @param method Outlier detection method (currently only "mad" supported)
#' @param n Numeric threshold(s) for MAD multiplier. Can be single value or vector matching features.
#' @param plot Logical indicating whether to generate QC diagnostic plots
#' @param filter Logical indicating whether to filter out flagged cells
#'
#' @return Seurat object with QC flags in metadata. If filter=TRUE, returns subsetted object.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' pbmc <- sn_filter_cells(pbmc, features = c("nCount_RNA", "nFeature_RNA"))
#'
#' # Grouped analysis
#' pbmc <- sn_filter_cells(pbmc, features = "percent.mt", group_by = "sample")
#'
#' # Custom threshold with plotting
#' pbmc <- sn_filter_cells(pbmc, features = "nCount_RNA", n = 3, plot = TRUE)
#' }
#'
#' @export
sn_filter_cells <- function(
  x,
  features,
  group_by = NULL,
  method = "mad",
  n = 5,
  plot = TRUE,
  filter = TRUE
) {
  if (!inherits(x, "Seurat")) stop("Input must be a Seurat object")
  features <- unique(as.character(features))
  if (length(features) == 0) {
    stop("`features` must contain at least one metadata column.")
  }
  method <- rlang::arg_match(method, values = "mad")
  if (!is_null(group_by) && (length(group_by) != 1 || !is.character(group_by))) {
    stop("`group_by` must be NULL or a single metadata column name.")
  }
  if (!all(features %in% colnames(x[[]]))) {
    stop("Missing features in metadata: ",
      paste(setdiff(features, colnames(x[[]]))),
      collapse = ", "
    )
  }
  if (!is.numeric(n)) stop("n must be numeric")
  if (length(n) != 1 && length(n) != length(features)) {
    stop("n must be length 1 or match length of features")
  }

  if (length(n) == 1) n <- rep(n, length(features))

  for (i in seq_along(features)) {
    x <- .sn_filter_cells_one(
      x = x,
      feature = features[i],
      group_by = group_by,
      method = method,
      n = n[i]
    )
  }

  if (plot) {
    check_installed(c("aplot", "ggrastr"), reason = "to plot QC filtering diagnostics.")
    plot_data <- x[[]]
    qc_data <- Seurat::Misc(x, "qc")

    plots <- lapply(features, function(f) {
      .create_qc_plot(
        metadata = plot_data,
        qc_info = qc_data[[f]],
        feature = f,
        group_by = group_by
      )
    })

    print(aplot::plot_list(gglist = plots))
  }

  if (filter) {
    keep <- Reduce(`&`, lapply(features, function(f) {
      x[[paste0(f, "_qc")]] == "Passed"
    }))
    x <- x[, keep]
  }

  .sn_log_seurat_command(object = x, name = "sn_filter_cells")
}

.sn_filter_cells_one <- function(x, feature, group_by = NULL, method = "mad", n = 5) {
  if (!feature %in% colnames(x[[]])) {
    stop("Feature ", feature, " not found in metadata")
  }

  feature_data <- x[[feature]][, 1]
  if (!is.numeric(feature_data)) {
    stop("Feature ", feature, " must be numeric")
  }

  if (is_null(x = group_by)) {
    qc_df <- .sn_qc_bounds(feature_data, n = n)
    meta <- x[[]]
    x[[paste0(feature, "_qc")]] <- ifelse(
      .sn_is_within_qc_bounds(meta[[feature]], qc_df$l, qc_df$u),
      "Passed", "Failed"
    )
  } else {
    if (!group_by %in% colnames(x[[]])) {
      stop("Grouping variable ", group_by, " not found in metadata")
    }

    qc_df <- x[[]] |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) |>
      dplyr::summarise(
        median = stats::median(.data[[feature]], na.rm = TRUE),
        mad = stats::mad(.data[[feature]], na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        l = pmax(.data$median - n * .data$mad, 0),
        u = .data$median + n * .data$mad
      )

    meta <- x[[]] |>
      tibble::rownames_to_column("barcode") |>
      dplyr::left_join(qc_df, by = group_by) |>
      tibble::column_to_rownames("barcode")

    x[[paste0(feature, "_qc")]] <- ifelse(
      .sn_is_within_qc_bounds(meta[[feature]], meta$l, meta$u),
      "Passed", "Failed"
    )
  }

  qc_list <- Seurat::Misc(x, "qc") %||% list()
  qc_list[[feature]] <- qc_df
  misc_data <- methods::slot(x, "misc")
  misc_data[["qc"]] <- qc_list
  methods::slot(x, "misc") <- misc_data

  x
}

.sn_qc_bounds <- function(x, n = 5) {
  median_value <- stats::median(x, na.rm = TRUE)
  mad_value <- stats::mad(x, na.rm = TRUE)

  data.frame(
    median = median_value,
    mad = mad_value,
    l = max(median_value - n * mad_value, 0),
    u = median_value + n * mad_value
  )
}

.sn_is_within_qc_bounds <- function(x, lower, upper) {
  !is.na(x) & !is.na(lower) & !is.na(upper) & x >= lower & x <= upper
}

.create_qc_plot <- function(metadata, qc_info, feature, group_by) {
  failed_points <- metadata[metadata[[paste0(feature, "_qc")]] == "Failed", , drop = FALSE]
  base_plot <- ggplot2::ggplot(metadata, ggplot2::aes(
    x = if (!is_null(group_by)) .data[[group_by]] else "All",
    y = .data[[feature]]
  )) +
    ggplot2::geom_violin(scale = "width", trim = TRUE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, title = feature)

  if (!is_null(group_by)) {
    base_plot <- base_plot +
      ggplot2::geom_errorbar(
        data = qc_info,
        ggplot2::aes(y = .data$median, ymin = .data$l, ymax = .data$u),
        color = "red", width = 0.2
      )
  } else {
    base_plot <- base_plot +
      ggplot2::geom_errorbar(
        data = qc_info,
        ggplot2::aes(
          x = "All",
          y = .data$median, ymin = .data$l, ymax = .data$u
        ),
        color = "red", width = 0.1, inherit.aes = FALSE
      )
  }

  base_plot +
    ggrastr::geom_jitter_rast(
      data = failed_points,
      size = 0.3, alpha = 0.5, color = "red", width = 0.2
    )
}

.sn_run_scDblFinder <- function(sce, ...) {
  scDblFinder::scDblFinder(sce = sce, ...)
}

#' Find doublets using scDblFinder
#'
#' This function identifies potential doublets in a Seurat object by converting
#' it to a SingleCellExperiment and using the \code{scDblFinder} package.
#'
#' @param object A \code{Seurat} object.
#' @param clusters Optional cluster assignments. If not provided, scDblFinder will attempt automatic clustering.
#' @param group_by An optional metadata column used as the donor or sample grouping.
#' @param dbr_sd A numeric value for adjusting the doublet rate; see \code{scDblFinder} documentation.
#' @param ncores Number of cores to use (for parallel processing).
#' @param assay Assay used for doublet detection. Defaults to \code{"RNA"}.
#' @param layer Layer used as the input count matrix. Defaults to \code{"counts"}.
#' @param min_features Minimum number of detected features required for a cell
#'   to be passed to \code{scDblFinder()}. Defaults to \code{200}. Cells below
#'   this threshold are skipped and retain \code{NA} in the output columns.
#'
#' @return The input Seurat object with two new columns in \code{meta.data}:
#'   \code{scDblFinder.class} and \code{scDblFinder.score} when
#'   \code{layer = "counts"}, or \code{scDblFinder.class_corrected} and
#'   \code{scDblFinder.score_corrected} for non-default corrected layers. Cells
#'   whose selected layer sums to zero or whose detected-feature count is below
#'   \code{min_features} are skipped and retain \code{NA} in the corresponding
#'   output columns.
#' @examples
#' \dontrun{
#' seurat_obj <- sn_find_doublets(
#'   seurat_obj,
#'   clusters = NULL,
#'   group_by = NULL,
#'   dbr_sd = NULL,
#'   ncores = 4
#' )
#' }
#' @export
sn_find_doublets <- function(
  object,
  clusters = NULL,
  group_by = NULL,
  dbr_sd = NULL,
  ncores = 1,
  assay = "RNA",
  layer = "counts",
  min_features = 200
) {
  check_installed("scDblFinder", reason = "to run doublet detection.")
  check_installed("SingleCellExperiment")
  stopifnot(is.numeric(min_features), length(min_features) == 1, min_features >= 0)

  counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  metadata <- object[[]]
  layer_sums <- Matrix::colSums(counts)
  feature_counts <- Matrix::colSums(counts > 0)
  positive_cells <- colnames(counts)[layer_sums > 0]
  keep_cells <- colnames(counts)[layer_sums > 0 & feature_counts >= min_features]
  skipped_zero_cells <- setdiff(colnames(counts), positive_cells)
  skipped_low_feature_cells <- setdiff(positive_cells, keep_cells)

  if (length(keep_cells) == 0) {
    stop(
      glue(
        "No cells passed doublet-detection filtering in assay '{assay}' layer '{layer}'. ",
        "Check zero-count cells or lower `min_features` (current: {min_features})."
      ),
      call. = FALSE
    )
  }

  if (length(skipped_zero_cells) > 0) {
    .sn_log_warn(
      "Skipping {length(skipped_zero_cells)} zero-count cell(s) in assay '{assay}' layer '{layer}' ",
      "before running `scDblFinder()`."
    )
  }
  if (length(skipped_low_feature_cells) > 0) {
    .sn_log_warn(
      "Skipping {length(skipped_low_feature_cells)} cell(s) in assay '{assay}' layer '{layer}' ",
      "with fewer than {min_features} detected feature(s) before running `scDblFinder()`."
    )
  }

  counts <- counts[, keep_cells, drop = FALSE]
  metadata <- metadata[keep_cells, , drop = FALSE]

  .sn_log_info("Converting the Seurat object to SingleCellExperiment for doublet detection.")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = metadata
  )

  if (is.character(clusters) && length(clusters) == 1 && clusters %in% colnames(object[[]])) {
    clusters <- object[[clusters]][keep_cells, 1]
  } else if (!is_null(clusters) && length(clusters) == ncol(object)) {
    clusters <- clusters[match(keep_cells, colnames(object))]
  }

  if (!is_null(group_by) && !group_by %in% colnames(object[[]])) {
    stop(glue("Grouping column '{group_by}' was not found in object metadata."))
  }

  if (is_null(group_by)) {
    .sn_log_info("Running `scDblFinder()` without donor grouping.")
    sce <- .sn_run_scDblFinder(
      sce = sce,
      clusters  = clusters,
      dbr.sd    = dbr_sd
    )
  } else {
    .sn_log_info("Running `scDblFinder()` with donor grouping (parallel with {ncores} cores).")
    sce <- .sn_run_scDblFinder(
      sce = sce,
      samples  = group_by,
      dbr.sd   = dbr_sd,
      BPPARAM  = BiocParallel::MulticoreParam(ncores)
    )
  }

  class_col <- "scDblFinder.class"
  score_col <- "scDblFinder.score"
  if (!identical(layer, "counts")) {
    class_col <- paste0(class_col, "_corrected")
    score_col <- paste0(score_col, "_corrected")
  }

  .sn_format_doublet_classes <- function(values, unresolved = "unresolved") {
    values <- as.character(values)
    values[is.na(values) | !nzchar(values)] <- unresolved
    preferred_levels <- c("singlet", "doublet")
    other_levels <- setdiff(unique(values), preferred_levels)
    factor(values, levels = c(preferred_levels, sort(other_levels)))
  }

  class_values <- rep("unresolved", ncol(object))
  names(class_values) <- colnames(object)
  class_values[keep_cells] <- as.character(sce$scDblFinder.class)

  score_values <- rep(NA_real_, ncol(object))
  names(score_values) <- colnames(object)
  score_values[keep_cells] <- as.numeric(sce$scDblFinder.score)

  object[[class_col]] <- .sn_format_doublet_classes(class_values[colnames(object)])
  object[[score_col]] <- score_values[colnames(object)]

  .sn_log_info("Doublet detection complete.")
  .sn_log_seurat_command(object = object, assay = assay, name = "sn_find_doublets")
}

.sn_resolve_counts_input <- function(x, arg = "x") {
  if (inherits(x, "Seurat")) {
    return(list(
      object = x,
      counts = .sn_as_sparse_matrix(SeuratObject::LayerData(object = x, layer = "counts"))
    ))
  }

  if (is.character(x)) {
    return(list(
      object = NULL,
      counts = .sn_as_sparse_matrix(sn_read(path = x))
    ))
  }

  if (inherits(x, c("Matrix", "matrix", "data.frame", "MatrixDir", "IterableMatrix", "RenameDims"))) {
    return(list(object = NULL, counts = .sn_as_sparse_matrix(x)))
  }

  stop(glue("`{arg}` must be a Seurat object, matrix-like object, or path."))
}

.sn_parse_10x_metrics_summary <- function(path) {
  if (is_null(path) || !file.exists(path)) {
    return(NULL)
  }

  metrics <- tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      if (grepl("no lines available in input", conditionMessage(e), fixed = TRUE)) {
        return(NULL)
      }
      stop(e)
    }
  )
  if (is_null(metrics)) {
    return(NULL)
  }
  if (nrow(metrics) == 0) {
    return(NULL)
  }
  metrics
}

.sn_locate_10x_outs_paths <- function(path) {
  if (!is.character(path) || length(path) != 1 || is.na(path) || !dir.exists(path)) {
    return(NULL)
  }

  normalized_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  candidates <- c(
    normalized_path,
    file.path(normalized_path, "outs")
  )
  candidates <- unique(candidates[dir.exists(candidates)])

  for (candidate in candidates) {
    filtered_dir <- file.path(candidate, "filtered_feature_bc_matrix")
    filtered_h5 <- file.path(candidate, "filtered_feature_bc_matrix.h5")
    raw_dir <- file.path(candidate, "raw_feature_bc_matrix")
    raw_h5 <- file.path(candidate, "raw_feature_bc_matrix.h5")
    metrics_csv <- file.path(candidate, "metrics_summary.csv")
    web_summary <- file.path(candidate, "web_summary.html")

    has_filtered <- dir.exists(filtered_dir) || file.exists(filtered_h5)
    if (!has_filtered) {
      next
    }

    sample_name <- if (identical(basename(candidate), "outs")) {
      basename(dirname(candidate))
    } else {
      basename(candidate)
    }

    return(list(
      sample_name = sample_name,
      outs_path = candidate,
      filtered_path = if (dir.exists(filtered_dir)) filtered_dir else if (file.exists(filtered_h5)) filtered_h5 else NULL,
      raw_path = if (dir.exists(raw_dir)) raw_dir else if (file.exists(raw_h5)) raw_h5 else NULL,
      metrics_path = if (file.exists(metrics_csv)) metrics_csv else NULL,
      web_summary_path = if (file.exists(web_summary)) web_summary else NULL
    ))
  }

  NULL
}

.sn_detect_10x_outs_source <- function(path) {
  outs_info <- .sn_locate_10x_outs_paths(path)
  if (is_null(outs_info)) {
    return(NULL)
  }

  metrics <- .sn_parse_10x_metrics_summary(outs_info$metrics_path)

  list(
    type = "10x_outs",
    sample_name = outs_info$sample_name %||% NULL,
    outs_path = outs_info$outs_path,
    filtered_path = outs_info$filtered_path,
    raw_path = outs_info$raw_path,
    metrics_path = outs_info$metrics_path,
    web_summary_path = outs_info$web_summary_path,
    metrics = metrics
  )
}

.sn_resolve_stored_raw_input <- function(object, method, verbose = FALSE) {
  if (is_null(object) || !inherits(object, "Seurat")) {
    return(NULL)
  }

  source_info <- Seurat::Misc(object, slot = "input_source")
  raw_path <- source_info$raw_path %||% NULL
  if (is_null(raw_path) || !file.exists(raw_path)) {
    return(NULL)
  }

  .sn_log_info(
    "No `raw` argument supplied; using the stored raw matrix at '{raw_path}' for `{method}`."
  )
  .sn_resolve_counts_input(raw_path, arg = "raw")
}

.sn_restore_count_shape <- function(original_counts, corrected_counts) {
  if (identical(rownames(original_counts), rownames(corrected_counts)) &&
      identical(colnames(original_counts), colnames(corrected_counts))) {
    return(corrected_counts)
  }

  restored_counts <- original_counts
  restored_counts[rownames(corrected_counts), colnames(corrected_counts)] <- corrected_counts
  restored_counts
}

.sn_handle_zero_count_cells <- function(original_counts,
                                        corrected_counts,
                                        remove_zero_count_cells = FALSE) {
  original_sums <- Matrix::colSums(original_counts)
  corrected_sums <- Matrix::colSums(corrected_counts)
  zero_cells <- colnames(corrected_counts)[original_sums > 0 & corrected_sums == 0]

  if (length(zero_cells) == 0) {
    return(list(
      counts = corrected_counts,
      zero_cells = character(0),
      removed_cells = character(0)
    ))
  }

  if (remove_zero_count_cells) {
    keep_cells <- setdiff(colnames(corrected_counts), zero_cells)
    .sn_log_warn(
      "decontX produced {length(zero_cells)} zero-count cell(s); removing them because ",
      "`remove_zero_count_cells = TRUE`."
    )
    return(list(
      counts = corrected_counts[, keep_cells, drop = FALSE],
      zero_cells = zero_cells,
      removed_cells = zero_cells
    ))
  }

  corrected_counts[, zero_cells] <- original_counts[, zero_cells, drop = FALSE]
  .sn_log_warn(
    "decontX produced {length(zero_cells)} zero-count cell(s); restored the original counts ",
    "for those cells. Set `remove_zero_count_cells = TRUE` to drop them instead."
  )

  list(
    counts = corrected_counts,
    zero_cells = zero_cells,
    removed_cells = character(0)
  )
}

.sn_resolve_ambient_clusters <- function(x_info, cluster, verbose = FALSE) {
  counts <- x_info$counts

  if (!is_null(cluster)) {
    if (is.character(cluster) && length(cluster) == 1 && !is_null(x_info$object)) {
      if (!cluster %in% colnames(x_info$object[[]])) {
        stop(glue("Cluster column '{cluster}' was not found in object metadata."))
      }
      cluster <- x_info$object[[cluster]][, 1]
    }

    if (length(cluster) != ncol(counts)) {
      stop("`cluster` must have one value per cell.")
    }

    cluster <- as.character(cluster)
    names(cluster) <- colnames(counts)
    return(cluster)
  }

  cluster_object <- x_info$object %||% sn_initialize_seurat_object(x = counts)
  cluster <- sn_run_cluster(
    object = cluster_object,
    return_cluster = TRUE,
    resolution = 2,
    block_genes = NULL,
    verbose = verbose
  )
  cluster <- as.character(cluster)
  names(cluster) <- colnames(counts)
  cluster
}

.sn_remove_ambient_soupx <- function(x_info,
                                     raw_info,
                                     cluster = NULL,
                                     force_accept = FALSE,
                                     contamination_range = c(0.01, 0.8),
                                     verbose = FALSE) {
  if (is_null(raw_info)) {
    stop("`raw` is required when `method = \"soupx\"`.")
  }

  check_installed(pkg = "SoupX", reason = "to remove ambient RNA contamination with SoupX.")

  tod <- raw_info$counts
  toc <- x_info$counts

  common_genes <- intersect(rownames(tod), rownames(toc))
  tod_common <- tod[common_genes, , drop = FALSE]
  toc_common <- toc[common_genes, , drop = FALSE]

  if (verbose) {
    .sn_log_info(
      "SoupX gene statistics: raw = {nrow(tod)} genes; filtered = {nrow(toc)} genes; ",
      "common = {length(common_genes)} genes."
    )
  }

  sc <- SoupX::SoupChannel(
    tod = tod_common,
    toc = toc_common,
    calcSoupProfile = FALSE
  )

  soup_profile <- data.frame(
    row.names = rownames(toc_common),
    est = Matrix::rowSums(toc_common) / sum(toc_common),
    counts = Matrix::rowSums(toc_common)
  )
  sc <- SoupX::setSoupProfile(sc, soup_profile)

  soupx_groups <- .sn_resolve_ambient_clusters(
    x_info = list(object = x_info$object, counts = toc_common),
    cluster = cluster,
    verbose = verbose
  )
  sc <- SoupX::setClusters(sc, soupx_groups)

  sc <- SoupX::autoEstCont(
    sc,
    doPlot = FALSE,
    forceAccept = force_accept,
    contaminationRange = contamination_range
  )
  out <- SoupX::adjustCounts(sc, roundToInt = TRUE)

  tod_sum <- sum(Matrix::rowSums(tod_common))
  toc_sum <- sum(Matrix::rowSums(toc_common))
  out_sum <- sum(Matrix::rowSums(out))
  if (verbose) {
    .sn_log_info(
      "SoupX count summary: raw = {format(tod_sum, big.mark = ',')}; ",
      "filtered = {format(toc_sum, big.mark = ',')}; output = {format(out_sum, big.mark = ',')}."
    )
  }

  list(
    counts = .sn_restore_count_shape(toc, out),
    metadata = NULL
  )
}

.sn_remove_ambient_decontx <- function(x_info,
                                       raw_info = NULL,
                                       cluster = NULL,
                                       remove_zero_count_cells = FALSE,
                                       verbose = FALSE,
                                       ...) {
  check_installed(pkg = "celda", reason = "to remove ambient RNA contamination with decontX.")
  check_installed(pkg = "SingleCellExperiment", reason = "to run decontX.")

  counts <- x_info$counts
  background_counts <- NULL

  if (!is_null(raw_info)) {
    common_genes <- intersect(rownames(counts), rownames(raw_info$counts))
    counts <- counts[common_genes, , drop = FALSE]
    background_counts <- raw_info$counts[common_genes, , drop = FALSE]
  }

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  background <- NULL
  if (!is_null(background_counts)) {
    background <- SingleCellExperiment::SingleCellExperiment(list(counts = background_counts))
  }

  z <- .sn_resolve_ambient_clusters(
    x_info = list(object = x_info$object, counts = counts),
    cluster = cluster,
    verbose = verbose
  )

  sce <- celda::decontX(
    x = sce,
    z = z,
    background = background,
    ...
  )
  out <- round(celda::decontXcounts(sce))
  out <- .sn_handle_zero_count_cells(
    original_counts = counts,
    corrected_counts = out,
    remove_zero_count_cells = remove_zero_count_cells
  )
  metadata <- as.data.frame(
    SummarizedExperiment::colData(sce)[, c("decontX_contamination", "decontX_clusters"), drop = FALSE]
  )
  metadata <- metadata[colnames(out$counts), , drop = FALSE]
  metadata$nCount_corrected <- Matrix::colSums(out$counts)
  metadata$nFeature_corrected <- Matrix::colSums(out$counts > 0)
  base_counts <- counts
  if (length(out$removed_cells) > 0) {
    base_counts <- base_counts[, colnames(out$counts), drop = FALSE]
  }

  list(
    counts = .sn_restore_count_shape(base_counts, out$counts),
    metadata = metadata,
    zero_cells = out$zero_cells,
    removed_cells = out$removed_cells
  )
}

.sn_apply_ambient_result_to_object <- function(object, out, assay, layer) {
  ncount_col <- paste0("nCount_", assay, "_corrected")
  nfeature_col <- paste0("nFeature_", assay, "_corrected")
  zero_flag_col <- paste0(layer, "_zero_count")

  if (length(out$removed_cells) > 0) {
    object <- object[, colnames(out$counts), drop = FALSE]
  }

  if (!is_null(out$metadata)) {
    if (all(c("nCount_corrected", "nFeature_corrected") %in% colnames(out$metadata))) {
      colnames(out$metadata)[colnames(out$metadata) == "nCount_corrected"] <- ncount_col
      colnames(out$metadata)[colnames(out$metadata) == "nFeature_corrected"] <- nfeature_col
    }
    out$metadata[[zero_flag_col]] <- rownames(out$metadata) %in% (out$zero_cells %||% character(0))
    for (col_name in colnames(out$metadata)) {
      values <- out$metadata[[col_name]][match(colnames(object), rownames(out$metadata))]
      object[[col_name]] <- values
    }
  }

  SeuratObject::LayerData(object = object, layer = layer) <- out$counts
  object
}

#' Remove ambient RNA contamination from counts
#'
#' This function provides a unified interface for ambient RNA correction using
#' either \code{SoupX} or \code{decontX}. The input can be a Seurat object, a
#' count matrix-like object, or a path that \code{sn_read()} can import. When a
#' Seurat object was initialized from a detected 10x `outs` directory, stored
#' raw matrix metadata is reused automatically if \code{raw = NULL}.
#'
#' @param x A Seurat object, count matrix-like object, or path to filtered data.
#' @param raw Optional raw/background counts. Required for \code{method = "soupx"}
#'   unless recoverable from stored initialization metadata. If supplied for
#'   \code{decontx}, it is used as the background matrix.
#' @param method One of \code{"decontx"} or \code{"soupx"}.
#' @param cluster Optional cluster labels. This can be a vector with one value
#'   per cell, or a metadata column name when \code{x} is a Seurat object. If
#'   \code{NULL}, clusters are inferred with \code{sn_run_cluster()}.
#' @param remove_zero_count_cells Logical; if \code{TRUE}, remove cells whose
#'   decontX-corrected counts sum to zero. If \code{FALSE}, keep those cells by
#'   restoring their original counts and emit a warning.
#' @param layer Layer name used when writing corrected counts back to a Seurat
#'   object. Defaults to \code{"decontaminated_counts"}. Use \code{"counts"} to
#'   overwrite the original counts layer explicitly.
#' @param return_object If \code{TRUE} and \code{x} is a Seurat object, return
#'   the updated Seurat object. Otherwise return the corrected counts matrix.
#' @param verbose Logical; whether to print progress from helper clustering.
#' @param ... Additional method-specific arguments passed to
#'   \code{celda::decontX()} when \code{method = "decontx"}, or to
#'   \code{SoupX::autoEstCont()} when \code{method = "soupx"}.
#'
#' @return A corrected counts matrix, or an updated Seurat object when
#'   \code{return_object = TRUE} and \code{x} is a Seurat object. For
#'   decontX-based Seurat returns, the \code{decontX_contamination} and
#'   \code{decontX_clusters} columns are added to \code{meta.data}.
#'
#' @examples
#' \dontrun{
#' data("pbmc_small", package = "Shennong")
#' data("pbmc_small_raw", package = "Shennong")
#' pbmc <- sn_remove_ambient_contamination(pbmc, method = "decontx")
#'
#' corrected <- sn_remove_ambient_contamination(
#'   x = SeuratObject::LayerData(pbmc, layer = "counts"),
#'   raw = pbmc_small_raw,
#'   method = "soupx",
#'   return_object = FALSE
#' )
#' }
#'
#' @export
sn_remove_ambient_contamination <- function(
  x,
  raw = NULL,
  method = c("decontx", "soupx"),
  cluster = NULL,
  remove_zero_count_cells = FALSE,
  layer = "decontaminated_counts",
  return_object = TRUE,
  verbose = FALSE,
  ...
) {
  check_installed(pkg = "SeuratObject", reason = "to handle Seurat objects.")

  method <- match.arg(method)
  x_info <- .sn_resolve_counts_input(x, arg = "x")
  raw_info <- if (is_null(raw)) NULL else .sn_resolve_counts_input(raw, arg = "raw")
  if (is_null(raw_info)) {
    raw_info <- .sn_resolve_stored_raw_input(
      object = x_info$object,
      method = method,
      verbose = verbose
    )
  }

  out <- switch(
    method,
    soupx = .sn_remove_ambient_soupx(
      x_info = x_info,
      raw_info = raw_info,
      cluster = cluster,
      verbose = verbose,
      ...
    ),
    decontx = .sn_remove_ambient_decontx(
      x_info = x_info,
      raw_info = raw_info,
      cluster = cluster,
      remove_zero_count_cells = remove_zero_count_cells,
      verbose = verbose,
      ...
    )
  )

  if (return_object && !is_null(x_info$object)) {
    assay <- SeuratObject::DefaultAssay(x_info$object)
    x_info$object <- .sn_apply_ambient_result_to_object(
      object = x_info$object,
      out = out,
      assay = assay,
      layer = layer
    )
    return(.sn_log_seurat_command(object = x_info$object, name = "sn_remove_ambient_contamination"))
  }

  out$counts
}
