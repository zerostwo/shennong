.sn_cnv_reference_cells <- function(object,
                                    reference_cells = NULL,
                                    reference_by = NULL,
                                    reference_cat = NULL) {
  cells <- colnames(object)
  if (!is_null(reference_cells)) {
    if (is.logical(reference_cells)) {
      if (length(reference_cells) != length(cells)) {
        stop("Logical `reference_cells` must have one value per cell.", call. = FALSE)
      }
      reference_cells <- cells[reference_cells]
    }
    reference_cells <- unique(as.character(reference_cells))
    missing <- setdiff(reference_cells, cells)
    if (length(missing) > 0L) {
      stop("`reference_cells` contains cells absent from `object`: ", paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE)
    }
  } else if (!is_null(reference_by)) {
    if (!reference_by %in% colnames(object[[]])) {
      stop("`reference_by` column '", reference_by, "' was not found.", call. = FALSE)
    }
    observed <- as.character(object[[reference_by, drop = TRUE]])
    if (is_null(reference_cat) || length(reference_cat) == 0L) {
      stop("`reference_cat` is required when `reference_by` is supplied.", call. = FALSE)
    }
    missing <- setdiff(as.character(reference_cat), unique(observed))
    if (length(missing) > 0L) {
      stop("Reference categories were not observed: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    reference_cells <- cells[observed %in% reference_cat]
  }
  if (length(reference_cells) == 0L) {
    stop("Supply at least one normal cell through `reference_cells` or `reference_by`/`reference_cat`.", call. = FALSE)
  }
  reference_cells
}

.sn_cnv_long_chromosomes <- function(table) {
  if (is_null(table) || nrow(table) == 0L) return(tibble::tibble())
  table <- as.data.frame(table, check.names = FALSE)
  if ("cell" %in% names(table) && all(c("chromosome", "cnv") %in% names(table))) {
    return(tibble::as_tibble(table))
  }
  if (!"cell" %in% names(table)) {
    table$cell <- rownames(table)
  }
  chromosomes <- setdiff(names(table), "cell")
  dplyr::bind_rows(lapply(chromosomes, function(chromosome_name) {
    values <- suppressWarnings(as.numeric(table[[chromosome_name]]))
    tibble::tibble(
      cell = as.character(table$cell),
      chromosome = chromosome_name,
      cnv = values
    )
  }))
}

.sn_cnv_scores_from_chromosomes <- function(chromosomes) {
  if (nrow(chromosomes) == 0L) return(numeric())
  scores <- vapply(split(chromosomes, chromosomes$cell), function(table) {
    values <- table$cnv[is.finite(table$cnv)]
    if (length(values) == 0L) NA_real_ else mean(abs(values))
  }, numeric(1))
  scores
}

.sn_cnv_malignancy <- function(scores, reference_cells, threshold = 2) {
  reference_scores <- scores[intersect(names(scores), reference_cells)]
  reference_scores <- reference_scores[is.finite(reference_scores)]
  if (length(reference_scores) == 0L) {
    stop("No finite CNV scores were available for reference cells.", call. = FALSE)
  }
  center <- stats::median(reference_scores)
  spread <- stats::mad(reference_scores, center = center)
  if (!is.finite(spread) || spread <= .Machine$double.eps) spread <- stats::sd(reference_scores)
  if (!is.finite(spread) || spread <= .Machine$double.eps) spread <- stats::sd(scores, na.rm = TRUE)
  if (!is.finite(spread) || spread <= .Machine$double.eps) spread <- 1
  malignant_score <- (scores - center) / spread
  call <- ifelse(is.finite(malignant_score) & malignant_score >= threshold, "malignant", "non_malignant")
  call[names(scores) %in% reference_cells] <- "reference"
  list(score = malignant_score, call = call, center = center, spread = spread, threshold = threshold)
}

.sn_cnv_subclones <- function(chromosomes,
                              malignant_calls,
                              existing = NULL,
                              k = 2L) {
  cells <- names(malignant_calls)
  subclone <- stats::setNames(rep(NA_character_, length(cells)), cells)
  malignant <- cells[malignant_calls == "malignant"]
  if (!is_null(existing)) {
    existing <- as.character(existing[cells])
    keep <- malignant[!is.na(existing[malignant]) & nzchar(existing[malignant])]
    subclone[keep] <- paste0("subclone_", existing[keep])
    malignant <- setdiff(malignant, keep)
  }
  if (length(malignant) == 1L) subclone[malignant] <- "subclone_1"
  if (length(malignant) >= 2L && nrow(chromosomes) > 0L) {
    wide <- stats::reshape(
      as.data.frame(chromosomes[chromosomes$cell %in% malignant, c("cell", "chromosome", "cnv")]),
      idvar = "cell", timevar = "chromosome", direction = "wide"
    )
    rownames(wide) <- wide$cell
    matrix <- as.matrix(wide[, setdiff(names(wide), "cell"), drop = FALSE])
    matrix[!is.finite(matrix)] <- 0
    groups <- stats::cutree(stats::hclust(stats::dist(matrix)), k = min(as.integer(k), nrow(matrix)))
    subclone[names(groups)] <- paste0("subclone_", groups)
  }
  subclone[malignant_calls == "reference"] <- "reference"
  subclone
}

.sn_cnv_expression_association <- function(object,
                                           scores,
                                           assay = NULL,
                                           layer = "data",
                                           n_features = 50L) {
  expression <- .sn_annotation_expression(object, assay = assay, layer = layer)$matrix
  cells <- intersect(colnames(expression), names(scores)[is.finite(scores)])
  if (length(cells) < 4L) return(tibble::tibble())
  matrix <- expression[, cells, drop = FALSE]
  means <- Matrix::rowMeans(matrix)
  variances <- Matrix::rowMeans(matrix ^ 2) - means ^ 2
  features <- names(utils::head(sort(variances, decreasing = TRUE), as.integer(n_features)))
  values <- as.matrix(matrix[features, , drop = FALSE])
  associations <- lapply(seq_len(nrow(values)), function(index) {
    estimate <- suppressWarnings(stats::cor(values[index, ], scores[cells], method = "spearman", use = "complete.obs"))
    p_value <- if (sum(is.finite(values[index, ]) & is.finite(scores[cells])) >= 4L) {
      suppressWarnings(stats::cor.test(values[index, ], scores[cells], method = "spearman", exact = FALSE)$p.value)
    } else {
      NA_real_
    }
    tibble::tibble(feature = rownames(values)[[index]], correlation = estimate, p_value = p_value)
  })
  table <- dplyr::bind_rows(associations)
  table$adjusted_p_value <- stats::p.adjust(table$p_value, method = "BH")
  table[order(abs(table$correlation), decreasing = TRUE), , drop = FALSE]
}

.sn_run_cnv_infercnvpy <- function(object,
                                   reference_cells,
                                   genome,
                                   store_name,
                                   assay,
                                   layer,
                                   sample_by,
                                   backend_control) {
  if (is.function(backend_control$runner)) {
    return(backend_control$runner(
      object = object, reference_cells = reference_cells, genome = genome,
      store_name = store_name, assay = assay, layer = layer, sample_by = sample_by
    ))
  }
  object$.sn_cnv_reference <- ifelse(colnames(object) %in% reference_cells, "reference", "query")
  prefix <- backend_control$metadata_prefix %||% "infercnvpy_"
  controls <- backend_control
  controls$runner <- NULL
  controls$metadata_prefix <- NULL
  controls <- utils::modifyList(list(
    object = object,
    assay = assay,
    layer = layer,
    species = genome,
    reference_by = ".sn_cnv_reference",
    reference_cat = "reference",
    run_umap = TRUE,
    metadata_prefix = prefix,
    result_name = store_name,
    return_object = TRUE
  ), controls, keep.null = TRUE)
  updated <- do.call(sn_run_infercnvpy, controls)
  updated$.sn_cnv_reference <- NULL
  metadata <- updated[[]]
  score_column <- .sn_communication_column(metadata, c(paste0(prefix, "cnv_score"), paste0(prefix, "score")))
  cluster_column <- .sn_communication_column(metadata, c(paste0(prefix, "cnv_leiden"), paste0(prefix, "leiden")))
  manifest <- updated@misc$infercnvpy[[store_name]]
  chromosome_path <- file.path(manifest$output_dir, "cnv_chromosome.csv")
  chromosomes <- if (file.exists(chromosome_path)) {
    .sn_cnv_long_chromosomes(utils::read.csv(chromosome_path, row.names = 1, check.names = FALSE))
  } else {
    tibble::tibble()
  }
  scores <- if (!is_null(score_column)) {
    stats::setNames(as.numeric(metadata[[score_column]]), rownames(metadata))
  } else {
    .sn_cnv_scores_from_chromosomes(chromosomes)
  }
  existing <- if (is_null(cluster_column)) NULL else stats::setNames(metadata[[cluster_column]], rownames(metadata))
  embedding <- NULL
  umap_name <- paste0(prefix, "cnv_umap")
  if (umap_name %in% names(updated@reductions)) {
    embedding <- SeuratObject::Embeddings(updated[[umap_name]])
  }
  list(
    object = updated, scores = scores, chromosomes = chromosomes,
    subclone = existing, prediction = NULL, embedding = embedding,
    artifacts = manifest
  )
}

.sn_copykat_chromosomes <- function(cna, cells = NULL) {
  cna <- as.data.frame(cna, check.names = FALSE)
  chromosome_column <- .sn_communication_column(cna, c("chrom", "chromosome", "chr"))
  if (is_null(chromosome_column)) return(tibble::tibble())
  coordinate_columns <- intersect(c("chrom", "chromosome", "chr", "start", "end", "gene.name", "gene", "cytoband", "band"), names(cna))
  cell_columns <- if (is_null(cells)) setdiff(names(cna), coordinate_columns) else intersect(cells, names(cna))
  if (length(cell_columns) == 0L) return(tibble::tibble())
  rows <- lapply(unique(as.character(cna[[chromosome_column]])), function(chromosome) {
    keep <- as.character(cna[[chromosome_column]]) == chromosome
    values <- as.matrix(cna[keep, cell_columns, drop = FALSE])
    storage.mode(values) <- "numeric"
    tibble::tibble(cell = cell_columns, chromosome = chromosome, cnv = Matrix::colMeans(values))
  })
  dplyr::bind_rows(rows)
}

.sn_run_cnv_copykat <- function(object,
                                reference_cells,
                                genome,
                                store_name,
                                assay,
                                layer,
                                sample_by,
                                backend_control) {
  if (is.function(backend_control$runner)) {
    return(backend_control$runner(
      object = object, reference_cells = reference_cells, genome = genome,
      store_name = store_name, assay = assay, layer = layer, sample_by = sample_by
    ))
  }
  check_installed("copykat", reason = "to run CopyKAT CNV inference.")
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  counts <- SeuratObject::LayerData(object, assay = assay, layer = layer %||% "counts")
  samples <- if (is_null(sample_by)) stats::setNames(rep("all", ncol(object)), colnames(object)) else {
    if (!sample_by %in% colnames(object[[]])) stop("`sample_by` column was not found.", call. = FALSE)
    stats::setNames(as.character(object[[sample_by, drop = TRUE]]), colnames(object))
  }
  controls <- backend_control
  controls$runner <- NULL
  output_dir <- controls$output_dir %||% tempfile("shennong-copykat-")
  controls$output_dir <- NULL
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results <- lapply(unique(samples), function(sample) {
    cells <- names(samples)[samples == sample]
    normal <- intersect(cells, reference_cells)
    defaults <- list(
      rawmat = as.matrix(counts[, cells, drop = FALSE]), id.type = "S",
      sam.name = gsub("[^[:alnum:]_-]", "_", paste(store_name, sample, sep = "_")),
      norm.cell.names = normal, genome = if (tolower(genome %||% "human") %in% c("mouse", "mm10")) "mm10" else "hg20",
      n.cores = 1, plot.genes = FALSE, output.seg = FALSE
    )
    old <- setwd(output_dir)
    on.exit(setwd(old), add = TRUE)
    do.call(copykat::copykat, utils::modifyList(defaults, controls, keep.null = TRUE))
  })
  prediction <- dplyr::bind_rows(lapply(results, function(result) tibble::as_tibble(result$prediction)))
  cna <- lapply(results, function(result) result$CNAmat)
  chromosomes <- dplyr::bind_rows(lapply(cna, .sn_copykat_chromosomes, cells = colnames(object)))
  scores <- .sn_cnv_scores_from_chromosomes(chromosomes)
  prediction_cell <- .sn_communication_column(prediction, c("cell.names", "cell", "cell_name"))
  prediction_call <- .sn_communication_column(prediction, c("copykat.pred", "prediction", "call"))
  calls <- if (is_null(prediction_cell) || is_null(prediction_call)) NULL else {
    stats::setNames(as.character(prediction[[prediction_call]]), as.character(prediction[[prediction_cell]]))
  }
  list(
    object = object, scores = scores, chromosomes = chromosomes,
    subclone = NULL, prediction = calls, embedding = NULL,
    artifacts = list(output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE))
  )
}

.sn_cnv_sample_summary <- function(primary) {
  dplyr::bind_rows(lapply(split(primary, primary$sample), function(table) {
    tibble::tibble(
      sample = table$sample[[1]], n_cells = nrow(table),
      mean_cnv_score = mean(table$cnv_score, na.rm = TRUE),
      median_malignant_score = stats::median(table$malignant_score, na.rm = TRUE),
      malignant_fraction = mean(table$malignant_call == "malignant", na.rm = TRUE),
      reference_fraction = mean(table$malignant_call == "reference", na.rm = TRUE)
    )
  }))
}

#' Run copy-number and malignancy analysis
#'
#' Runs an optional CNV backend, standardizes per-cell and chromosome-level
#' evidence, derives a reference-calibrated malignancy score, assigns
#' subclones, and stores sample-aware diagnostics.
#'
#' @param object A Seurat object.
#' @param method CNV backend: `"infercnvpy"` or `"copykat"`.
#' @param reference_cells Normal reference cell names or a logical vector.
#' @param genome Species/genome label. Human and mouse are supported by the
#'   bundled inferCNVpy positions and CopyKAT adapter.
#' @param store_name Stored result name.
#' @param reference_by,reference_cat Alternative metadata-based reference definition.
#' @param sample_by Optional sample/patient metadata column.
#' @param assay,layer Expression assay and layer. CopyKAT should use counts;
#'   inferCNVpy should use normalized data.
#' @param malignant_threshold Reference-scaled CNV threshold for malignant calls.
#' @param subclones Maximum number of fallback hierarchical subclones.
#' @param association_features Number of variable genes tested for CNV-expression association.
#' @param backend_control Named list forwarded to the selected backend. A
#'   `runner` function can provide a custom backend adapter.
#' @param return_object Return the modified object or the unified result.
#'
#' @return A Seurat object or unified CNV result.
#' @examples
#' \dontrun{
#' object <- sn_run_cnv(
#'   object, method = "infercnvpy", reference_by = "cell_type",
#'   reference_cat = c("T cell", "Myeloid"), sample_by = "patient"
#' )
#' cnv <- sn_get_result(object, "cnv", "cnv")
#' cnv$tables$sample_summary
#' }
#' @export
sn_run_cnv <- function(object,
                       method = c("infercnvpy", "copykat"),
                       reference_cells = NULL,
                       genome = NULL,
                       store_name = "cnv",
                       reference_by = NULL,
                       reference_cat = NULL,
                       sample_by = NULL,
                       assay = NULL,
                       layer = NULL,
                       malignant_threshold = 2,
                       subclones = 2L,
                       association_features = 50L,
                       backend_control = list(),
                       return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  if (!is_null(sample_by) && !sample_by %in% colnames(object[[]])) {
    stop("`sample_by` column '", sample_by, "' was not found.", call. = FALSE)
  }
  reference_cells <- .sn_cnv_reference_cells(object, reference_cells, reference_by, reference_cat)
  genome <- genome %||% tryCatch(sn_get_species(object), error = function(e) "human")
  resolved_assay <- assay %||% SeuratObject::DefaultAssay(object)
  association_layer <- layer %||% if (identical(method, "copykat")) {
    "counts"
  } else {
    .sn_select_infercnvpy_layer(object, resolved_assay)
  }
  backend <- switch(
    method,
    infercnvpy = .sn_run_cnv_infercnvpy,
    copykat = .sn_run_cnv_copykat
  )(
    object = object, reference_cells = reference_cells, genome = genome,
    store_name = store_name, assay = assay, layer = layer,
    sample_by = sample_by, backend_control = backend_control
  )
  object <- backend$object %||% object
  scores <- backend$scores
  if (is.null(names(scores))) stop("CNV backend scores must be named by cell.", call. = FALSE)
  scores <- scores[colnames(object)]
  malignancy <- .sn_cnv_malignancy(scores, reference_cells, threshold = malignant_threshold)
  if (!is_null(backend$prediction)) {
    predictions <- tolower(as.character(backend$prediction[colnames(object)]))
    malignancy$call[predictions == "aneuploid"] <- "malignant"
    malignancy$call[predictions == "diploid"] <- "non_malignant"
    malignancy$call[colnames(object) %in% reference_cells] <- "reference"
  }
  subclone <- .sn_cnv_subclones(
    backend$chromosomes %||% tibble::tibble(), malignancy$call,
    existing = backend$subclone, k = subclones
  )
  metadata <- object[[]]
  sample <- if (is_null(sample_by)) rep("all", ncol(object)) else as.character(metadata[[sample_by]])
  primary <- tibble::tibble(
    cell = colnames(object), sample = sample,
    cnv_score = as.numeric(scores[colnames(object)]),
    malignant_score = as.numeric(malignancy$score[colnames(object)]),
    malignant_call = unname(malignancy$call[colnames(object)]),
    subclone = unname(subclone[colnames(object)]), method = method,
    is_reference = colnames(object) %in% reference_cells
  )
  prefix <- gsub("[^[:alnum:]_]+", "_", store_name)
  cell_metadata <- data.frame(
    stats::setNames(list(primary$cnv_score), paste0(prefix, "_cnv_score")),
    stats::setNames(list(primary$malignant_score), paste0(prefix, "_malignant_score")),
    stats::setNames(list(primary$malignant_call), paste0(prefix, "_malignant_call")),
    stats::setNames(list(primary$subclone), paste0(prefix, "_subclone")),
    row.names = primary$cell, check.names = FALSE
  )
  object <- SeuratObject::AddMetaData(object, metadata = cell_metadata)
  associations <- .sn_cnv_expression_association(
    object, stats::setNames(primary$malignant_score, primary$cell),
    assay = resolved_assay, layer = association_layer, n_features = association_features
  )
  embeddings <- list()
  if (!is_null(backend$embedding)) embeddings$cnv_umap <- backend$embedding
  result <- list(
    schema_version = "1.0", analysis_type = "cnv", name = store_name,
    method = method, backend = method,
    input = list(
      assay = resolved_assay, layer = association_layer,
      genome = genome, reference_cells = reference_cells, sample_by = sample_by
    ),
    parameters = list(malignant_threshold = malignant_threshold, subclones = subclones),
    tables = list(
      primary = primary,
      chromosome = backend$chromosomes %||% tibble::tibble(),
      sample_summary = .sn_cnv_sample_summary(primary),
      expression_association = associations
    ),
    embeddings = embeddings, graphs = list(),
    models = list(backend_artifacts = backend$artifacts %||% list()),
    diagnostics = list(
      reference_cells = length(reference_cells), finite_scores = sum(is.finite(primary$cnv_score)),
      reference_center = malignancy$center, reference_spread = malignancy$spread,
      malignant_cells = sum(primary$malignant_call == "malignant", na.rm = TRUE),
      samples = length(unique(primary$sample))
    ),
    warnings = character(), provenance = .sn_analysis_provenance()
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "cnv", store_name, result)
  object <- .sn_log_seurat_command(object, assay = resolved_assay, name = "sn_run_cnv")
  if (isTRUE(return_object)) object else sn_get_result(object, "cnv", store_name)
}
