#' Run bulk RNA-seq deconvolution with single-cell references
#'
#' This helper supports two bulk deconvolution workflows that combine a
#' single-cell reference with bulk RNA-seq mixtures:
#'
#' \itemize{
#'   \item \code{"bayesprism"} runs the \pkg{BayesPrism} R package locally.
#'   \item \code{"cibersortx"} prepares inputs and runs the local CIBERSORTx
#'         container workflow, or imports an existing fractions result file.
#' }
#'
#' @param x A \code{Seurat} object used as the single-cell reference, or a
#'   gene-by-cell matrix-like object.
#' @param bulk A bulk expression matrix-like object.
#' @param method One of \code{"bayesprism"} or \code{"cibersortx"}.
#' @param cell_type_by Metadata column containing cell-type labels when
#'   \code{x} is a \code{Seurat} object. If \code{x} is a matrix, supply
#'   \code{cell_type_labels} instead.
#' @param cell_state_by Optional metadata column containing cell-state labels.
#'   Defaults to \code{cell_type_by}.
#' @param cell_type_labels Optional vector of cell-type labels for matrix
#'   references.
#' @param cell_state_labels Optional vector of cell-state labels for matrix
#'   references.
#' @param assay Assay used to extract the single-cell reference matrix.
#' @param layer Layer used to extract the single-cell reference matrix.
#' @param bulk_gene_axis Orientation of genes in \code{bulk}. Use
#'   \code{"rows"} for the common gene-by-sample layout, \code{"columns"} for
#'   sample-by-gene input, or \code{"auto"} to infer it from the overlap with
#'   the reference genes.
#' @param key Optional malignant-cell label_by passed to BayesPrism.
#' @param outdir Output directory used by CIBERSORTx file export.
#' @param prefix Prefix used for exported files and stored result names.
#' @param result_id Stable identifier for the stored deconvolution result
#'   when \code{x} is a \code{Seurat} object and a fraction table is available.
#' @param cibersortx_result Optional path to a completed CIBERSORTx fractions
#'   result file to import instead of running the local container.
#' @param cibersortx_email Optional CIBERSORTx account email.
#' @param cibersortx_token Optional CIBERSORTx access token.
#' @param cibersortx_container Container runtime used for local execution. One
#'   of \code{"docker"} or \code{"apptainer"}.
#' @param cibersortx_container_path Optional Apptainer image path.
#' @param cibersortx_dry_run If \code{TRUE}, prepare files and return local
#'   commands without running the container.
#' @param cibersortx_rmbatch_b_mode Whether to enable B-mode batch_by correction.
#' @param cibersortx_rmbatch_s_mode Whether to enable S-mode batch_by correction.
#' @param cibersortx_perm Number of permutations used by the fractions module.
#' @param cibersortx_qn Whether to enable quantile normalization.
#' @param cibersortx_absolute Whether to enable absolute mode.
#' @param cibersortx_abs_method CIBERSORTx absolute-mode method.
#' @param cibersortx_k_max Maximum condition number used when constructing the
#'   signature matrix.
#' @param gibbs_control Optional BayesPrism Gibbs-sampler control list.
#' @param opt_control Optional BayesPrism optimization control list.
#' @param n_cores Number of cores passed to BayesPrism.
#' @param update_gibbs Whether BayesPrism should run the final Gibbs update.
#' @param max_dense_gb Maximum estimated size, in GiB, of any expression matrix
#'   that may be materialized as a dense double matrix. Large sparse inputs fail
#'   before allocation unless the caller explicitly raises this budget.
#' @param return_object If \code{TRUE} and \code{x} is a \code{Seurat} object,
#'   return the updated object when a result table is available.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#'
#' @return A stored-result list, an export-bundle list, or an updated
#'   \code{Seurat} object depending on the selected backend and
#'   \code{return_object}.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(20 * 18, lambda = 3), nrow = 20, ncol = 18)
#'   rownames(counts) <- paste0("gene", seq_len(20))
#'   colnames(counts) <- paste0("cell", seq_len(18))
#'   ref <- sn_initialize_seurat_object(counts, species = "human")
#'   ref$cell_type <- rep(c("Tcell", "Bcell", "Mono"), each = 6)
#'   bulk <- cbind(
#'     sample_a = rowSums(counts[, 1:12, drop = FALSE]),
#'     sample_b = rowSums(counts[, 7:18, drop = FALSE])
#'   )
#'   bundle <- sn_run_bulk_deconvolution(
#'     ref,
#'     bulk = bulk,
#'     method = "cibersortx",
#'     cell_type_by = "cell_type",
#'     outdir = tempdir(),
#'     cibersortx_email = "demo@example.org",
#'     cibersortx_token = "fake-token",
#'     cibersortx_dry_run = TRUE,
#'     return_object = FALSE
#'   )
#'   names(bundle$files)
#' }
#'
#' @export
sn_run_bulk_deconvolution <- function(x,
                               bulk,
                               method = c("bayesprism", "cibersortx"),
                               cell_type_by = NULL,
                               cell_state_by = NULL,
                               cell_type_labels = NULL,
                               cell_state_labels = NULL,
                               assay = "RNA",
                               layer = "counts",
                               bulk_gene_axis = c("auto", "rows", "columns"),
                               key = NULL,
                               outdir = NULL,
                               prefix = "deconvolution",
                               result_id = "default",
                               cibersortx_result = NULL,
                               cibersortx_email = NULL,
                               cibersortx_token = NULL,
                               cibersortx_container = c("docker", "apptainer"),
                               cibersortx_container_path = NULL,
                               cibersortx_dry_run = FALSE,
                               cibersortx_rmbatch_b_mode = FALSE,
                               cibersortx_rmbatch_s_mode = FALSE,
                               cibersortx_perm = 0,
                               cibersortx_qn = FALSE,
                               cibersortx_absolute = FALSE,
                               cibersortx_abs_method = "sig.score",
                               cibersortx_k_max = 999,
                               gibbs_control = list(),
                               opt_control = list(),
                               n_cores = 1,
                               update_gibbs = TRUE,
                               max_dense_gb = 2,
                               return_object = TRUE,
                               object = NULL) {
  result_id <- .sn_validate_result_id(result_id)
  x <- .sn_resolve_object_alias(x, object, missing(x))
  method <- match.arg(method)
  bulk_gene_axis <- match.arg(bulk_gene_axis)
  cibersortx_container <- match.arg(cibersortx_container)
  if (!is.numeric(max_dense_gb) || length(max_dense_gb) != 1L ||
      !is.finite(max_dense_gb) || max_dense_gb <= 0) {
    stop("`max_dense_gb` must be one positive finite number.", call. = FALSE)
  }
  if (inherits(x, "Seurat") && identical(method, "bayesprism") &&
      !.sn_name_declares_count_scale(layer)) {
    stop(
      "BayesPrism requires a raw/count-like reference layer; `", layer,
      "` is not declared count-scale.",
      call. = FALSE
    )
  }
  if (inherits(x, "Seurat") && identical(method, "cibersortx") &&
      .sn_name_declares_log_expression_scale(layer)) {
    stop(
      "CIBERSORTx accepts raw counts or non-log linear expression, but Seurat layer `",
      layer, "` is declared log/normalized. Supply a count-like or explicitly linear layer.",
      call. = FALSE
    )
  }

  reference_info <- .sn_prepare_deconvolution_reference(
    x = x,
    cell_type_col = cell_type_by,
    cell_state_col = cell_state_by,
    cell_type_labels = cell_type_labels,
    cell_state_labels = cell_state_labels,
    assay = assay,
    layer = layer,
    max_dense_gb = max_dense_gb
  )
  bulk_samples_by_gene <- .sn_prepare_bulk_matrix(
    bulk = bulk,
    reference_genes = colnames(reference_info$reference_cells_by_gene),
    gene_axis = bulk_gene_axis,
    max_dense_gb = max_dense_gb
  )

  if (identical(method, "bayesprism")) {
    result <- .sn_run_bayesprism(
      reference_cells_by_gene = reference_info$reference_cells_by_gene,
      cell_type_labels = reference_info$cell_type_labels,
      cell_state_labels = reference_info$cell_state_labels,
      bulk_samples_by_gene = bulk_samples_by_gene,
      key = key,
      gibbs_control = gibbs_control,
      opt_control = opt_control,
      n_cores = n_cores,
      update_gibbs = update_gibbs,
      max_dense_gb = max_dense_gb
    )
  } else {
    result <- .sn_run_cibersortx_local(
      reference_genes_by_cells = if (inherits(reference_info$reference_cells_by_gene, "Matrix")) {
        Matrix::t(reference_info$reference_cells_by_gene)
      } else {
        t(reference_info$reference_cells_by_gene)
      },
      cell_type_labels = reference_info$cell_type_labels,
      bulk_genes_by_samples = t(bulk_samples_by_gene),
      outdir = outdir,
      prefix = prefix,
      result_path = cibersortx_result,
      email = cibersortx_email,
      token = cibersortx_token,
      container = cibersortx_container,
      container_path = cibersortx_container_path,
      dry_run = cibersortx_dry_run,
      rmbatch_b_mode = cibersortx_rmbatch_b_mode,
      rmbatch_s_mode = cibersortx_rmbatch_s_mode,
      perm = cibersortx_perm,
      qn = cibersortx_qn,
      absolute = cibersortx_absolute,
      abs_method = cibersortx_abs_method,
      k_max = cibersortx_k_max,
      max_dense_gb = max_dense_gb
    )
  }

  if (is.list(result$scale_provenance)) {
    result$scale_provenance$reference_source <- if (inherits(x, "Seurat")) "Seurat" else "matrix"
    result$scale_provenance$reference_assay <- if (inherits(x, "Seurat")) assay else NULL
    result$scale_provenance$reference_layer <- if (inherits(x, "Seurat")) layer else NULL
    result$scale_provenance$bulk_gene_axis <- bulk_gene_axis
    result$artifacts <- result$artifacts %||% list()
    result$artifacts$scale_provenance <- result$scale_provenance
  }

  if (!inherits(x, "Seurat") || is.null(result$table)) {
    return(result)
  }

  object <- sn_store_deconvolution(
    object = x,
    result = result$table,
    result_id = result_id,
    method = method,
    bulk_samples = rownames(bulk_samples_by_gene),
    reference_label = reference_info$reference_label,
    artifacts = result$artifacts %||% result$files
  )

  if (isTRUE(return_object)) {
    return(.sn_log_seurat_command(object = object, name = "sn_run_bulk_deconvolution"))
  }

  sn_get_result(
    object = object,
    type = "deconvolution",
    result_id = result_id
  )
}

.sn_cibersortx_env <- new.env(parent = emptyenv())

#' Store local CIBERSORTx credentials for container execution
#'
#' @param email Email registered with CIBERSORTx.
#' @param token Access token issued by CIBERSORTx.
#'
#' @return Invisibly returns \code{TRUE}.
#' @export
sn_set_cibersortx_credentials <- function(email, token) {
  assign("email", email, envir = .sn_cibersortx_env)
  assign("token", token, envir = .sn_cibersortx_env)
  invisible(TRUE)
}

.sn_prepare_deconvolution_reference <- function(x,
                                                cell_type_col = NULL,
                                                cell_state_col = NULL,
                                                cell_type_labels = NULL,
                                                cell_state_labels = NULL,
                                                assay = "RNA",
                                                layer = "counts",
                                                max_dense_gb = 2) {
  if (inherits(x, "Seurat")) {
    if (is.null(cell_type_col) || !cell_type_col %in% colnames(x[[]])) {
      stop("`cell_type_col` must be supplied and present in `x@meta.data`.", call. = FALSE)
    }
    if (!is.null(cell_state_col) && !cell_state_col %in% colnames(x[[]])) {
      stop(glue("Metadata column '{cell_state_col}' was not found."), call. = FALSE)
    }

    counts <- .sn_get_seurat_layer_data(object = x, assay = assay, layer = layer)
    # Keep sparse/BPCells-backed inputs sparse for streamable backends such as
    # CIBERSORTx. Backends that require a dense matrix enforce their own budget
    # immediately before conversion.
    reference <- if (inherits(counts, "Matrix")) Matrix::t(counts) else t(counts)
    cell_type_labels <- as.character(x[[cell_type_col, drop = TRUE]])
    cell_state_labels <- if (is.null(cell_state_col)) {
      cell_type_labels
    } else {
      as.character(x[[cell_state_col, drop = TRUE]])
    }
    reference_label <- cell_type_col
  } else {
    reference <- if (inherits(x, "Matrix") || is.matrix(x)) {
      x
    } else if (is.data.frame(x)) {
      as.matrix(x)
    } else {
      stop("`x` must be matrix-like.", call. = FALSE)
    }
    reference_values <- if (inherits(reference, "Matrix") &&
        "x" %in% methods::slotNames(reference)) {
      methods::slot(reference, "x")
    } else {
      reference
    }
    if (!is.numeric(reference_values)) {
      stop("`x` must contain numeric expression values.", call. = FALSE)
    }
    if (is.null(cell_type_labels)) {
      stop("`cell_type_labels` must be supplied when `x` is not a Seurat object.", call. = FALSE)
    }
    cell_type_labels <- as.character(cell_type_labels)
    cell_state_labels <- as.character(cell_state_labels %||% cell_type_labels)
    if (length(cell_type_labels) != nrow(reference)) {
      stop("`cell_type_labels` must have length equal to the number of reference cells.", call. = FALSE)
    }
    if (length(cell_state_labels) != nrow(reference)) {
      stop("`cell_state_labels` must have length equal to the number of reference cells.", call. = FALSE)
    }
    reference_label <- "cell_type"
  }

  if (is.null(rownames(reference))) {
    rownames(reference) <- paste0("cell_", seq_len(nrow(reference)))
  }
  if (is.null(colnames(reference))) {
    stop("Reference genes must be supplied as column names.", call. = FALSE)
  }

  if (length(cell_type_labels) != nrow(reference) ||
      length(cell_state_labels) != nrow(reference)) {
    stop("Cell-type and cell-state labels must have one value per reference cell.", call. = FALSE)
  }
  cell_type_labels <- trimws(as.character(cell_type_labels))
  cell_state_labels <- trimws(as.character(cell_state_labels))
  if (anyNA(cell_type_labels) || any(!nzchar(cell_type_labels)) ||
      anyNA(cell_state_labels) || any(!nzchar(cell_state_labels))) {
    stop("Cell-type and cell-state labels must be non-missing and non-empty.", call. = FALSE)
  }

  list(
    reference_cells_by_gene = reference,
    cell_type_labels = cell_type_labels,
    cell_state_labels = cell_state_labels,
    reference_label = reference_label
  )
}

.sn_as_matrix <- function(x, name = "input", max_dense_gb = 2) {
  .sn_assert_dense_materialization_budget(
    x, max_dense_gb = max_dense_gb, name = name
  )
  if (inherits(x, "Matrix")) {
    return(as.matrix(x))
  }
  if (inherits(x, "data.frame")) {
    return(as.matrix(x))
  }
  if (inherits(x, "matrix")) {
    return(x)
  }
  stop(glue("`{name}` must be matrix-like."), call. = FALSE)
}

.sn_prepare_bulk_matrix <- function(bulk,
                                    reference_genes,
                                    gene_axis = c("auto", "rows", "columns"),
                                    max_dense_gb = 2) {
  gene_axis <- match.arg(gene_axis)
  bulk <- .sn_as_matrix(bulk, name = "bulk", max_dense_gb = max_dense_gb)
  if (is.null(rownames(bulk)) && is.null(colnames(bulk))) {
    stop("Bulk input must have gene names either on rows or on columns.", call. = FALSE)
  }

  if (identical(gene_axis, "auto")) {
    row_overlap <- if (is.null(rownames(bulk))) 0L else sum(rownames(bulk) %in% reference_genes)
    col_overlap <- if (is.null(colnames(bulk))) 0L else sum(colnames(bulk) %in% reference_genes)
    gene_axis <- if (row_overlap >= col_overlap) "rows" else "columns"
  }

  bulk_samples_by_gene <- if (identical(gene_axis, "rows")) {
    if (is.null(rownames(bulk))) {
      stop("Bulk genes were requested on rows but row names are missing.", call. = FALSE)
    }
    t(bulk)
  } else {
    if (is.null(colnames(bulk))) {
      stop("Bulk genes were requested on columns but column names are missing.", call. = FALSE)
    }
    bulk
  }

  if (is.null(rownames(bulk_samples_by_gene))) {
    rownames(bulk_samples_by_gene) <- paste0("sample_", seq_len(nrow(bulk_samples_by_gene)))
  }
  if (is.null(colnames(bulk_samples_by_gene))) {
    stop("Bulk genes could not be identified after orientation.", call. = FALSE)
  }

  bulk_samples_by_gene
}

.sn_name_declares_log_expression_scale <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    return(FALSE)
  }
  key <- tolower(name)
  grepl("(^|[._-])(data|log|scaled?|normalized?)($|[._-])", key)
}

.sn_deconvolution_matrix_values <- function(x, label) {
  if (!(inherits(x, "Matrix") || is.matrix(x))) {
    stop(label, " must be a numeric matrix-like object.", call. = FALSE)
  }
  values <- if (inherits(x, "Matrix") && "x" %in% methods::slotNames(x)) {
    methods::slot(x, "x")
  } else {
    x
  }
  if (!is.numeric(values)) {
    stop(label, " must contain numeric expression values.", call. = FALSE)
  }
  as.numeric(values)
}

.sn_validate_deconvolution_expression <- function(x,
                                                  gene_axis = c("rows", "columns"),
                                                  label,
                                                  require_counts = FALSE) {
  gene_axis <- match.arg(gene_axis)
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(label, " must contain at least one gene and one observation.", call. = FALSE)
  }
  gene_ids <- if (identical(gene_axis, "rows")) rownames(x) else colnames(x)
  normalized_ids <- trimws(as.character(gene_ids %||% character()))
  expected_genes <- if (identical(gene_axis, "rows")) nrow(x) else ncol(x)
  if (length(normalized_ids) != expected_genes ||
      anyNA(normalized_ids) || any(!nzchar(normalized_ids)) || anyDuplicated(normalized_ids)) {
    stop(label, " gene identifiers must be unique and non-empty.", call. = FALSE)
  }
  values <- .sn_deconvolution_matrix_values(x, label)
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
    stop(label, " must contain finite, non-negative expression values.", call. = FALSE)
  }
  integer_like <- length(values) == 0L || all(abs(values - round(values)) <= 1e-8)
  if (isTRUE(require_counts) && !integer_like) {
    stop(label, " must contain integer-like raw counts; values are never rounded silently.", call. = FALSE)
  }
  if (integer_like) "counts" else "linear_nonlog"
}

.sn_validate_deconvolution_labels <- function(labels, n_cells, label) {
  labels <- trimws(as.character(labels))
  if (length(labels) != n_cells || anyNA(labels) || any(!nzchar(labels))) {
    stop(label, " must contain one non-missing, non-empty value per reference cell.", call. = FALSE)
  }
  labels
}

.sn_deconvolution_scale_provenance <- function(reference_scale, bulk_scale) {
  list(
    reference = reference_scale,
    bulk = bulk_scale,
    matched = identical(reference_scale, bulk_scale),
    contract = if (identical(reference_scale, "counts")) {
      "finite_nonnegative_integer_like_raw_counts"
    } else {
      "finite_nonnegative_nonlog_linear_expression"
    }
  )
}

.sn_fraction_matrix_to_long <- function(fraction_matrix,
                                        sample_by = "sample",
                                        cell_type_by = "cell_type",
                                        value_col = "fraction") {
  fraction_matrix <- as.matrix(fraction_matrix)
  sample_names <- rownames(fraction_matrix) %||% paste0("sample_", seq_len(nrow(fraction_matrix)))
  cell_types <- colnames(fraction_matrix) %||% paste0("cell_type_", seq_len(ncol(fraction_matrix)))

  table <- expand.grid(
    sample = sample_names,
    cell_type = cell_types,
    stringsAsFactors = FALSE
  )
  table$fraction <- as.numeric(fraction_matrix[cbind(match(table$sample, sample_names), match(table$cell_type, cell_types))])
  names(table) <- c(sample_by, cell_type_by, value_col)
  tibble::as_tibble(table)
}

.sn_restore_directory <- function(path) {
  if (!dir.exists(path)) {
    created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(created) && !dir.exists(path)) {
      stop("Could not restore the directory removed by an external backend: ", path, call. = FALSE)
    }
  }
  invisible(path)
}

.sn_with_preserved_directory <- function(path, code) {
  on.exit(.sn_restore_directory(path), add = TRUE)
  value <- force(code)
  .sn_restore_directory(path)
  value
}

.sn_run_bayesprism <- function(reference_cells_by_gene,
                               cell_type_labels,
                               cell_state_labels,
                               bulk_samples_by_gene,
                               key = NULL,
                               gibbs_control = list(),
                               opt_control = list(),
                               n_cores = 1,
                               update_gibbs = TRUE,
                               max_dense_gb = 2) {
  reference_scale <- .sn_validate_deconvolution_expression(
    reference_cells_by_gene,
    gene_axis = "columns",
    label = "BayesPrism reference",
    require_counts = TRUE
  )
  bulk_scale <- .sn_validate_deconvolution_expression(
    bulk_samples_by_gene,
    gene_axis = "columns",
    label = "BayesPrism mixture",
    require_counts = TRUE
  )
  cell_type_labels <- .sn_validate_deconvolution_labels(
    cell_type_labels, nrow(reference_cells_by_gene), "BayesPrism cell-type labels"
  )
  cell_state_labels <- .sn_validate_deconvolution_labels(
    cell_state_labels, nrow(reference_cells_by_gene), "BayesPrism cell-state labels"
  )
  scale_provenance <- .sn_deconvolution_scale_provenance(reference_scale, bulk_scale)
  check_installed_github(
    pkg = "BayesPrism",
    repo = "Danko-Lab/BayesPrism/BayesPrism",
    reason = "to run BayesPrism bulk deconvolution."
  )
  .sn_assert_dense_materialization_budget(
    reference_cells_by_gene, max_dense_gb, "BayesPrism reference"
  )
  .sn_assert_dense_materialization_budget(
    bulk_samples_by_gene, max_dense_gb, "BayesPrism mixture"
  )
  reference_cells_by_gene <- .sn_as_matrix(
    reference_cells_by_gene,
    name = "BayesPrism reference",
    max_dense_gb = max_dense_gb
  )
  bulk_samples_by_gene <- .sn_as_matrix(
    bulk_samples_by_gene,
    name = "BayesPrism mixture",
    max_dense_gb = max_dense_gb
  )

  prism <- BayesPrism::new.prism(
    reference = reference_cells_by_gene,
    input.type = "count.matrix",
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_state_labels,
    key = key,
    mixture = bulk_samples_by_gene
  )
  session_tempdir <- tempdir(check = TRUE)
  bp_fit <- .sn_with_preserved_directory(
    session_tempdir,
    BayesPrism::run.prism(
      prism = prism,
      n.cores = n_cores,
      update.gibbs = update_gibbs,
      gibbs.control = gibbs_control,
      opt.control = opt_control
    )
  )

  which_theta <- if (isTRUE(update_gibbs)) "final" else "first"
  fractions <- BayesPrism::get.fraction(
    bp = bp_fit,
    which.theta = which_theta,
    state.or.type = "type"
  )

  list(
    table = .sn_fraction_matrix_to_long(fractions),
    artifacts = list(
      fraction_matrix = fractions,
      fit = bp_fit,
      prism = prism,
      scale_provenance = scale_provenance
    ),
    scale_provenance = scale_provenance,
    method = "bayesprism"
  )
}

.sn_get_cibersortx_credentials <- function(email = NULL, token = NULL) {
  email <- email %||% if (exists("email", envir = .sn_cibersortx_env, inherits = FALSE)) get("email", envir = .sn_cibersortx_env) else NULL
  token <- token %||% if (exists("token", envir = .sn_cibersortx_env, inherits = FALSE)) get("token", envir = .sn_cibersortx_env) else NULL
  email <- email %||% Sys.getenv("SHENNONG_CIBERSORTX_EMAIL", unset = "")
  token <- token %||% Sys.getenv("SHENNONG_CIBERSORTX_TOKEN", unset = "")
  if (!nzchar(email) || !nzchar(token)) {
    stop(
      "CIBERSORTx credentials are required. Supply `cibersortx_email` and `cibersortx_token`, ",
      "call `sn_set_cibersortx_credentials()`, or set SHENNONG_CIBERSORTX_EMAIL / SHENNONG_CIBERSORTX_TOKEN.",
      call. = FALSE
    )
  }
  if (length(email) != 1L || length(token) != 1L ||
      grepl("[\r\n]", email) || grepl("[\r\n]", token)) {
    stop("CIBERSORTx credentials must be scalar values without line breaks.", call. = FALSE)
  }
  list(email = email, token = token)
}

.sn_write_cibersortx_credentials <- function(email, token) {
  path <- tempfile("shennong-cibersortx-credentials-")
  writeLines(c(email, token), con = path, useBytes = TRUE)
  Sys.chmod(path, mode = "0600")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.sn_check_container_binary <- function(container = c("docker", "apptainer")) {
  container <- match.arg(container)
  binary <- Sys.which(container)
  if (!nzchar(binary)) {
    stop(glue("Could not find the `{container}` executable required for local CIBERSORTx runs."), call. = FALSE)
  }
  binary
}

.sn_cibersortx_display_command <- function(executable,
                                           args,
                                           sensitive_flags = c("--username", "--token")) {
  display_args <- as.character(args)
  for (flag in sensitive_flags) {
    idx <- which(display_args == flag)
    value_idx <- idx + 1L
    value_idx <- value_idx[value_idx <= length(display_args)]
    display_args[value_idx] <- "<redacted>"
  }
  paste(c(shQuote(executable), shQuote(display_args)), collapse = " ")
}

.sn_cibersortx_command_spec <- function(executable, args) {
  list(
    executable = executable,
    args = as.character(args),
    display = .sn_cibersortx_display_command(executable, args)
  )
}

.sn_run_command <- function(command, dry_run = FALSE, verbose = FALSE) {
  if (is.list(command) && !is.null(command$executable) && !is.null(command$args)) {
    display <- command$display %||% .sn_cibersortx_display_command(command$executable, command$args)
    if (isTRUE(dry_run)) {
      return(list(status = 0L, command = display))
    }

    output_target <- if (isTRUE(verbose)) "" else FALSE
    status <- system2(command = command$executable, args = command$args, stdout = output_target, stderr = output_target)
    return(list(status = status, command = display))
  }

  stop("Internal command execution requires a structured command spec.", call. = FALSE)
}

.sn_transform_and_save_cibersortx_single_cell <- function(reference_genes_by_cells,
                                                          cell_type_labels,
                                                          path,
                                                          max_dense_gb = 2) {
  output_file <- file.path(path, "sample_file_for_cibersort.txt")
  if (length(cell_type_labels) != ncol(reference_genes_by_cells) ||
      anyNA(cell_type_labels) || any(!nzchar(as.character(cell_type_labels)))) {
    stop("CIBERSORTx cell-type labels must be non-empty and match the reference columns.", call. = FALSE)
  }
  .sn_write_cibersortx_matrix(
    x = reference_genes_by_cells,
    output_file = output_file,
    first_column = "GeneSymbol",
    column_labels = as.character(cell_type_labels),
    max_dense_gb = max_dense_gb,
    label = "CIBERSORTx single-cell reference"
  )
  output_file
}

.sn_transform_and_save_cibersortx_bulk <- function(bulk_genes_by_samples, path,
                                                   max_dense_gb = 2) {
  output_file <- file.path(path, "mixture_file_for_cibersort.txt")
  if (is.null(colnames(bulk_genes_by_samples)) ||
      anyNA(colnames(bulk_genes_by_samples)) ||
      any(!nzchar(colnames(bulk_genes_by_samples)))) {
    stop("CIBERSORTx bulk samples require non-empty column names.", call. = FALSE)
  }
  .sn_write_cibersortx_matrix(
    x = bulk_genes_by_samples,
    output_file = output_file,
    first_column = "Gene",
    column_labels = colnames(bulk_genes_by_samples),
    max_dense_gb = max_dense_gb,
    label = "CIBERSORTx bulk mixture"
  )
  output_file
}

.sn_write_cibersortx_matrix <- function(x,
                                        output_file,
                                        first_column,
                                        column_labels,
                                        max_dense_gb = 2,
                                        label = "CIBERSORTx matrix") {
  if (!(inherits(x, "Matrix") || is.matrix(x))) {
    stop("`x` must be a numeric matrix-like object.", call. = FALSE)
  }
  matrix_values <- if (inherits(x, "Matrix") && "x" %in% methods::slotNames(x)) {
    methods::slot(x, "x")
  } else {
    x
  }
  if (!is.numeric(matrix_values)) {
    stop(label, " must contain numeric values.", call. = FALSE)
  }
  if (is.null(rownames(x)) || anyNA(rownames(x)) || any(!nzchar(rownames(x))) ||
      anyDuplicated(rownames(x))) {
    stop(label, " requires unique, non-empty gene row names.", call. = FALSE)
  }
  values <- if (inherits(x, "Matrix") && "x" %in% methods::slotNames(x)) {
    methods::slot(x, "x")
  } else {
    as.numeric(x)
  }
  if (any(!is.finite(values)) || any(values < 0)) {
    stop(label, " must contain finite, non-negative expression values.", call. = FALSE)
  }
  if (!is.numeric(max_dense_gb) || length(max_dense_gb) != 1L ||
      is.na(max_dense_gb) || !is.finite(max_dense_gb) || max_dense_gb <= 0) {
    stop("`max_dense_gb` must be one positive finite number.", call. = FALSE)
  }
  if (length(column_labels) != ncol(x)) {
    stop("CIBERSORTx column labels must match the number of matrix columns.", call. = FALSE)
  }

  # write.table() formats numeric cells as strings. Budget for that expensive
  # representation and cap chunks so the complete matrix is never coerced to a
  # giant character object (the old rbind/data.frame path did exactly that).
  budget_bytes <- max_dense_gb * 1024^3
  estimated_row_bytes <- max(1, ncol(x)) * 128
  chunk_rows <- min(256L, floor(budget_bytes / estimated_row_bytes))
  if (!is.finite(chunk_rows) || chunk_rows < 1L) {
    stop(
      "Streaming one row of ", label,
      " is estimated to exceed `max_dense_gb = ", max_dense_gb, "`.",
      call. = FALSE
    )
  }
  chunk_rows <- as.integer(min(chunk_rows, max(1L, nrow(x))))
  connection <- file(output_file, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  writeLines(
    paste(c(first_column, as.character(column_labels)), collapse = "\t"),
    connection,
    useBytes = TRUE
  )
  if (nrow(x) == 0L) return(invisible(output_file))
  starts <- seq.int(1L, nrow(x), by = chunk_rows)
  for (start in starts) {
    stop <- min(nrow(x), start + chunk_rows - 1L)
    chunk <- as.matrix(x[start:stop, , drop = FALSE])
    output <- data.frame(
      gene = rownames(chunk),
      chunk,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    utils::write.table(
      output,
      file = connection,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
  invisible(output_file)
}

.sn_create_cibersortx_command <- function(input_dir,
                                          output_dir,
                                          credentials_path,
                                          container = c("docker", "apptainer"),
                                          container_path = NULL,
                                          method = c("create_sig", "impute_cell_fractions"),
                                          verbose = FALSE,
                                          refsample = "sample_file_for_cibersort.txt",
                                          sigmatrix = "signature_matrix.txt",
                                          mixture = "mixture_file_for_cibersort.txt",
                                          label_by = "shennong",
                                          rmbatch_b_mode = FALSE,
                                          rmbatch_s_mode = FALSE,
                                          perm = 0,
                                          qn = FALSE,
                                          absolute = FALSE,
                                          abs_method = "sig.score",
                                          k_max = 999) {
  container <- match.arg(container)
  method <- match.arg(method)
  if (!is.character(credentials_path) || length(credentials_path) != 1L ||
      !file.exists(credentials_path)) {
    stop("`credentials_path` must identify the ephemeral CIBERSORTx credential file.", call. = FALSE)
  }
  secret_path <- "/run/secrets/shennong_cibersortx"
  credential_script <- paste(
    "IFS= read -r CIBERSORTX_USERNAME <", secret_path, ";",
    "CIBERSORTX_TOKEN=$(sed -n '2p'", secret_path, ");",
    "exec /src/CIBERSORTxFractions --single_cell TRUE",
    "--username \"$CIBERSORTX_USERNAME\" --token \"$CIBERSORTX_TOKEN\" \"$@\""
  )

  if (identical(container, "docker")) {
    executable <- "docker"
    base_args <- c(
      "run",
      "--mount", paste0(
        "type=bind,source=", normalizePath(credentials_path),
        ",target=", secret_path, ",readonly"
      ),
      "-v", paste0(normalizePath(input_dir), ":/src/data:z"),
      "-v", paste0(normalizePath(output_dir), ":/src/outdir:z"),
      "--entrypoint", "/bin/sh",
      "cibersortx/fractions",
      "-c", shQuote(credential_script), "shennong-cibersortx"
    )
  } else {
    if (is.null(container_path) || !file.exists(path.expand(container_path))) {
      stop("`cibersortx_container_path` must point to an existing Apptainer image.", call. = FALSE)
    }
    executable <- "apptainer"
    base_args <- c(
      "exec",
      "--no-home",
      "-c",
      "-B", paste0(normalizePath(credentials_path), ":", secret_path, ":ro"),
      "-B", paste0(normalizePath(input_dir), "/:/src/data"),
      "-B", paste0(normalizePath(output_dir), "/:/src/outdir"),
      normalizePath(path.expand(container_path)),
      "/bin/sh", "-c", shQuote(credential_script), "shennong-cibersortx"
    )
  }

  if (isTRUE(verbose)) {
    base_args <- c(base_args, "--verbose", "TRUE")
  }
  options <- if (identical(method, "create_sig")) {
    c(
      "--refsample", refsample,
      "--G.min", "300",
      "--G.max", "500",
      "--q.value", "0.01",
      "--filter", "FALSE",
      "--k.max", as.character(as.integer(k_max)),
      "--remake", "FALSE",
      "--replicates", "5",
      "--sampling", "0.5",
      "--fraction", "0.75"
    )
  } else {
    option_string <- c(
      "--mixture", mixture,
      "--sigmatrix", sigmatrix,
      "--perm", as.character(as.integer(perm)),
      "--label", label_by,
      "--rmbatchBmode", toupper(as.character(isTRUE(rmbatch_b_mode))),
      "--rmbatchSmode", toupper(as.character(isTRUE(rmbatch_s_mode))),
      "--sourceGEPs", sigmatrix,
      "--QN", toupper(as.character(isTRUE(qn))),
      "--absolute", toupper(as.character(isTRUE(absolute))),
      "--abs_method", abs_method
    )
    if (isTRUE(rmbatch_b_mode) || isTRUE(rmbatch_s_mode)) {
      option_string <- c(option_string, "--refsample", refsample)
    }
    option_string
  }

  .sn_cibersortx_command_spec(
    executable = executable,
    args = c(base_args, options)
  )
}

.sn_run_cibersortx_local <- function(reference_genes_by_cells,
                                     cell_type_labels,
                                     bulk_genes_by_samples,
                                     outdir = NULL,
                                     prefix,
                                     result_path = NULL,
                                     email = NULL,
                                     token = NULL,
                                     container = c("docker", "apptainer"),
                                     container_path = NULL,
                                     dry_run = FALSE,
                                     verbose = FALSE,
                                     rmbatch_b_mode = FALSE,
                                     rmbatch_s_mode = FALSE,
                                     perm = 0,
                                     qn = FALSE,
                                     absolute = FALSE,
                                     abs_method = "sig.score",
                                     k_max = 999,
                                     max_dense_gb = 2) {
  if (!is.null(result_path)) {
    if (!is.character(result_path) || length(result_path) != 1L ||
        is.na(result_path) || !nzchar(result_path)) {
      stop("`cibersortx_result` must be one existing fractions-result file.", call. = FALSE)
    }
    result_path <- path.expand(result_path)
    if (!file.exists(result_path) || dir.exists(result_path)) {
      stop("`cibersortx_result` does not exist or is not a file: ", result_path, call. = FALSE)
    }
    result_path <- normalizePath(result_path, winslash = "/", mustWork = TRUE)
    scale_provenance <- list(
      reference = "not_exported",
      bulk = "not_exported",
      matched = NA,
      contract = "import_only"
    )
    return(list(
      table = .sn_import_cibersortx_fractions(result_path = result_path),
      files = list(result = result_path),
      artifacts = list(import_only = TRUE, scale_provenance = scale_provenance),
      scale_provenance = scale_provenance,
      method = "cibersortx"
    ))
  }

  reference_scale <- .sn_validate_deconvolution_expression(
    reference_genes_by_cells,
    gene_axis = "rows",
    label = "CIBERSORTx reference",
    require_counts = FALSE
  )
  bulk_scale <- .sn_validate_deconvolution_expression(
    bulk_genes_by_samples,
    gene_axis = "rows",
    label = "CIBERSORTx mixture",
    require_counts = FALSE
  )
  cell_type_labels <- .sn_validate_deconvolution_labels(
    cell_type_labels, ncol(reference_genes_by_cells), "CIBERSORTx cell-type labels"
  )
  if (!identical(reference_scale, bulk_scale)) {
    stop(
      "CIBERSORTx reference and bulk inputs must use the same expression scale; detected `",
      reference_scale, "` and `", bulk_scale, "`.",
      call. = FALSE
    )
  }
  scale_provenance <- .sn_deconvolution_scale_provenance(reference_scale, bulk_scale)

  retain_run_dir <- !is.null(outdir)
  if (retain_run_dir && (!is.character(outdir) || length(outdir) != 1L ||
      is.na(outdir) || !nzchar(outdir))) {
    stop("`outdir` must be NULL or one non-empty parent directory.", call. = FALSE)
  }
  parent <- outdir %||% file.path(tempdir(), "shennong-cibersortx")
  run_dir <- .sn_create_owned_run_dir(parent = parent, prefix = "cibersortx_")
  work_dir <- file.path(run_dir, "input")
  if (!dir.create(work_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(work_dir)) {
    stop("Could not create the isolated CIBERSORTx work directory: ", work_dir, call. = FALSE)
  }

  complete <- FALSE
  on.exit({
    if (!complete && dir.exists(run_dir) && .sn_is_owned_run_dir(run_dir)) {
      try(.sn_sanitize_failed_python_run(run_dir, method = "cibersortx", stage = "execute"), silent = TRUE)
    }
  }, add = TRUE)
  result <- tryCatch(
    .sn_execute_cibersortx_run(
      reference_genes_by_cells = reference_genes_by_cells,
      cell_type_labels = cell_type_labels,
      bulk_genes_by_samples = bulk_genes_by_samples,
      work_dir = work_dir,
      prefix = prefix,
      email = email,
      token = token,
      container = container,
      container_path = container_path,
      dry_run = dry_run,
      verbose = verbose,
      rmbatch_b_mode = rmbatch_b_mode,
      rmbatch_s_mode = rmbatch_s_mode,
      perm = perm,
      qn = qn,
      absolute = absolute,
      abs_method = abs_method,
      k_max = k_max,
      max_dense_gb = max_dense_gb
    ),
    error = identity
  )
  if (inherits(result, "error")) {
    sanitization_complete <- tryCatch(
      .sn_sanitize_failed_python_run(run_dir, method = "cibersortx", stage = "execute"),
      error = function(...) FALSE
    )
    complete <- TRUE
    stop(
      conditionMessage(result),
      .sn_python_failure_suffix(run_dir, sanitization_complete),
      call. = FALSE
    )
  }

  result$scale_provenance <- scale_provenance
  result$artifacts$scale_provenance <- scale_provenance
  result$artifacts$run_dir_retained <- retain_run_dir
  result$artifacts$run_dir <- if (retain_run_dir) run_dir else NULL
  if (!retain_run_dir) {
    result$files <- list()
    .sn_cleanup_owned_run_dir(run_dir, label = "CIBERSORTx temporary run directory")
  }
  complete <- TRUE
  result
}

.sn_execute_cibersortx_run <- function(reference_genes_by_cells,
                                       cell_type_labels,
                                       bulk_genes_by_samples,
                                       work_dir,
                                       prefix,
                                       email = NULL,
                                       token = NULL,
                                       container = c("docker", "apptainer"),
                                       container_path = NULL,
                                       dry_run = FALSE,
                                       verbose = FALSE,
                                       rmbatch_b_mode = FALSE,
                                       rmbatch_s_mode = FALSE,
                                       perm = 0,
                                       qn = FALSE,
                                       absolute = FALSE,
                                       abs_method = "sig.score",
                                       k_max = 999,
                                       max_dense_gb = 2) {
  container <- match.arg(container)
  if (!isTRUE(dry_run)) {
    .sn_check_container_binary(container)
  }
  creds <- .sn_get_cibersortx_credentials(email = email, token = token)
  credentials_path <- .sn_write_cibersortx_credentials(creds$email, creds$token)
  on.exit(unlink(credentials_path, force = TRUE), add = TRUE)

  reference_path <- .sn_transform_and_save_cibersortx_single_cell(
    reference_genes_by_cells,
    cell_type_labels,
    work_dir,
    max_dense_gb = max_dense_gb
  )
  mixture_path <- .sn_transform_and_save_cibersortx_bulk(
    bulk_genes_by_samples,
    work_dir,
    max_dense_gb = max_dense_gb
  )

  signature_command <- .sn_create_cibersortx_command(
    input_dir = work_dir,
    output_dir = work_dir,
    credentials_path = credentials_path,
    container = container,
    container_path = container_path,
    method = "create_sig",
    verbose = verbose,
    refsample = basename(reference_path),
    k_max = k_max
  )
  signature_run <- .sn_run_command(signature_command, dry_run = dry_run, verbose = verbose)
  if (!identical(signature_run$status, 0L)) {
    stop("CIBERSORTx signature-matrix creation failed.", call. = FALSE)
  }

  signature_filename <- paste0(
    "CIBERSORTx_sample_file_for_cibersort_inferred_phenoclasses.CIBERSORTx_sample",
    "_file_for_cibersort_inferred_refsample.bm.K", as.integer(k_max), ".txt"
  )
  signature_path <- file.path(work_dir, signature_filename)
  label_by <- prefix

  fractions_command <- .sn_create_cibersortx_command(
    input_dir = work_dir,
    output_dir = work_dir,
    credentials_path = credentials_path,
    container = container,
    container_path = container_path,
    method = "impute_cell_fractions",
    verbose = verbose,
    refsample = basename(reference_path),
    sigmatrix = basename(signature_path),
    mixture = basename(mixture_path),
    label_by = label_by,
    rmbatch_b_mode = rmbatch_b_mode,
    rmbatch_s_mode = rmbatch_s_mode,
    perm = perm,
    qn = qn,
    absolute = absolute,
    abs_method = abs_method
  )
  fractions_run <- .sn_run_command(fractions_command, dry_run = dry_run, verbose = verbose)
  if (!identical(fractions_run$status, 0L)) {
    stop("CIBERSORTx fractions estimation failed.", call. = FALSE)
  }

  default_result_file <- if (isTRUE(rmbatch_b_mode) || isTRUE(rmbatch_s_mode)) {
    paste0("CIBERSORTx_", label_by, "_Adjusted.txt")
  } else {
    paste0("CIBERSORTx_", label_by, "_Results.txt")
  }
  result_path <- file.path(work_dir, default_result_file)
  if (!isTRUE(dry_run) && !file.exists(result_path)) {
    stop("CIBERSORTx completed without the expected fractions result: ", result_path, call. = FALSE)
  }
  imported <- if (!isTRUE(dry_run)) {
    .sn_import_cibersortx_fractions(result_path = result_path)
  } else {
    NULL
  }

  list(
    table = imported,
    files = list(
      single_cell_reference = reference_path,
      mixture = mixture_path,
      signature_matrix = signature_path,
      result = result_path
    ),
    artifacts = list(
      commands = list(
        create_signature = signature_run$command,
        deconvolve = fractions_run$command
      ),
      commands_redacted = TRUE,
      credential_transport = "ephemeral_read_only_bind_file",
      host_argv_contains_credentials = FALSE,
      dry_run = dry_run,
      container = container
    ),
    method = "cibersortx"
  )
}

.sn_import_cibersortx_fractions <- function(result_path) {
  result <- sn_read(result_path)
  result <- as.data.frame(result, stringsAsFactors = FALSE)

  if ("Mixture" %in% colnames(result)) {
    sample_names <- result$Mixture
    fraction_matrix <- as.matrix(result[, setdiff(colnames(result), c("Mixture", "P-value", "Correlation", "RMSE")) , drop = FALSE])
    rownames(fraction_matrix) <- sample_names
    return(.sn_fraction_matrix_to_long(fraction_matrix))
  }

  first_col <- colnames(result)[1]
  fraction_matrix <- as.matrix(result[, -1, drop = FALSE])
  rownames(fraction_matrix) <- result[[first_col]]
  .sn_fraction_matrix_to_long(t(fraction_matrix))
}

#' Store a deconvolution result on a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param result A deconvolution table.
#' @param result_id Stable identifier for the stored deconvolution result.
#' @param method Deconvolution backend, for example \code{"bayesprism"}.
#' @param bulk_samples Optional bulk sample identifiers.
#' @param reference_label Metadata column or label_by set used as the reference.
#' @param artifacts Optional backend-specific artifacts or file paths.
#' @param random_seed Optional random seed recorded in the result provenance.
#' @param return_object If \code{TRUE}, return the updated object.
#'
#' @return A \code{Seurat} object or stored-result list.
#' @export
sn_store_deconvolution <- function(object,
                                   result,
                                   result_id = "default",
                                   method = "bayesprism",
                                   bulk_samples = NULL,
                                   reference_label = NULL,
                                   artifacts = NULL,
                                   random_seed = NULL,
                                   return_object = TRUE) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)

  stored_result <- list(
    schema_version = .sn_analysis_result_schema_version(),
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = tibble::as_tibble(result),
    analysis = "deconvolution",
    method = method,
    bulk_samples = bulk_samples,
    reference_label = reference_label,
    artifacts = artifacts,
    provenance = .sn_analysis_provenance(random_seed = random_seed %||% NA_integer_)
  )

  object <- sn_store_result(
    object = object,
    type = "deconvolution",
    result_id = result_id,
    result = stored_result
  )

  if (isTRUE(return_object)) {
    return(.sn_log_seurat_command(object = object, name = "sn_store_deconvolution"))
  }

  sn_get_result(
    object = object,
    type = "deconvolution",
    result_id = result_id
  )
}

#' Retrieve a stored deconvolution result from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param result_id Name of the stored result.
#' @param samples Optional subset of bulk samples to keep.
#' @param cell_types Optional subset of cell types to keep.
#' @param with_metadata If \code{TRUE}, return the full stored-result list.
#'
#' @return A tibble or stored-result list.
#' @export
sn_get_deconvolution_result <- function(object,
                                        result_id = "default",
                                        samples = NULL,
                                        cell_types = NULL,
                                        with_metadata = FALSE) {
  .sn_validate_seurat_object(object)

  stored <- sn_get_result(
    object = object,
    type = "deconvolution",
    result_id = result_id
  )
  if (isTRUE(with_metadata)) {
    return(stored)
  }

  table <- tibble::as_tibble(stored$tables$primary)
  if (!is.null(samples) && "sample" %in% colnames(table)) {
    table <- dplyr::filter(table, .data$sample %in% samples)
  }
  if (!is.null(cell_types) && "cell_type" %in% colnames(table)) {
    table <- dplyr::filter(table, .data$cell_type %in% cell_types)
  }
  table
}

#' Deprecated alias of `sn_run_bulk_deconvolution()`
#'
#' `sn_deconvolve_bulk()` is a deprecated compatibility alias. Use [sn_run_bulk_deconvolution()] directly;
#' the alias will be removed in a future release.
#'
#' @param ... Named arguments passed on to [sn_run_bulk_deconvolution()].
#'
#' @return Result of `sn_run_bulk_deconvolution(...)`.
#'
#' @export
sn_deconvolve_bulk <- function(...) {
  .Deprecated("sn_run_bulk_deconvolution", package = "Shennong")
  sn_run_bulk_deconvolution(...)
}
