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
#' @param cell_type_col Metadata column containing cell-type labels when
#'   \code{x} is a \code{Seurat} object. If \code{x} is a matrix, supply
#'   \code{cell_type_labels} instead.
#' @param cell_state_col Optional metadata column containing cell-state labels.
#'   Defaults to \code{cell_type_col}.
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
#' @param key Optional malignant-cell label passed to BayesPrism.
#' @param outdir Output directory used by CIBERSORTx file export.
#' @param prefix Prefix used for exported files and stored result names.
#' @param store_name Name used under \code{object@misc$deconvolution_results}
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
#' @param cibersortx_rmbatch_b_mode Whether to enable B-mode batch correction.
#' @param cibersortx_rmbatch_s_mode Whether to enable S-mode batch correction.
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
#' @param return_object If \code{TRUE} and \code{x} is a \code{Seurat} object,
#'   return the updated object when a result table is available.
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
#'   bundle <- sn_deconvolve_bulk(
#'     ref,
#'     bulk = bulk,
#'     method = "cibersortx",
#'     cell_type_col = "cell_type",
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
sn_deconvolve_bulk <- function(x,
                               bulk,
                               method = c("bayesprism", "cibersortx"),
                               cell_type_col = NULL,
                               cell_state_col = NULL,
                               cell_type_labels = NULL,
                               cell_state_labels = NULL,
                               assay = "RNA",
                               layer = "counts",
                               bulk_gene_axis = c("auto", "rows", "columns"),
                               key = NULL,
                               outdir = tempdir(),
                               prefix = "deconvolution",
                               store_name = "default",
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
                               return_object = TRUE) {
  method <- match.arg(method)
  bulk_gene_axis <- match.arg(bulk_gene_axis)
  cibersortx_container <- match.arg(cibersortx_container)

  reference_info <- .sn_prepare_deconvolution_reference(
    x = x,
    cell_type_col = cell_type_col,
    cell_state_col = cell_state_col,
    cell_type_labels = cell_type_labels,
    cell_state_labels = cell_state_labels,
    assay = assay,
    layer = layer
  )
  bulk_samples_by_gene <- .sn_prepare_bulk_matrix(
    bulk = bulk,
    reference_genes = colnames(reference_info$reference_cells_by_gene),
    gene_axis = bulk_gene_axis
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
      update_gibbs = update_gibbs
    )
  } else {
    result <- .sn_run_cibersortx_local(
      reference_genes_by_cells = t(reference_info$reference_cells_by_gene),
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
      k_max = cibersortx_k_max
    )
  }

  if (!inherits(x, "Seurat") || is.null(result$table)) {
    return(result)
  }

  object <- sn_store_deconvolution(
    object = x,
    result = result$table,
    store_name = store_name,
    method = method,
    bulk_samples = rownames(bulk_samples_by_gene),
    reference_label = reference_info$reference_label,
    artifacts = result$artifacts %||% result$files
  )

  if (isTRUE(return_object)) {
    return(.sn_log_seurat_command(object = object, name = "sn_deconvolve_bulk"))
  }

  .sn_get_misc_result(
    object = object,
    collection = "deconvolution_results",
    store_name = store_name
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
                                                layer = "counts") {
  if (inherits(x, "Seurat")) {
    if (is.null(cell_type_col) || !cell_type_col %in% colnames(x[[]])) {
      stop("`cell_type_col` must be supplied and present in `x@meta.data`.", call. = FALSE)
    }
    if (!is.null(cell_state_col) && !cell_state_col %in% colnames(x[[]])) {
      stop(glue("Metadata column '{cell_state_col}' was not found."), call. = FALSE)
    }

    counts <- .sn_get_seurat_layer_data(object = x, assay = assay, layer = layer)
    reference <- t(as.matrix(counts))
    cell_type_labels <- as.character(x[[cell_type_col, drop = TRUE]])
    cell_state_labels <- if (is.null(cell_state_col)) {
      cell_type_labels
    } else {
      as.character(x[[cell_state_col, drop = TRUE]])
    }
    reference_label <- cell_type_col
  } else {
    reference <- .sn_as_matrix(x, name = "x")
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

  keep <- !is.na(cell_type_labels) & !is.na(cell_state_labels)
  reference <- reference[keep, , drop = FALSE]
  cell_type_labels <- cell_type_labels[keep]
  cell_state_labels <- cell_state_labels[keep]

  list(
    reference_cells_by_gene = reference,
    cell_type_labels = cell_type_labels,
    cell_state_labels = cell_state_labels,
    reference_label = reference_label
  )
}

.sn_as_matrix <- function(x, name = "input") {
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
                                    gene_axis = c("auto", "rows", "columns")) {
  gene_axis <- match.arg(gene_axis)
  bulk <- .sn_as_matrix(bulk, name = "bulk")
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

.sn_fraction_matrix_to_long <- function(fraction_matrix,
                                        sample_col = "sample",
                                        cell_type_col = "cell_type",
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
  names(table) <- c(sample_col, cell_type_col, value_col)
  tibble::as_tibble(table)
}

.sn_run_bayesprism <- function(reference_cells_by_gene,
                               cell_type_labels,
                               cell_state_labels,
                               bulk_samples_by_gene,
                               key = NULL,
                               gibbs_control = list(),
                               opt_control = list(),
                               n_cores = 1,
                               update_gibbs = TRUE) {
  check_installed_github(
    pkg = "BayesPrism",
    repo = "Danko-Lab/BayesPrism/BayesPrism",
    reason = "to run BayesPrism bulk deconvolution."
  )

  prism <- BayesPrism::new.prism(
    reference = reference_cells_by_gene,
    input.type = "count.matrix",
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_state_labels,
    key = key,
    mixture = bulk_samples_by_gene
  )
  bp_fit <- BayesPrism::run.prism(
    prism = prism,
    n.cores = n_cores,
    update.gibbs = update_gibbs,
    gibbs.control = gibbs_control,
    opt.control = opt_control
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
      prism = prism
    ),
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
  list(email = email, token = token)
}

.sn_check_container_binary <- function(container = c("docker", "apptainer")) {
  container <- match.arg(container)
  binary <- Sys.which(container)
  if (!nzchar(binary)) {
    stop(glue("Could not find the `{container}` executable required for local CIBERSORTx runs."), call. = FALSE)
  }
  binary
}

.sn_run_command <- function(command, dry_run = FALSE, verbose = FALSE) {
  if (isTRUE(dry_run)) {
    return(list(status = 0L, command = command))
  }

  status <- system(command, ignore.stdout = !verbose, ignore.stderr = !verbose)
  list(status = status, command = command)
}

.sn_transform_and_save_cibersortx_single_cell <- function(reference_genes_by_cells, cell_type_labels, path) {
  output_file <- file.path(path, "sample_file_for_cibersort.txt")
  reference_genes_by_cells <- as.matrix(reference_genes_by_cells)
  colnames(reference_genes_by_cells) <- as.character(cell_type_labels)
  output <- rbind(colnames(reference_genes_by_cells), reference_genes_by_cells)
  rownames(output) <- c("GeneSymbol", rownames(reference_genes_by_cells))
  output <- data.frame("GeneSymbol" = rownames(output), output, stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.table(output, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  output_file
}

.sn_transform_and_save_cibersortx_bulk <- function(bulk_genes_by_samples, path) {
  output_file <- file.path(path, "mixture_file_for_cibersort.txt")
  output <- data.frame("Gene" = rownames(bulk_genes_by_samples), bulk_genes_by_samples, stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.table(output, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  output_file
}

.sn_create_cibersortx_command <- function(input_dir,
                                          output_dir,
                                          email,
                                          token,
                                          container = c("docker", "apptainer"),
                                          container_path = NULL,
                                          method = c("create_sig", "impute_cell_fractions"),
                                          verbose = FALSE,
                                          refsample = "sample_file_for_cibersort.txt",
                                          sigmatrix = "signature_matrix.txt",
                                          mixture = "mixture_file_for_cibersort.txt",
                                          label = "shennong",
                                          rmbatch_b_mode = FALSE,
                                          rmbatch_s_mode = FALSE,
                                          perm = 0,
                                          qn = FALSE,
                                          absolute = FALSE,
                                          abs_method = "sig.score",
                                          k_max = 999) {
  container <- match.arg(container)
  method <- match.arg(method)

  if (identical(container, "docker")) {
    base <- paste0(
      "docker run -v ", shQuote(normalizePath(input_dir)), ":/src/data:z -v ",
      shQuote(normalizePath(output_dir)), ":/src/outdir:z cibersortx/fractions --single_cell TRUE"
    )
  } else {
    if (is.null(container_path) || !file.exists(path.expand(container_path))) {
      stop("`cibersortx_container_path` must point to an existing Apptainer image.", call. = FALSE)
    }
    base <- paste0(
      "apptainer exec --no-home -c -B ", shQuote(paste0(normalizePath(input_dir), "/:/src/data")),
      " -B ", shQuote(paste0(normalizePath(output_dir), "/:/src/outdir")),
      " ", shQuote(normalizePath(path.expand(container_path))), " /src/CIBERSORTxFractions --single_cell TRUE"
    )
  }

  if (isTRUE(verbose)) {
    base <- paste(base, "--verbose TRUE")
  }
  credentials <- paste("--username", shQuote(email), "--token", shQuote(token))

  options <- if (identical(method, "create_sig")) {
    paste(
      "--refsample", shQuote(refsample),
      "--G.min 300 --G.max 500 --q.value 0.01 --filter FALSE",
      "--k.max", as.integer(k_max),
      "--remake FALSE --replicates 5 --sampling 0.5 --fraction 0.75"
    )
  } else {
    option_string <- paste(
      "--mixture", shQuote(mixture),
      "--sigmatrix", shQuote(sigmatrix),
      "--perm", as.integer(perm),
      "--label", shQuote(label),
      "--rmbatchBmode", toupper(as.character(isTRUE(rmbatch_b_mode))),
      "--rmbatchSmode", toupper(as.character(isTRUE(rmbatch_s_mode))),
      "--sourceGEPs", shQuote(sigmatrix),
      "--QN", toupper(as.character(isTRUE(qn))),
      "--absolute", toupper(as.character(isTRUE(absolute))),
      "--abs_method", shQuote(abs_method)
    )
    if (isTRUE(rmbatch_b_mode) || isTRUE(rmbatch_s_mode)) {
      option_string <- paste(option_string, "--refsample", shQuote(refsample))
    }
    option_string
  }

  paste(base, credentials, options)
}

.sn_run_cibersortx_local <- function(reference_genes_by_cells,
                                     cell_type_labels,
                                     bulk_genes_by_samples,
                                     outdir,
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
                                     k_max = 999) {
  if (!is.null(result_path) && file.exists(result_path) && !isTRUE(dry_run)) {
    return(list(
      table = .sn_import_cibersortx_fractions(result_path = result_path),
      files = list(result = result_path),
      artifacts = list(import_only = TRUE),
      method = "cibersortx"
    ))
  }

  outdir <- sn_set_path(outdir)
  container <- match.arg(container)
  if (!isTRUE(dry_run)) {
    .sn_check_container_binary(container)
  }
  creds <- .sn_get_cibersortx_credentials(email = email, token = token)

  reference_path <- .sn_transform_and_save_cibersortx_single_cell(reference_genes_by_cells, cell_type_labels, outdir)
  mixture_path <- .sn_transform_and_save_cibersortx_bulk(bulk_genes_by_samples, outdir)

  signature_command <- .sn_create_cibersortx_command(
    input_dir = outdir,
    output_dir = outdir,
    email = creds$email,
    token = creds$token,
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
  signature_path <- file.path(outdir, signature_filename)
  label <- prefix

  fractions_command <- .sn_create_cibersortx_command(
    input_dir = outdir,
    output_dir = outdir,
    email = creds$email,
    token = creds$token,
    container = container,
    container_path = container_path,
    method = "impute_cell_fractions",
    verbose = verbose,
    refsample = basename(reference_path),
    sigmatrix = basename(signature_path),
    mixture = basename(mixture_path),
    label = label,
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

  default_result_name <- if (isTRUE(rmbatch_b_mode) || isTRUE(rmbatch_s_mode)) {
    paste0("CIBERSORTx_", label, "_Adjusted.txt")
  } else {
    paste0("CIBERSORTx_", label, "_Results.txt")
  }
  result_path <- result_path %||% file.path(outdir, default_result_name)
  imported <- if (!isTRUE(dry_run) && file.exists(result_path)) .sn_import_cibersortx_fractions(result_path = result_path) else NULL

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
        create_signature = signature_command,
        deconvolve = fractions_command
      ),
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
#' @param store_name Name used under \code{object@misc$deconvolution_results}.
#' @param method Deconvolution backend, for example \code{"bayesprism"}.
#' @param bulk_samples Optional bulk sample identifiers.
#' @param reference_label Metadata column or label set used as the reference.
#' @param artifacts Optional backend-specific artifacts or file paths.
#' @param return_object If \code{TRUE}, return the updated object.
#'
#' @return A \code{Seurat} object or stored-result list.
#' @export
sn_store_deconvolution <- function(object,
                                   result,
                                   store_name = "default",
                                   method = "bayesprism",
                                   bulk_samples = NULL,
                                   reference_label = NULL,
                                   artifacts = NULL,
                                   return_object = TRUE) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  stored_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = tibble::as_tibble(result),
    analysis = "deconvolution",
    method = method,
    bulk_samples = bulk_samples,
    reference_label = reference_label,
    artifacts = artifacts
  )

  object <- .sn_store_misc_result(
    object = object,
    collection = "deconvolution_results",
    store_name = store_name,
    result = stored_result
  )

  if (isTRUE(return_object)) {
    return(.sn_log_seurat_command(object = object, name = "sn_store_deconvolution"))
  }

  stored_result
}

#' Retrieve a stored deconvolution result from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param deconvolution_name Name of the stored result.
#' @param samples Optional subset of bulk samples to keep.
#' @param cell_types Optional subset of cell types to keep.
#' @param with_metadata If \code{TRUE}, return the full stored-result list.
#'
#' @return A tibble or stored-result list.
#' @export
sn_get_deconvolution_result <- function(object,
                                        deconvolution_name = "default",
                                        samples = NULL,
                                        cell_types = NULL,
                                        with_metadata = FALSE) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.")
  }

  stored <- .sn_get_misc_result(
    object = object,
    collection = "deconvolution_results",
    store_name = deconvolution_name
  )
  if (isTRUE(with_metadata)) {
    return(stored)
  }

  table <- tibble::as_tibble(stored$table)
  if (!is.null(samples) && "sample" %in% colnames(table)) {
    table <- dplyr::filter(table, .data$sample %in% samples)
  }
  if (!is.null(cell_types) && "cell_type" %in% colnames(table)) {
    table <- dplyr::filter(table, .data$cell_type %in% cell_types)
  }
  table
}
