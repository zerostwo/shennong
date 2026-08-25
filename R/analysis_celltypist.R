.sn_write_celltypist_sparse_input <- function(counts,
                                               outdir,
                                               transpose_input = TRUE,
                                               gene_file = NULL,
                                               cell_file = NULL) {
  counts <- .sn_as_sparse_matrix(counts)
  if (!inherits(counts, "Matrix")) {
    stop("CellTypist export requires a matrix-like counts layer.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("CellTypist export requires feature and cell names.", call. = FALSE)
  }

  # CellTypist consumes a cell-by-gene Matrix Market file unless its
  # --transpose-input flag is supplied. Preserve the established flag while
  # keeping the in-memory representation sparse in either orientation.
  export_counts <- if (isTRUE(transpose_input)) counts else Matrix::t(counts)
  input_data <- file.path(outdir, "counts.mtx")
  Matrix::writeMM(obj = export_counts, file = input_data)

  if (is.null(gene_file)) {
    gene_file <- file.path(outdir, "genes.csv")
    utils::write.table(
      data.frame(feature = rownames(counts), check.names = FALSE),
      file = gene_file,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      quote = TRUE
    )
  }
  if (is.null(cell_file)) {
    cell_file <- file.path(outdir, "cells.csv")
    utils::write.table(
      data.frame(cell = colnames(counts), check.names = FALSE),
      file = cell_file,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      quote = TRUE
    )
  }

  list(
    input_data = input_data,
    gene_file = gene_file,
    cell_file = cell_file
  )
}


#' Run CellTypist for automated cell type annotation
#'
#' @param x A Seurat object or a path to a count matrix / AnnData file that CellTypist can consume.
#' @param celltypist Path to the `celltypist` binary. When `NULL`, use
#'   `getOption("shennong.celltypist_path")` and otherwise search `PATH` with
#'   `Sys.which("celltypist")`.
#' @param model Model used for predictions. Defaults to "Immune_All_Low.pkl".
#' @param outdir Directory to store the output files. If NULL, use a temporary directory.
#' @param prefix Prefix for the output files. By default, use the model name plus a dot.
#' @param mode Choose the cell type with the largest score/probability (`"best_match"`) or enable multi-label classification (`"prob_match"`).
#' @param p_thres Probability threshold for the multi-label classification. Ignored if `mode = "best_match"`.
#' @param majority_voting Logical. Whether to refine labels using majority voting after over-clustering.
#' @param over_clustering Input file or a string key specifying an existing metadata column in the AnnData object, or "auto".
#' @param min_prop For the dominant cell type within a subcluster, the minimum proportion of cells required to name the subcluster by this cell type.
#' @param transpose_input Logical. For Seurat input, `TRUE` exports counts in
#'   gene-by-cell orientation and `FALSE` exports a sparse cell-by-gene
#'   transpose. For an existing path, the input file is not rewritten. In both
#'   cases, `TRUE` adds CellTypist's `--transpose-input` flag and `FALSE` does
#'   not, so path inputs must set this to match their stored orientation.
#' @param gene_file If the provided input is in the `mtx` format, path to the
#'   file storing gene information. For Seurat input, a sidecar is generated
#'   from the feature names when this is `NULL`.
#' @param cell_file If the provided input is in the `mtx` format, path to the
#'   file storing cell information. For Seurat input, a sidecar is generated
#'   from the cell names when this is `NULL`.
#' @param assay Assay used when exporting Seurat counts to CellTypist. Defaults to \code{"RNA"}.
#' @param layer Raw or count-like layer used as the input matrix for Seurat
#'   objects. Defaults to \code{"counts"}. MatrixMarket/CSV input is normalized
#'   by CellTypist, so passing a log-normalized `data` layer would normalize it
#'   twice and is rejected; negative values are also rejected.
#' @param xlsx Logical. If `TRUE`, ask CellTypist for its combined
#'   `annotation_result.xlsx` workbook and import the first (prediction) sheet
#'   through the optional `rio` dependency. Defaults to `FALSE`.
#' @param plot_results Logical. If `TRUE`, plot the prediction results. Defaults to `FALSE`.
#' @param quiet Logical. If `TRUE`, hide the banner and config info from `celltypist`. Defaults to `FALSE`.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#'
#' @return When \code{x} is a Seurat object, a Seurat object with prediction
#'   columns added to metadata. When \code{x} is a path, the CellTypist
#'   prediction table is returned.
#'
#' @examples
#' \dontrun{
#' pbmc <- qs2::qs_read(file.path(Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"))
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
                              quiet = FALSE,
                              object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  check_installed(c("logger", "glue"),
    reason = "to run CellTypist analysis."
  )
  x_is_seurat <- inherits(x, "Seurat")
  if (x_is_seurat) {
    check_installed("Seurat", reason = "to annotate a Seurat object with CellTypist.")
  }

  mode <- match.arg(mode)
  celltypist <- celltypist %||% getOption("shennong.celltypist_path", Sys.which("celltypist"))
  sn_check_file(celltypist)
  if (isTRUE(xlsx)) {
    check_installed("rio", reason = "to import CellTypist's Excel result workbook.")
  }

  .sn_log_info("Starting CellTypist analysis with model = {model}.")
  tictoc::tic("Total CellTypist runtime")

  temporary_outdir <- is_null(outdir)
  if (temporary_outdir) {
    outdir <- tempfile("celltypist_")
    dir.create(outdir, recursive = TRUE)
    .sn_log_info("Using a temporary output directory: {outdir}.")
  } else {
    outdir <- sn_set_path(outdir)
    .sn_log_info("Using the user-specified output directory: {outdir}.")
  }

  on.exit(
    {
      if (temporary_outdir && dir.exists(outdir)) {
        unlink(outdir, recursive = TRUE)
        log_debug("Cleaned temporary directory: {outdir}")
      }
    },
    add = TRUE
  )

  if (x_is_seurat) {
    .sn_log_info("Converting the Seurat object to sparse CellTypist input format.")

    counts <- tryCatch(
      .sn_get_seurat_layer_data(object = x, assay = assay, layer = layer),
      error = function(e) {
        .sn_log_error("Failed to extract the counts layer: {e$message}.")
        stop("Counts layer extraction failed")
      }
    )
    counts <- .sn_as_sparse_matrix(counts)
    nonzero_values <- if ("x" %in% methods::slotNames(counts)) {
      methods::slot(counts, "x")
    } else {
      numeric()
    }
    if (any(!is.finite(nonzero_values)) || any(nonzero_values < 0)) {
      stop(
        "CellTypist MatrixMarket input must contain finite, non-negative count-like values.",
        call. = FALSE
      )
    }
    if (tolower(layer) %in% c("data", "scale.data", "scaled.data")) {
      stop(
        "CellTypist treats MatrixMarket input as raw counts and normalizes it internally; ",
        "use a raw/count-like layer instead of `", layer, "`.",
        call. = FALSE
      )
    }

    exported_input <- .sn_write_celltypist_sparse_input(
      counts = counts,
      outdir = outdir,
      transpose_input = transpose_input,
      gene_file = gene_file,
      cell_file = cell_file
    )
    input_data <- exported_input$input_data
    gene_file <- exported_input$gene_file
    cell_file <- exported_input$cell_file
    log_debug("Count matrix written to {input_data} ({file.size(input_data)} bytes)")
  } else {
    input_data <- x
    .sn_log_info("Using the precomputed input matrix: {input_data}.")
  }

  over_clustering_path <- NULL

  if (!is_null(over_clustering) && !identical(over_clustering, "auto")) {
    if (isTRUE(x_is_seurat)) {
      if (over_clustering %in% colnames(x@meta.data)) {
        .sn_log_info("Using the existing clustering column: {over_clustering}.")
        over_clustering_path <- file.path(outdir, "over_clustering.txt")
        writeLines(
          as.character(x[[over_clustering, drop = TRUE]]),
          over_clustering_path
        )
      } else if (file.exists(over_clustering)) {
        over_clustering_path <- over_clustering
      } else {
        stop(
          "`over_clustering` must be a metadata column or existing file when `x` is a Seurat object.",
          call. = FALSE
        )
      }
    } else {
      over_clustering_path <- over_clustering
    }
  } else if (identical(over_clustering, "auto") && isTRUE(majority_voting)) {
    .sn_log_info("Using CellTypist's automatic over-clustering for majority voting.")
  }

  model_name <- tools::file_path_sans_ext(basename(model))
  prefix <- prefix %||% glue("{model_name}.")
  prediction_path <- if (isTRUE(xlsx)) {
    file.path(outdir, glue("{prefix}annotation_result.xlsx"))
  } else {
    file.path(outdir, glue("{prefix}predicted_labels.csv"))
  }

  cmd_args <- c(
    "--indata", shQuote(input_data),
    "--model", shQuote(model),
    "--mode", mode,
    "--outdir", shQuote(outdir),
    "--prefix", shQuote(prefix)
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

  if (!file.exists(prediction_path)) {
    .sn_log_error("Prediction output is missing: {prediction_path}.")
    stop("CellTypist did not generate expected output files")
  }

  predicted_labels <- if (isTRUE(xlsx)) {
    sn_read(
      prediction_path,
      format = "xlsx",
      which = 1,
      row_names = 1
    )
  } else {
    utils::read.csv(prediction_path, row.names = 1, check.names = FALSE)
  }
  predicted_labels <- as.data.frame(predicted_labels, check.names = FALSE)
  normalized_columns <- tolower(gsub("[^A-Za-z0-9]+", "_", colnames(predicted_labels)))
  normalized_columns <- gsub("^_+|_+$", "", normalized_columns)
  if (!"predicted_labels" %in% normalized_columns) {
    stop(
      "CellTypist prediction output does not contain a `predicted_labels` column.",
      call. = FALSE
    )
  }

  probability_matrix <- if (isTRUE(xlsx)) {
    tryCatch(
      sn_read(
        prediction_path,
        format = "xlsx",
        which = 3,
        row_names = 1
      ),
      error = function(e) NULL
    )
  } else {
    probability_path <- file.path(outdir, glue("{prefix}probability_matrix.csv"))
    if (file.exists(probability_path)) {
      utils::read.csv(probability_path, row.names = 1, check.names = FALSE)
    } else {
      NULL
    }
  }
  if (!is_null(probability_matrix)) {
    probability_matrix <- as.data.frame(probability_matrix, check.names = FALSE)
    label_column <- match("majority_voting", normalized_columns)
    if (is.na(label_column)) label_column <- match("predicted_labels", normalized_columns)
    selected_labels <- as.character(predicted_labels[[label_column]])
    probability_columns <- tolower(gsub(
      "[^A-Za-z0-9]+", "_", colnames(probability_matrix)
    ))
    probability_columns <- gsub("^_+|_+$", "", probability_columns)
    probability_rows <- match(rownames(predicted_labels), rownames(probability_matrix))
    confidence <- vapply(seq_along(selected_labels), function(index) {
      row_index <- probability_rows[[index]]
      if (is.na(row_index)) return(NA_real_)
      label_parts <- trimws(unlist(strsplit(selected_labels[[index]], "[|,]")))
      label_parts <- tolower(gsub("[^A-Za-z0-9]+", "_", label_parts))
      label_parts <- gsub("^_+|_+$", "", label_parts)
      column_indices <- match(label_parts, probability_columns)
      column_indices <- column_indices[!is.na(column_indices)]
      if (!length(column_indices)) return(NA_real_)
      values <- suppressWarnings(as.numeric(unlist(
        probability_matrix[row_index, column_indices, drop = FALSE],
        use.names = FALSE
      )))
      values <- values[is.finite(values)]
      if (length(values)) max(values) else NA_real_
    }, numeric(1))
    predicted_labels$confidence <- confidence
    normalized_columns <- c(normalized_columns, "confidence")
  }

  model_key <- gsub("\\.pkl$", "", model_name)
  colnames(predicted_labels) <- make.unique(
    paste0(model_key, "_", normalized_columns),
    sep = "_"
  )

  if (!isTRUE(x_is_seurat)) {
    tictoc::toc()
    .sn_log_info("CellTypist analysis completed successfully.")
    return(tibble::as_tibble(predicted_labels, rownames = "cell"))
  }

  .sn_log_info("Adding {ncol(predicted_labels)} metadata columns to the Seurat object.")
  x <- SeuratObject::AddMetaData(x, metadata = predicted_labels)

  tictoc::toc()
  .sn_log_info("CellTypist analysis completed successfully.")

  .sn_log_seurat_command(object = x, assay = assay, name = "sn_run_celltypist")
}
