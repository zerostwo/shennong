.make_conformance_seurat_object <- function(seed = 2402L, features = 40L, cells = 24L) {
  .conformance_require_package("Seurat")
  .conformance_require_package("SeuratObject")
  set.seed(seed)
  counts <- matrix(
    stats::rpois(features * cells, lambda = 3),
    nrow = features,
    dimnames = list(
      paste0("gene", seq_len(features)),
      paste0("cell_", seq_len(cells))
    )
  )
  SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
}

.seurat_normalize_projection <- function(object) {
  list(
    counts = SeuratObject::LayerData(object, assay = "RNA", layer = "counts"),
    data = SeuratObject::LayerData(object, assay = "RNA", layer = "data"),
    default_assay = SeuratObject::DefaultAssay(object),
    features = rownames(object),
    cells = colnames(object),
    metadata = object[[]],
    variable_features = SeuratObject::VariableFeatures(object[["RNA"]])
  )
}

test_that("sn_normalize_data matches Seurat::NormalizeData in its pilot envelope", {
  contract <- .conformance_contract("preprocessing::seurat_normalize")
  object <- .make_conformance_seurat_object()
  object$batch <- rep(c("a", "b"), each = ncol(object) / 2L)
  before <- .conformance_fingerprint(object)
  cases <- list(
    list(
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      margin = 1,
      verbose = FALSE
    ),
    list(
      normalization.method = "LogNormalize",
      scale.factor = 2345,
      margin = 1,
      verbose = FALSE
    )
  )

  for (arguments in cases) {
    upstream <- .conformance_without_acceleration(do.call(
      Seurat::NormalizeData,
      c(list(object = object), arguments)
    ))
    candidate <- .conformance_without_acceleration(do.call(
      sn_normalize_data,
      c(list(object = object, method = "seurat"), arguments)
    ))

    .conformance_expect_equal(
      .seurat_normalize_projection(candidate),
      .seurat_normalize_projection(upstream),
      contract,
      info = paste(
        "Seurat", as.character(utils::packageVersion("Seurat")),
        "scale.factor", arguments$scale.factor
      )
    )
    .conformance_expect_unchanged(object, before, "Seurat normalization input")
  }
})

.make_conformance_silhouette_object <- function() {
  object <- .make_conformance_seurat_object(seed = 2401L, features = 20L, cells = 24L)
  set.seed(2401L)
  embeddings <- matrix(
    stats::rnorm(ncol(object) * 5L),
    nrow = ncol(object),
    dimnames = list(colnames(object), paste0("PC_", seq_len(5L)))
  )
  object[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = embeddings,
    key = "PC_",
    assay = "RNA"
  )
  object$group <- rep(c("T", "B", "Mono"), each = 8L)
  object
}

.silhouette_oracle <- function(object, candidate, dims = NULL) {
  cells <- candidate$cell_id
  embeddings <- SeuratObject::Embeddings(object, reduction = "pca")
  if (!is.null(dims)) {
    embeddings <- embeddings[, dims, drop = FALSE]
  }
  embeddings <- embeddings[cells, , drop = FALSE]
  labels <- factor(object[[]][cells, "group"])
  result <- cluster::silhouette(
    x = as.integer(labels),
    dist = stats::dist(embeddings)
  )
  data.frame(
    cell_id = cells,
    group = as.character(labels),
    silhouette_width = unname(result[, "sil_width"]),
    stringsAsFactors = FALSE
  )
}

test_that("sn_calculate_silhouette matches cluster::silhouette across parameter cases", {
  contract <- .conformance_contract("metrics::silhouette")
  object <- .make_conformance_silhouette_object()
  before <- .conformance_fingerprint(object)
  cases <- list(
    list(dims = NULL, cells = NULL, max_cells = NULL, stratify_by = "group", seed = 717L),
    list(
      dims = c(1L, 3L, 5L),
      cells = paste0("cell_", c(1:7, 10:24)),
      max_cells = NULL,
      stratify_by = "group",
      seed = 82L
    ),
    list(dims = 1:4, cells = NULL, max_cells = 15L, stratify_by = "group", seed = 91L)
  )

  for (arguments in cases) {
    candidate <- do.call(
      sn_calculate_silhouette,
      c(list(x = object, label_by = "group", reduction = "pca"), arguments)
    )
    upstream <- .silhouette_oracle(object, candidate, dims = arguments$dims)
    .conformance_expect_equal(
      as.data.frame(candidate),
      upstream,
      contract,
      info = paste("cluster", as.character(utils::packageVersion("cluster")))
    )
    repeat_candidate <- do.call(
      sn_calculate_silhouette,
      c(list(x = object, label_by = "group", reduction = "pca"), arguments)
    )
    expect_identical(candidate, repeat_candidate)
    .conformance_expect_unchanged(object, before, "silhouette input")
  }
})

.make_conformance_edger_fixture <- function() {
  set.seed(2403L)
  samples <- paste0("sample_", seq_len(12L))
  metadata <- data.frame(
    condition = factor(
      rep(c("control", "treated"), each = 6L),
      levels = c("control", "treated")
    ),
    batch = factor(rep(seq_len(3L), times = 4L)),
    row.names = samples
  )
  counts <- matrix(
    stats::rnbinom(120L * 12L, mu = 40, size = 5),
    nrow = 120L,
    dimnames = list(paste0("gene_", seq_len(120L)), samples)
  )
  counts[seq_len(12L), metadata$condition == "treated"] <-
    counts[seq_len(12L), metadata$condition == "treated"] + 45L
  list(counts = counts, metadata = metadata)
}

.edger_contrast_vector <- function(design_matrix, contrast) {
  coefficient <- paste0(contrast[[1]], contrast[[2]])
  index <- match(coefficient, colnames(design_matrix))
  if (is.na(index)) {
    stop("The independent edgeR oracle could not identify coefficient ", coefficient, ".")
  }
  vector <- numeric(ncol(design_matrix))
  vector[[index]] <- 1
  vector
}

.edger_oracle <- function(counts, metadata, design, contrast, robust) {
  design_matrix <- stats::model.matrix(design, data = metadata)
  y <- edgeR::DGEList(counts = round(counts))
  keep <- edgeR::filterByExpr(y, design = design_matrix)
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (exists("normLibSizes", envir = asNamespace("edgeR"), inherits = FALSE)) {
    y <- edgeR::normLibSizes(y)
  } else {
    y <- edgeR::calcNormFactors(y)
  }
  y <- edgeR::estimateDisp(y, design_matrix, robust = robust)
  fit <- edgeR::glmQLFit(y, design_matrix, robust = robust)
  test <- edgeR::glmQLFTest(
    fit,
    contrast = .edger_contrast_vector(design_matrix, contrast)
  )
  table <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  output <- data.frame(
    gene = rownames(table),
    log2_fold_change = as.numeric(table$logFC),
    statistic = as.numeric(table$F),
    p_value = as.numeric(table$PValue),
    adjusted_p_value = as.numeric(table$FDR),
    base_mean = as.numeric(table$logCPM),
    stringsAsFactors = FALSE
  )
  output <- output[order(output$adjusted_p_value, output$p_value), , drop = FALSE]
  rownames(output) <- NULL
  list(
    table = output,
    keep = keep,
    contrast = .edger_contrast_vector(design_matrix, contrast)
  )
}

test_that("sn_find_bulk_de matches a direct edgeR QL pipeline", {
  .conformance_require_package("edgeR")
  contract <- .conformance_contract("bulk_de::edger")
  fixture <- .make_conformance_edger_fixture()
  before <- .conformance_fingerprint(fixture)
  contrast <- c("condition", "treated", "control")
  cases <- list(
    list(design = ~condition, robust = TRUE, method = "edger"),
    list(design = ~batch + condition, robust = FALSE, method = "edger"),
    list(design = ~batch + condition, robust = TRUE, method = "auto")
  )

  for (case in cases) {
    upstream <- .edger_oracle(
      fixture$counts,
      fixture$metadata,
      case$design,
      contrast,
      robust = case$robust
    )
    candidate <- sn_find_bulk_de(
      fixture$counts,
      metadata = fixture$metadata,
      design = case$design,
      contrast = contrast,
      method = case$method,
      backend_control = list(robust = case$robust)
    )
    projection <- unlist(contract$equivalence$projection, use.names = FALSE)
    actual <- as.data.frame(candidate$tables$primary[, projection, drop = FALSE])
    rownames(actual) <- NULL

    .conformance_expect_equal(
      actual,
      upstream$table,
      contract,
      info = paste(
        "edgeR", as.character(utils::packageVersion("edgeR")),
        "method", case$method,
        "robust", case$robust
      )
    )
    expect_identical(candidate$method, "edger")
    expect_equal(candidate$diagnostics$retained_features, sum(upstream$keep))
    expect_equal(
      as.numeric(candidate$diagnostics$contrast_vector),
      as.numeric(upstream$contrast)
    )
    expect_identical(
      names(candidate$diagnostics$contrast_vector),
      colnames(stats::model.matrix(case$design, data = fixture$metadata))
    )
    .conformance_expect_unchanged(fixture, before, "edgeR fixture")
  }
})
