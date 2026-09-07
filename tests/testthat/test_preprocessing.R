library(testthat)

test_that("sn_initialize_seurat_object accepts base matrices and adds metadata", {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    c(
      5, 0, 3, 1,
      2, 4, 1, 0,
      0, 1, 0, 2,
      3, 3, 3, 3
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- c("MT-CO1", "RPLP0", "HBB", "CD3D")
  colnames(counts) <- paste0("cell", 1:4)
  metadata <- data.frame(batch = c("a", "a", "b", "b"), row.names = colnames(counts))

  object <- sn_initialize_seurat_object(
    x = counts,
    metadata = metadata,
    project = "prep",
    sample_name = "sample_a",
    study = "study_a",
    species = "human"
  )

  expect_s4_class(object, "Seurat")
  expect_equal(unname(object$sample), rep("sample_a", 4))
  expect_equal(unname(object$study), rep("study_a", 4))
  expect_true(all(c("percent.mt", "percent.ribo", "percent.hb", "batch") %in% colnames(object[[]])))
  expect_true("sn_initialize_seurat_object" %in% names(object@commands))
})

test_that("sn_initialize_seurat_object inherits sample name from named 10x paths", {
  skip_if_not_installed("Seurat")

  root <- tempfile("inputs-")
  dir.create(root)
  outs_dir <- file.path(root, "sampleA", "outs")
  dir.create(file.path(outs_dir, "filtered_feature_bc_matrix"), recursive = TRUE)
  dir.create(file.path(outs_dir, "raw_feature_bc_matrix"), recursive = TRUE)
  writeLines(
    c(
      "%%MatrixMarket matrix coordinate integer general",
      "%",
      "2 2 2",
      "1 1 5",
      "2 2 3"
    ),
    file.path(outs_dir, "filtered_feature_bc_matrix", "matrix.mtx")
  )
  writeLines(c("cell1", "cell2"), file.path(outs_dir, "filtered_feature_bc_matrix", "barcodes.tsv"))
  writeLines(c("gene1\tgene1\tGene Expression", "gene2\tgene2\tGene Expression"), file.path(outs_dir, "filtered_feature_bc_matrix", "features.tsv"))

  tenx_paths <- sn_list_10x_paths(root)
  object <- sn_initialize_seurat_object(
    x = tenx_paths[1],
    project = "prep-named-path",
    species = "human"
  )

  expect_equal(unique(as.character(object$sample)), "sampleA")
  expect_equal(Seurat::Misc(object, "input_source")$sample_name, "sampleA")
})

test_that("sn_initialize_seurat_object accepts multiple 10x paths from sn_list_10x_paths", {
  skip_if_not_installed("Seurat")

  root <- tempfile("inputs-multi-")
  dir.create(root)

  write_tenx_sample <- function(sample_name, entries) {
    outs_dir <- file.path(root, sample_name, "outs")
    dir.create(file.path(outs_dir, "filtered_feature_bc_matrix"), recursive = TRUE)
    writeLines(
      c(
        "%%MatrixMarket matrix coordinate integer general",
        "%",
        "2 2 2",
        entries[[1]],
        entries[[2]]
      ),
      file.path(outs_dir, "filtered_feature_bc_matrix", "matrix.mtx")
    )
    writeLines(
      c("cell1", "cell2"),
      file.path(outs_dir, "filtered_feature_bc_matrix", "barcodes.tsv")
    )
    writeLines(
      c("gene1\tgene1\tGene Expression", "gene2\tgene2\tGene Expression"),
      file.path(outs_dir, "filtered_feature_bc_matrix", "features.tsv")
    )
  }

  write_tenx_sample("sampleA", c("1 1 5", "2 2 3"))
  write_tenx_sample("sampleB", c("1 2 4", "2 1 2"))

  tenx_paths <- sn_list_10x_paths(root)
  objects <- sn_initialize_seurat_object(
    x = tenx_paths,
    project = "prep-multi",
    species = "human"
  )

  expect_type(objects, "list")
  expect_length(objects, 2)
  expect_named(objects, c("sampleA", "sampleB"))
  expect_s4_class(objects[[1]], "Seurat")
  expect_equal(unique(as.character(objects$sampleA$sample)), "sampleA")
  expect_equal(unique(as.character(objects$sampleB$sample)), "sampleB")
})

test_that("sn_initialize_seurat_object infers species and computes QC metrics", {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    c(
      5, 0, 3, 1,
      2, 4, 1, 0,
      0, 1, 0, 2,
      3, 3, 3, 3
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- c("MT-CO1", "RPLP0", "HBB", "CD3D")
  colnames(counts) <- paste0("cell", 1:4)

  object <- sn_initialize_seurat_object(
    x = counts,
    project = "prep-inferred"
  )

  expect_equal(sn_get_species(object), "human")
  expect_true(all(c("percent.mt", "percent.ribo", "percent.hb") %in% colnames(object[[]])))
})

test_that("sn_filter_genes drops genes below the minimum cell threshold", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 1, 1, 1,
        1, 0, 0, 0,
        0, 0, 0, 1,
        2, 2, 2, 2
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", 1:4)
  colnames(counts) <- paste0("cell", 1:4)
  object <- sn_initialize_seurat_object(x = counts, project = "genes")

  filtered <- sn_filter_genes(object, min_cells = 2, plot = FALSE, filter = TRUE)

  expect_equal(rownames(filtered), c("gene1", "gene4"))
})

test_that("sn_filter_genes can use a non-default layer", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      c(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", 1:4)
  colnames(counts) <- paste0("cell", 1:4)

  alt_counts <- Matrix::Matrix(
    matrix(
      c(
        1, 1, 1, 1,
        1, 0, 0, 0,
        0, 0, 0, 1,
        2, 2, 2, 2
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(alt_counts) <- rownames(counts)
  colnames(alt_counts) <- colnames(counts)

  object <- sn_initialize_seurat_object(x = counts, project = "genes-layer")
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- alt_counts

  filtered <- sn_filter_genes(
    object,
    min_cells = 2,
    plot = FALSE,
    filter = TRUE,
    assay = "RNA",
    layer = "decontaminated_counts"
  )

  expect_equal(rownames(filtered), c("gene1", "gene4"))
})

test_that("sn_filter_genes filters only the selected non-default assay", {
  skip_if_not_installed("Seurat")

  rna_counts <- Matrix::Matrix(
    matrix(seq_len(12), nrow = 3, dimnames = list(
      c("CD3D", "MS4A1", "MALAT1"), paste0("cell", 1:4)
    )),
    sparse = TRUE
  )
  adt_counts <- Matrix::Matrix(
    matrix(c(1, 1, 1, 1, 1, 0, 0, 0), nrow = 2, byrow = TRUE,
      dimnames = list(c("CD3", "CD19"), colnames(rna_counts))),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(rna_counts)
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = adt_counts)
  rna_before <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")

  filtered <- sn_filter_genes(
    object,
    min_cells = 2,
    plot = FALSE,
    assay = "ADT",
    layer = "counts"
  )

  expect_identical(rownames(filtered[["ADT"]]), "CD3")
  expect_identical(rownames(filtered[["RNA"]]), rownames(rna_counts))
  expect_equal(
    SeuratObject::LayerData(filtered, assay = "RNA", layer = "counts"),
    rna_before
  )
  expect_identical(SeuratObject::DefaultAssay(filtered), "RNA")
})

test_that("sn_filter_genes can retain only coding genes from bundled annotations", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 1, 1, 1,
        1, 1, 0, 0,
        2, 2, 2, 2,
        3, 0, 0, 0
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- c("CD3D", "MALAT1", "TRBC1", "LINC01409")
  colnames(counts) <- paste0("cell", 1:4)
  object <- sn_initialize_seurat_object(x = counts, project = "genes-coding", species = "human")

  filtered <- sn_filter_genes(
    object,
    min_cells = 1,
    plot = FALSE,
    filter = TRUE,
    gene_class = "coding"
  )

  expect_equal(rownames(filtered), c("CD3D", "TRBC1"))
})

test_that("sn_filter_genes can retain exact GENCODE gene types", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 1, 1, 1,
        1, 0, 0, 0,
        2, 2, 2, 2,
        3, 3, 0, 0
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- c("CD3D", "MALAT1", "TRBC1", "RNU6-1")
  colnames(counts) <- paste0("cell", 1:4)
  object <- sn_initialize_seurat_object(x = counts, project = "genes-type", species = "human")

  filtered <- sn_filter_genes(
    object,
    min_cells = 1,
    plot = FALSE,
    filter = TRUE,
    gene_type = c("lncRNA", "snRNA")
  )

  expect_equal(rownames(filtered), c("MALAT1", "RNU6-1"))
})

test_that("sn_filter_cells adds qc flags and can subset outliers", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(80, lambda = 3), nrow = 20, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:20)
  colnames(counts) <- paste0("cell", 1:4)

  object <- sn_initialize_seurat_object(x = counts, project = "cells")
  object$nCount_RNA[4] <- 1000
  object$nFeature_RNA[4] <- 1000

  flagged <- sn_filter_cells(
    x = object,
    features = c("nCount_RNA", "nFeature_RNA"),
    n = 1,
    plot = FALSE,
    filter = FALSE
  )

  expect_true(all(c("nCount_RNA_qc", "nFeature_RNA_qc") %in% colnames(flagged[[]])))
  expect_equal(unname(flagged$nCount_RNA_qc[4]), "Failed")

  filtered <- sn_filter_cells(
    x = object,
    features = c("nCount_RNA", "nFeature_RNA"),
    n = 1,
    plot = FALSE,
    filter = TRUE
  )

  expect_lt(ncol(filtered), ncol(object))
})

test_that("Seurat-returning helpers record commands in the object history", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(80, lambda = 3), nrow = 20, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:20)
  colnames(counts) <- paste0("cell", 1:4)

  object <- sn_initialize_seurat_object(x = counts, project = "commands")
  object <- sn_filter_genes(object, min_cells = 1, plot = FALSE, filter = FALSE)

  command_names <- names(object@commands)
  expect_true(any(grepl("^sn_initialize_seurat_object", command_names)))
  expect_true(any(grepl("^sn_filter_genes", command_names)))
})

test_that("sn_normalize_data can normalize from a non-default layer", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scran")
  skip_if_not_installed("SingleCellExperiment")

  counts <- Matrix::Matrix(0, nrow = 50, ncol = 60, sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:50)
  colnames(counts) <- paste0("cell", 1:60)

  alt_counts <- Matrix::Matrix(matrix(rpois(50 * 60, lambda = 3), nrow = 50), sparse = TRUE)
  rownames(alt_counts) <- rownames(counts)
  colnames(alt_counts) <- colnames(counts)

  object <- sn_initialize_seurat_object(x = counts, project = "normalize-layer")
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- alt_counts

  normalized <- sn_normalize_data(
    object = object,
    method = "scran",
    assay = "RNA",
    layer = "decontaminated_counts",
    clusters = rep(c("a", "b"), each = 30)
  )

  expect_gt(sum(SeuratObject::LayerData(normalized, assay = "RNA", layer = "data")), 0)
  expect_equal(
    as.matrix(SeuratObject::LayerData(normalized, assay = "RNA", layer = "counts")),
    as.matrix(counts)
  )
})

test_that("sn_normalize_data supports automatic scran clustering from BPCells counts", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("scran")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SeuratObject")

  set.seed(717)
  counts <- Matrix::Matrix(
    matrix(rpois(100 * 120, lambda = 2), nrow = 100),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  path <- tempfile("shennong-scran-bpcells-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  BPCells::write_matrix_dir(
    BPCells::convert_matrix_type(counts, "uint32_t"),
    dir = path
  )
  object <- sn_initialize_seurat_object(
    x = BPCells::open_matrix_dir(path),
    project = "scran-bpcells"
  )

  normalized <- sn_normalize_data(object = object, method = "scran")

  expect_true(Shennong:::.sn_is_iterable_matrix(
    SeuratObject::LayerData(normalized, assay = "RNA", layer = "counts")
  ))
  expect_s4_class(
    SeuratObject::LayerData(normalized, assay = "RNA", layer = "data"),
    "dgCMatrix"
  )
  expect_true(all(is.finite(normalized$size.factor)))
  expect_gt(sum(SeuratObject::LayerData(normalized, assay = "RNA", layer = "data")), 0)
})

test_that("sn_normalize_data SCTransform restores future globals option", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glmGamPoi")

  counts <- Matrix::Matrix(matrix(rpois(60 * 40, lambda = 4), nrow = 60), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- sn_initialize_seurat_object(x = counts, project = "normalize-sct")

  old_future <- getOption("future.globals.maxSize", NULL)
  on.exit(options(future.globals.maxSize = old_future), add = TRUE)
  options(future.globals.maxSize = 1024)

  normalized <- sn_normalize_data(
    object = object,
    method = "sctransform",
    assay = "RNA",
    variable.features.n = 30,
    verbose = FALSE
  )

  expect_true("SCT" %in% names(normalized@assays))
  expect_equal(getOption("future.globals.maxSize"), 1024)
})

test_that("sn_get_species reads stored species and can infer from Seurat object features", {
  skip_if_not_installed("Seurat")

  mouse_counts <- Matrix::Matrix(
    matrix(1, nrow = 3, ncol = 2, dimnames = list(c("Cd3d", "Ltb", "Ms4a1"), c("cell1", "cell2"))),
    sparse = TRUE
  )
  mouse_object <- SeuratObject::CreateSeuratObject(counts = mouse_counts)
  Seurat::Misc(mouse_object, slot = "species") <- "mouse"

  human_counts <- Matrix::Matrix(
    matrix(
      c(
        5, 0, 3, 1,
        2, 4, 1, 0,
        0, 1, 0, 2,
        3, 3, 3, 3
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(c("MT-CO1", "RPLP0", "HBB", "CD3D"), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )
  human_object <- SeuratObject::CreateSeuratObject(counts = human_counts)

  expect_equal(sn_get_species(mouse_object), "mouse")
  expect_equal(sn_get_species(human_object), "human")
})

test_that("sn_standardize_gene_symbols converts gene IDs and aggregates duplicates", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("HGNChelper")
  skip_if_not_installed("dplyr")

  annotations <- Shennong:::.sn_get_gene_annotation_table("human")
  dup_name <- annotations$gene_name[duplicated(annotations$gene_name)][[1]]
  dup_rows <- annotations[annotations$gene_name == dup_name, , drop = FALSE]
  unique_row <- annotations[
    !duplicated(annotations$gene_name) &
      !is.na(annotations$gene_name) &
      nzchar(annotations$gene_name),
    ,
    drop = FALSE
  ][1, , drop = FALSE]

  counts <- Matrix::Matrix(
    matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE),
    sparse = TRUE
  )
  rownames(counts) <- c(dup_rows$gene_id[[1]], dup_rows$gene_id[[2]], unique_row$gene_id[[1]])
  colnames(counts) <- c("cell1", "cell2")

  standardized <- suppressWarnings(
    sn_standardize_gene_symbols(
      counts,
      species = "human",
      is_gene_id = TRUE
    )
  )

  expect_equal(nrow(standardized), 2)
  expect_true(dup_name %in% rownames(standardized))
  expect_true(unique_row$gene_name[[1]] %in% rownames(standardized))
  expect_equal(as.numeric(standardized[dup_name, ]), c(4, 6))

  object <- sn_initialize_seurat_object(x = counts, project = "gene-standardize", species = "human")
  standardized_object <- suppressWarnings(sn_standardize_gene_symbols(
    object,
    species = "human",
    is_gene_id = TRUE
  ))

  expect_true(dup_name %in% rownames(standardized_object))
  expect_true("sn_standardize_gene_symbols" %in% names(standardized_object@commands))
})

test_that("sn_standardize_gene_symbols accepts character vectors", {
  skip_if_not_installed("HGNChelper")
  skip_if_not_installed("dplyr")

  standardized <- suppressWarnings(sn_standardize_gene_symbols(
    c("CD3D", "fakegene1", "1-Mar", ""),
    species = "human",
    is_gene_id = FALSE
  ))

  expect_type(standardized, "character")
  expect_equal(standardized, c("CD3D", "fakegene1", "1-Mar"))
  expect_false(any(is.na(standardized)))
  expect_false(any(!nzchar(standardized)))
})

test_that("sn_standardize_gene_symbols converts gene ID vectors", {
  skip_if_not_installed("HGNChelper")
  skip_if_not_installed("dplyr")

  annotations <- Shennong:::.sn_get_gene_annotation_table("human")
  dup_name <- annotations$gene_name[duplicated(annotations$gene_name)][[1]]
  dup_rows <- annotations[annotations$gene_name == dup_name, , drop = FALSE]
  unique_row <- annotations[
    !duplicated(annotations$gene_name) &
      !is.na(annotations$gene_name) &
      nzchar(annotations$gene_name),
    ,
    drop = FALSE
  ][1, , drop = FALSE]

  standardized <- suppressWarnings(sn_standardize_gene_symbols(
    c(dup_rows$gene_id[[1]], dup_rows$gene_id[[2]], unique_row$gene_id[[1]], "missing-id"),
    species = "human",
    is_gene_id = TRUE
  ))

  expect_equal(standardized, c(dup_name, dup_name, unique_row$gene_name[[1]]))
})

test_that("sn_standardize_gene_symbols preserves unresolved symbols without returning NA row names", {
  skip_if_not_installed("HGNChelper")
  skip_if_not_installed("dplyr")

  counts <- Matrix::Matrix(
    matrix(seq_len(8), nrow = 4),
    sparse = TRUE
  )
  rownames(counts) <- c("CD3D", "fakegene1", "1-Mar", "")
  colnames(counts) <- c("cell1", "cell2")

  standardized <- suppressWarnings(sn_standardize_gene_symbols(
    counts,
    species = "human",
    is_gene_id = FALSE
  ))

  expect_false(any(is.na(rownames(standardized))))
  expect_false(any(!nzchar(rownames(standardized))))
  expect_true("fakegene1" %in% rownames(standardized))
  expect_true("1-Mar" %in% rownames(standardized))
})

test_that("sn_standardize_gene_symbols targets RNA and preserves a non-default assay", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("HGNChelper")
  skip_if_not_installed("dplyr")

  rna_counts <- Matrix::Matrix(
    matrix(seq_len(12), nrow = 3, dimnames = list(
      c("CD3D", "MS4A1", "MALAT1"), paste0("cell", 1:4)
    )),
    sparse = TRUE
  )
  adt_counts <- Matrix::Matrix(
    matrix(seq_len(8), nrow = 2, dimnames = list(
      c("CD3", "CD19"), colnames(rna_counts)
    )),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(rna_counts)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  SeuratObject::VariableFeatures(object[["RNA"]]) <- c("CD3D", "MS4A1")
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = adt_counts)
  SeuratObject::DefaultAssay(object) <- "ADT"
  rna_data_before <- SeuratObject::LayerData(object, assay = "RNA", layer = "data")
  adt_before <- SeuratObject::LayerData(object, assay = "ADT", layer = "counts")

  standardized <- suppressWarnings(sn_standardize_gene_symbols(
    object,
    species = "human"
  ))

  expect_identical(SeuratObject::DefaultAssay(standardized), "ADT")
  expect_identical(rownames(standardized[["RNA"]]), rownames(rna_counts))
  expect_equal(
    SeuratObject::LayerData(standardized, assay = "RNA", layer = "data"),
    rna_data_before
  )
  expect_equal(
    SeuratObject::LayerData(standardized, assay = "ADT", layer = "counts"),
    adt_before
  )
  expect_identical(
    SeuratObject::VariableFeatures(standardized[["RNA"]]),
    c("CD3D", "MS4A1")
  )
})

test_that("gene-symbol standardization preserves assay provenance and invalidates RNA-derived state", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("HGNChelper")
  skip_if_not_installed("dplyr")

  set.seed(817L)
  counts <- matrix(
    stats::rpois(20L * 24L, lambda = 4),
    nrow = 20L,
    dimnames = list(c("CD3D", "MS4A1", "MALAT1", paste0("GENE", 4:20)), paste0("cell", 1:24))
  )
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  object[["RNA"]]@assay.orig <- "source_rna"
  object[["RNA"]]@misc <- list(source = "retained provenance")
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- suppressWarnings(Seurat::FindVariableFeatures(object, nfeatures = 10L, verbose = FALSE))
  object <- Seurat::ScaleData(object, verbose = FALSE)
  object <- suppressWarnings(Seurat::RunPCA(object, npcs = 5L, verbose = FALSE))
  object <- Seurat::FindNeighbors(object, dims = 1:5, verbose = FALSE)
  old_rna_graphs <- names(object@graphs)

  adt <- Matrix::Matrix(
    matrix(stats::rpois(2L * ncol(object), 2), nrow = 2L,
      dimnames = list(c("CD3", "CD19"), colnames(object))),
    sparse = TRUE
  )
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = adt)
  adt_embedding <- matrix(
    stats::rnorm(ncol(object) * 2L),
    ncol = 2L,
    dimnames = list(colnames(object), c("ADT_1", "ADT_2"))
  )
  object[["adt_dr"]] <- SeuratObject::CreateDimReducObject(
    embeddings = adt_embedding,
    key = "ADT_",
    assay = "ADT"
  )
  object <- Shennong:::.sn_log_seurat_command(
    object,
    assay = "ADT",
    name = "adt_command"
  )

  standardized <- suppressWarnings(sn_standardize_gene_symbols(
    object,
    species = "human"
  ))

  expect_identical(standardized[["RNA"]]@assay.orig, "source_rna")
  expect_identical(standardized[["RNA"]]@misc, list(source = "retained provenance"))
  expect_false("pca" %in% names(standardized@reductions))
  expect_false(any(old_rna_graphs %in% names(standardized@graphs)))
  expect_true("adt_dr" %in% names(standardized@reductions))
  expect_true("adt_command" %in% names(standardized@commands))
  expect_true("sn_standardize_gene_symbols" %in% names(standardized@commands))
  expect_false(any(vapply(
    standardized@commands[setdiff(names(standardized@commands), "sn_standardize_gene_symbols")],
    function(command) identical(methods::slot(command, "assay.used"), "RNA"),
    logical(1)
  )))
})

test_that("ambient correction writes back to the explicitly selected assay", {
  skip_if_not_installed("Seurat")

  rna_counts <- Matrix::Matrix(
    matrix(seq_len(12), nrow = 3, dimnames = list(
      c("CD3D", "MS4A1", "MALAT1"), paste0("cell", 1:4)
    )),
    sparse = TRUE
  )
  adt_counts <- Matrix::Matrix(
    matrix(seq_len(8), nrow = 2, dimnames = list(
      c("CD3", "CD19"), colnames(rna_counts)
    )),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(rna_counts)
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = adt_counts)
  corrected <- adt_counts
  corrected[] <- pmax(corrected[] - 1, 0)

  updated <- Shennong:::.sn_apply_ambient_result_to_object(
    object = object,
    out = list(
      counts = corrected,
      metadata = NULL,
      zero_cells = character(),
      removed_cells = character()
    ),
    assay = "ADT",
    layer = "decontaminated_counts"
  )

  expect_true("decontaminated_counts" %in% SeuratObject::Layers(updated[["ADT"]]))
  expect_false("decontaminated_counts" %in% SeuratObject::Layers(updated[["RNA"]]))
  expect_equal(
    SeuratObject::LayerData(updated, assay = "ADT", layer = "decontaminated_counts"),
    corrected
  )
})

test_that("SoupX estimates its soup profile from raw droplets", {
  skip_if_not_installed("SoupX")

  filtered <- Matrix::Matrix(
    matrix(c(9, 1, 9, 1), nrow = 2L,
      dimnames = list(c("G1", "G2"), c("cell1", "cell2"))),
    sparse = TRUE
  )
  raw <- Matrix::Matrix(
    matrix(c(1, 9, 1, 9), nrow = 2L,
      dimnames = list(c("G1", "G2"), c("drop1", "drop2"))),
    sparse = TRUE
  )
  captured <- new.env(parent = emptyenv())

  out <- with_mocked_bindings(
    with_mocked_bindings(
      Shennong:::.sn_remove_ambient_soupx(
        x_info = list(object = NULL, counts = filtered),
        raw_info = list(object = NULL, counts = raw),
        cluster = c(cell1 = "a", cell2 = "b"),
        seed = 19L
      ),
      SoupChannel = function(tod, toc, ...) list(tod = tod, toc = toc),
      setSoupProfile = function(sc, soupProfile) {
        captured$profile <- soupProfile
        sc
      },
      setClusters = function(sc, clusters) sc,
      autoEstCont = function(sc, ...) sc,
      .package = "SoupX"
    ),
    .sn_adjust_soupx_counts = function(sc, seed) {
      captured$seed <- seed
      sc$toc
    },
    .package = "Shennong"
  )

  raw_sums <- Matrix::rowSums(raw)
  expect_equal(captured$profile$counts, as.numeric(raw_sums))
  expect_equal(captured$profile$est, as.numeric(raw_sums / sum(raw_sums)))
  expect_identical(captured$seed, 19L)
  expect_equal(out$counts, filtered)
})

test_that("SoupX stochastic rounding is local and reproducible", {
  expect_identical(eval(formals(sn_remove_ambient_contamination)$seed), 717L)
  counts <- Matrix::Matrix(
    matrix(c(0.2, 1.8, 2.5, 3.1), nrow = 2L),
    sparse = TRUE
  )
  set.seed(1201L)
  rng_before <- .Random.seed
  first <- Shennong:::.sn_round_soupx_counts(counts, seed = 31L)
  expect_identical(.Random.seed, rng_before)
  second <- Shennong:::.sn_round_soupx_counts(counts, seed = 31L)
  expect_identical(first, second)
  expect_true(all(first@x == floor(first@x)))
})

test_that("decontX and decontPro retain filtered-only features", {
  skip_if_not_installed("decontX")

  filtered <- Matrix::Matrix(
    matrix(c(5, 4, 7, 5, 4, 7), nrow = 3L,
      dimnames = list(c("shared1", "shared2", "filtered_only"), c("cell1", "cell2"))),
    sparse = TRUE
  )
  raw <- Matrix::Matrix(
    matrix(c(9, 8, 3, 9, 8, 3), nrow = 3L,
      dimnames = list(c("shared1", "shared2", "raw_only"), c("drop1", "drop2"))),
    sparse = TRUE
  )

  decontx <- with_mocked_bindings(
    Shennong:::.sn_remove_ambient_decontx(
      x_info = list(object = NULL, counts = filtered),
      raw_info = list(object = NULL, counts = raw),
      cluster_backend = "native"
    ),
    decontX = function(x, background, ...) list(
      decontXcounts = x - 1,
      contamination = rep(0.1, ncol(x)),
      z = rep(1L, ncol(x))
    ),
    .package = "decontX"
  )
  expect_identical(rownames(decontx$counts), rownames(filtered))
  expect_equal(
    as.numeric(decontx$counts["filtered_only", ]),
    as.numeric(filtered["filtered_only", ])
  )

  decontpro <- with_mocked_bindings(
    Shennong:::.sn_remove_ambient_decontpro(
      x_info = list(object = NULL, counts = filtered),
      raw_info = list(object = NULL, counts = raw),
      cluster = c("a", "b")
    ),
    decontPro = function(filtered_counts, cell_type, ambient_counts, ...) list(
      decontaminated_counts = filtered_counts - 1,
      ambient_counts = filtered_counts * 0,
      background_counts = filtered_counts * 0
    ),
    .package = "decontX"
  )
  expect_identical(rownames(decontpro$counts), rownames(filtered))
  expect_equal(
    as.numeric(decontpro$counts["filtered_only", ]),
    as.numeric(filtered["filtered_only", ])
  )
})

test_that("sn_score_cell_cycle returns the object unchanged when markers do not overlap", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(1, nrow = 5, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:5)
  colnames(counts) <- paste0("cell", 1:4)

  object <- sn_initialize_seurat_object(x = counts, project = "cell-cycle-skip", species = "human")
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  updated <- sn_score_cell_cycle(object, species = "human")

  expect_false(any(c("S.Score", "G2M.Score", "Phase") %in% colnames(updated[[]])))
})

test_that("sn_score_cell_cycle scores an explicit assay and custom layer", {
  skip_if_not_installed("Seurat")

  s_features <- sn_get_signatures(species = "human", category = "Programs/cellCycle.G1S")
  g2m_features <- sn_get_signatures(species = "human", category = "Programs/cellCycle.G2M")
  features <- unique(c(s_features, g2m_features, paste0("BACKGROUND", seq_len(3000))))
  cells <- paste0("cell", seq_len(18))
  set.seed(717)
  counts <- Matrix::Matrix(
    matrix(rpois(length(features) * length(cells), lambda = 3), nrow = length(features),
      dimnames = list(features, cells)),
    sparse = TRUE
  )

  object <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA")
  object[["cycle"]] <- SeuratObject::CreateAssay5Object(counts = counts)
  object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
  object <- Seurat::NormalizeData(object, assay = "cycle", verbose = FALSE)
  original_cycle_data <- SeuratObject::LayerData(object, assay = "cycle", layer = "data")
  selected_data <- as.matrix(original_cycle_data)
  selected_data[s_features, seq_len(6)] <- selected_data[s_features, seq_len(6)] + 5
  selected_data[g2m_features, 7:12] <- selected_data[g2m_features, 7:12] + 5
  selected_data <- Matrix::Matrix(selected_data, sparse = TRUE)
  SeuratObject::LayerData(object, assay = "cycle", layer = "cycle_input") <- selected_data

  expected <- object
  SeuratObject::LayerData(expected, assay = "cycle", layer = "data") <- selected_data
  expected <- Seurat::CellCycleScoring(
    expected,
    s.features = s_features,
    g2m.features = g2m_features,
    assay = "cycle",
    slot = "data"
  )

  updated <- sn_score_cell_cycle(
    object,
    species = "human",
    assay = "cycle",
    layer = "cycle_input"
  )

  expect_equal(updated$S.Score, expected$S.Score)
  expect_equal(updated$G2M.Score, expected$G2M.Score)
  expect_identical(updated$Phase, expected$Phase)
  expect_equal(updated$CC.Difference, expected$S.Score - expected$G2M.Score)
  expect_identical(SeuratObject::DefaultAssay(updated), "RNA")
  expect_equal(
    SeuratObject::LayerData(updated, assay = "cycle", layer = "data"),
    original_cycle_data
  )
  expect_equal(
    SeuratObject::LayerData(updated, assay = "cycle", layer = "cycle_input"),
    selected_data
  )
  expect_error(
    sn_score_cell_cycle(object, species = "human", assay = "cycle", layer = "missing"),
    "Layer 'missing' was not found"
  )
})

test_that("count-input resolution supports Seurat objects, paths, and matrix-like inputs", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(1:12, nrow = 3, dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4))),
    sparse = TRUE
  )
  object <- sn_initialize_seurat_object(x = counts, project = "resolve-counts")

  from_object <- Shennong:::.sn_resolve_counts_input(object)
  from_matrix <- Shennong:::.sn_resolve_counts_input(as.matrix(counts))

  expect_s4_class(from_object$object, "Seurat")
  expect_s4_class(from_object$counts, "dgCMatrix")
  expect_equal(as.matrix(from_object$counts), as.matrix(counts))
  expect_null(from_matrix$object)
  expect_s4_class(from_matrix$counts, "dgCMatrix")
  expect_equal(as.matrix(from_matrix$counts), as.matrix(counts))

  if (rlang::is_installed("rio")) {
    path <- tempfile(fileext = ".csv")
    utils::write.csv(as.data.frame(as.matrix(counts)), path, row.names = FALSE)
    from_path <- Shennong:::.sn_resolve_counts_input(path)
    expect_null(from_path$object)
    expect_s4_class(from_path$counts, "dgCMatrix")
    expect_identical(dim(from_path$counts), dim(counts))
  }

  expect_error(
    Shennong:::.sn_resolve_counts_input(list(counts = counts), arg = "x"),
    "`x` must be a Seurat object, matrix-like object, or path"
  )
})

test_that("count-shape restoration preserves the original matrix extent", {
  original_counts <- Matrix::Matrix(
    matrix(1:9, nrow = 3, dimnames = list(paste0("gene", 1:3), paste0("cell", 1:3))),
    sparse = TRUE
  )
  corrected_counts <- Matrix::Matrix(
    matrix(c(10, 20, 30, 40), nrow = 2, dimnames = list(c("gene1", "gene3"), c("cell1", "cell3"))),
    sparse = TRUE
  )
  expected <- original_counts
  expected[rownames(corrected_counts), colnames(corrected_counts)] <- corrected_counts

  restored_subset <- Shennong:::.sn_restore_count_shape(original_counts, corrected_counts)
  restored_same_shape <- Shennong:::.sn_restore_count_shape(original_counts, original_counts)

  expect_equal(as.matrix(restored_subset), as.matrix(expected))
  expect_identical(dimnames(restored_subset), dimnames(original_counts))
  expect_identical(restored_same_shape, original_counts)
})

test_that("10x metric parsing and source detection handle empty and h5-backed outputs", {
  empty_metrics <- tempfile(fileext = ".csv")
  file.create(empty_metrics)

  expect_null(Shennong:::.sn_parse_10x_metrics_summary(empty_metrics))
  expect_null(Shennong:::.sn_parse_10x_metrics_summary(tempfile(fileext = ".csv")))

  sample_root <- tempfile("tenx-h5-")
  outs_dir <- file.path(sample_root, "outs")
  dir.create(outs_dir, recursive = TRUE)
  filtered_h5 <- file.path(outs_dir, "filtered_feature_bc_matrix.h5")
  raw_h5 <- file.path(outs_dir, "raw_feature_bc_matrix.h5")
  web_summary <- file.path(outs_dir, "web_summary.html")
  metrics_csv <- file.path(outs_dir, "metrics_summary.csv")
  file.create(filtered_h5)
  file.create(raw_h5)
  writeLines("<html></html>", web_summary)
  utils::write.csv(data.frame(metric = "cells", value = 10), metrics_csv, row.names = FALSE)

  source_info <- Shennong:::.sn_detect_10x_outs_source(outs_dir)

  expect_equal(source_info$type, "10x_outs")
  expect_equal(normalizePath(source_info$filtered_path, winslash = "/", mustWork = TRUE), normalizePath(filtered_h5, winslash = "/", mustWork = TRUE))
  expect_equal(normalizePath(source_info$raw_path, winslash = "/", mustWork = TRUE), normalizePath(raw_h5, winslash = "/", mustWork = TRUE))
  expect_equal(normalizePath(source_info$web_summary_path, winslash = "/", mustWork = TRUE), normalizePath(web_summary, winslash = "/", mustWork = TRUE))
  expect_s3_class(source_info$metrics, "data.frame")
})

test_that("sn_initialize_seurat_object validates metadata shape for multi-path imports", {
  expect_error(sn_initialize_seurat_object(
    x = c(sample1 = "a", sample2 = "b"),
    metadata = data.frame(a = 1)
  ))
})

test_that("ambient-cluster resolution supports metadata columns and explicit vectors", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(24, lambda = 2), nrow = 6, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:6)
  colnames(counts) <- paste0("cell", 1:4)
  object <- sn_initialize_seurat_object(x = counts, project = "ambient-clusters")
  object$cluster_id <- rep(c("a", "b"), each = 2)
  x_info <- list(
    object = object,
    counts = SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  )

  resolved_from_column <- Shennong:::.sn_resolve_ambient_clusters(x_info, cluster = "cluster_id")
  resolved_from_vector <- Shennong:::.sn_resolve_ambient_clusters(
    x_info,
    cluster = c("x", "x", "y", "y")
  )

  expect_equal(unname(resolved_from_column), as.character(object$cluster_id))
  expect_equal(names(resolved_from_column), colnames(object))
  expect_equal(unname(resolved_from_vector), c("x", "x", "y", "y"))
  expect_equal(names(resolved_from_vector), colnames(object))

  expect_error(
    Shennong:::.sn_resolve_ambient_clusters(x_info, cluster = "missing_cluster"),
    "Cluster column 'missing_cluster' was not found"
  )
  expect_error(
    Shennong:::.sn_resolve_ambient_clusters(x_info, cluster = c("a", "b")),
    "must have one value per cell"
  )
})

test_that("ambient cluster backends use method-specific defaults and explicit overrides", {
  expect_identical(
    Shennong:::.sn_resolve_ambient_cluster_backend("decontx"),
    "native"
  )
  expect_identical(
    Shennong:::.sn_resolve_ambient_cluster_backend("soupx"),
    "shennong"
  )
  expect_identical(
    Shennong:::.sn_resolve_ambient_cluster_backend("decontpro"),
    "shennong"
  )
  expect_identical(
    Shennong:::.sn_resolve_ambient_cluster_backend(
      "decontx",
      cluster_backend = "shennong"
    ),
    "shennong"
  )
  expect_identical(
    Shennong:::.sn_resolve_ambient_cluster_backend(
      "decontx",
      cluster_backend = "native",
      cluster = c("a", "b")
    ),
    "provided"
  )
  expect_error(
    Shennong:::.sn_resolve_ambient_cluster_backend(
      "soupx",
      cluster_backend = "native"
    ),
    "SoupX does not provide native clustering"
  )
  expect_error(
    Shennong:::.sn_resolve_ambient_cluster_backend(
      "decontpro",
      cluster_backend = "native"
    ),
    "decontPro does not provide native clustering"
  )
})

test_that("sn_filter_cells stores grouped QC thresholds in misc metadata", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(36, lambda = 3), nrow = 6, ncol = 6), sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:6)
  colnames(counts) <- paste0("cell", 1:6)

  object <- sn_initialize_seurat_object(x = counts, project = "cells-grouped")
  object$sample <- rep(c("A", "B"), each = 3)
  object$nCount_RNA <- c(10, 11, 30, 20, 21, 60)

  flagged <- sn_filter_cells(
    x = object,
    features = "nCount_RNA",
    group_by = "sample",
    n = 1,
    plot = FALSE,
    filter = FALSE
  )

  qc_info <- Seurat::Misc(flagged, "qc")$nCount_RNA

  expect_equal(nrow(qc_info), 2)
  expect_true(all(c("sample", "median", "mad", "l", "u") %in% colnames(qc_info)))
  expect_true("nCount_RNA_qc" %in% colnames(flagged[[]]))
  expect_true(all(flagged$nCount_RNA_qc %in% c("Passed", "Failed")))
})

test_that("sn_filter_cells keeps constant-valued cells when MAD is zero", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(1, nrow = 6, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:6)
  colnames(counts) <- paste0("cell", 1:4)

  object <- sn_initialize_seurat_object(x = counts, project = "cells-constant")
  object$nCount_RNA <- rep(100, 4)

  flagged <- sn_filter_cells(
    x = object,
    features = "nCount_RNA",
    n = 1,
    plot = FALSE,
    filter = FALSE
  )

  expect_true(all(flagged$nCount_RNA_qc == "Passed"))
})

test_that("sn_filter_cells validates the filtering method", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(1, nrow = 6, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:6)
  colnames(counts) <- paste0("cell", 1:4)

  object <- sn_initialize_seurat_object(x = counts, project = "cells-method")

  expect_error(
    sn_filter_cells(
      x = object,
      features = "nCount_RNA",
      method = "iqr",
      plot = FALSE
    ),
    "must be one of"
  )
})

test_that("sn_initialize_seurat_object records 10x outs source metadata", {
  skip_if_not_installed("Seurat")

  sample_dir <- tempfile("tenx-sample-")
  outs_dir <- file.path(sample_dir, "outs")
  filtered_dir <- file.path(outs_dir, "filtered_feature_bc_matrix")
  raw_dir <- file.path(outs_dir, "raw_feature_bc_matrix")
  dir.create(filtered_dir, recursive = TRUE)
  dir.create(raw_dir, recursive = TRUE)
  utils::write.csv(
    data.frame(
      `Estimated Number of Cells` = 1234,
      `Median Genes per Cell` = 567,
      check.names = FALSE
    ),
    file.path(outs_dir, "metrics_summary.csv"),
    row.names = FALSE
  )

  counts <- Matrix::Matrix(
    matrix(1:12, nrow = 3, dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4))),
    sparse = TRUE
  )

  object <- with_mocked_bindings(
    sn_initialize_seurat_object(x = outs_dir, project = "tenx-outs"),
    sn_read = function(path, ...) {
      expect_equal(
        normalizePath(path, winslash = "/", mustWork = TRUE),
        normalizePath(filtered_dir, winslash = "/", mustWork = TRUE)
      )
      counts
    },
    .package = "Shennong"
  )

  source_info <- Seurat::Misc(object, "input_source")
  expect_equal(source_info$type, "10x_outs")
  expect_equal(
    normalizePath(source_info$outs_path, winslash = "/", mustWork = TRUE),
    normalizePath(outs_dir, winslash = "/", mustWork = TRUE)
  )
  expect_equal(
    normalizePath(source_info$filtered_path, winslash = "/", mustWork = TRUE),
    normalizePath(filtered_dir, winslash = "/", mustWork = TRUE)
  )
  expect_equal(
    normalizePath(source_info$raw_path, winslash = "/", mustWork = TRUE),
    normalizePath(raw_dir, winslash = "/", mustWork = TRUE)
  )
  expect_s3_class(source_info$metrics, "data.frame")
  expect_equal(source_info$metrics$`Estimated Number of Cells`, 1234)
})

test_that("sn_filter_genes accepts large plotting thresholds without NA counts", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 1, 0, 0,
        1, 0, 0, 0,
        0, 0, 0, 1
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", 1:3)
  colnames(counts) <- paste0("cell", 1:4)
  object <- sn_initialize_seurat_object(x = counts, project = "genes-plot")

  expect_no_error(
    sn_filter_genes(
      object,
      min_cells = 10,
      plot = TRUE,
      filter = FALSE
    )
  )
})
