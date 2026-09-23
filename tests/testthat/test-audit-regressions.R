audit_de_object <- function() {
  set.seed(812)
  counts <- matrix(0, 100, 80, dimnames = list(paste0("g", 1:100), paste0("c", 1:80)))
  counts[1:10, ] <- rpois(800, 3)
  counts[1:5, rep(c(TRUE, FALSE), each = 20, times = 2)] <-
    counts[1:5, rep(c(TRUE, FALSE), each = 20, times = 2)] + 6
  object <- Seurat::NormalizeData(SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE)), verbose = FALSE)
  object$group <- rep(c("A", "B"), each = 20, times = 2)
  object$stratum <- rep(c("one", "two"), each = 40)
  object
}

test_that("DE stores post-filter backgrounds and subset markers preserve gene IDs (#9, #11)", {
  object <- audit_de_object()
  result <- sn_find_de(object, analysis = "contrast", ident_1 = "A", ident_2 = "B",
    group_by = "group", min_pct = .25, logfc_threshold = 0, verbose = FALSE)
  stored <- sn_get_result(result, "de", "default")
  expect_setequal(stored$input$tested_features, paste0("g", 1:10))
  markers <- sn_find_de(object, analysis = "markers", group_by = "group",
    subset_by = "stratum", min_pct = .25, logfc_threshold = 0, verbose = FALSE)
  stored_markers <- sn_get_result(markers, "de", "default")
  expect_setequal(stored_markers$tables$primary$stratum, c("one", "two"))
  expect_true(all(stored_markers$tables$primary$gene %in% rownames(object)))
  expect_false(anyDuplicated(names(stored_markers$tables$primary)) > 0L)
  expect_true(all(stored_markers$input$tested_features %in% paste0("g", 1:10)))
})

test_that("NicheNet split layers agree with a joined public call (#8)", {
  skip_if_not_installed("nichenetr")
  set.seed(717)
  m <- matrix(rpois(30 * 20, 4), 30,
    dimnames = list(paste0("g", 1:30), paste0("c", 1:20)))
  m[1, 1:10] <- 0
  object <- Seurat::NormalizeData(SeuratObject::CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE)), verbose = FALSE)
  object$batch <- rep(c("b1", "b2"), each = 10)
  object$group <- rep(c("T", "B"), 10)
  target <- matrix(runif(60), 30, dimnames = list(rownames(object), c("g1", "g2")))
  args <- list(method = "nichenet", group_by = "group", species = "human",
    sender = "T", receiver = "B", geneset = paste0("g", 5:15),
    background_genes = paste0("g", 5:30), ligand_target_matrix = target,
    lr_network = data.frame(ligand = c("g1", "g2"), receptor = c("g3", "g4")),
    expressed_pct = .1, return_object = FALSE)
  joined <- do.call(sn_run_cell_communication, c(list(object = object), args))
  object[["RNA"]] <- split(object[["RNA"]], f = object$batch)
  before <- SeuratObject::Layers(object[["RNA"]])
  split_result <- do.call(sn_run_cell_communication, c(list(object = object), args))
  expect_equal(split_result$tables$primary, joined$tables$primary)
  expect_identical(SeuratObject::Layers(object[["RNA"]]), before)
})

test_that("Milo annotation filtering works without an FDR filter (#12)", {
  object <- audit_de_object()
  table <- data.frame(Nhood = 1:2, logFC = c(1, -1), PValue = c(.01, .02),
    FDR = c(.02, .03), SpatialFDR = c(.03, .04), cell_type = c("Tcell", "Bcell"))
  object <- sn_store_milo(object, table, sample_by = "sample", group_by = "group", annotation_by = "cell_type")
  expect_identical(sn_get_milo_result(object, annotation = "Tcell")$cell_type, "Tcell")
  expect_equal(nrow(sn_get_milo_result(object, annotation = "absent")), 0L)
})

test_that("native matrix serialization and BPCells writers preserve matrices (#13)", {
  m <- matrix(c(1, 0, 2, 3), 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))
  for (format in c("rds", "qs2")) {
    if (format == "qs2" && !requireNamespace("qs2", quietly = TRUE)) next
    path <- tempfile(fileext = paste0(".", format))
    sn_write(m, path, auto_install = FALSE)
    actual <- if (format == "rds") readRDS(path) else qs2::qs_read(path)
    expect_identical(actual, m)
    unlink(path)
  }
  skip_if_not_installed("BPCells")
  for (matrix in list(m, methods::as(Matrix::Matrix(m, sparse = TRUE), "dgCMatrix"))) {
    for (format in c("h5", "bpcells")) {
      path <- tempfile(fileext = paste0(".", format))
      suppressWarnings(sn_write(matrix, path, auto_install = FALSE))
      readback <- if (format == "h5") BPCells::open_matrix_10x_hdf5(path) else BPCells::open_matrix_dir(path)
      expect_equal(as.matrix(readback), m)
      unlink(path, recursive = TRUE)
    }
  }
})

test_that("program comparisons recognize complete donor pairs (#10)", {
  baseline <- seq(10, 60, 10)
  values <- as.vector(rbind(baseline, baseline + c(1, 2, 3, 1, 2, 3)))
  m <- rbind(g1 = values, g2 = rep(1, 12)); colnames(m) <- paste0("c", 1:12)
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE))
  object$donor <- rep(paste0("d", 1:6), each = 2)
  object$condition <- rep(c("A", "B"), 6)
  object <- sn_score_programs(object, list(P = "g1"), method = "mean", layer = "counts", result_id = "s")
  actual <- sn_test_programs(object, "s", "condition", sample_by = "donor",
    contrast = c("B", "A"), return_object = FALSE)
  expected <- stats::wilcox.test(values[seq(2, 12, 2)], values[seq(1, 11, 2)],
    paired = TRUE, exact = FALSE)$p.value
  expect_equal(actual$tables$primary$p_value, expected)
  expect_true(actual$tables$primary$paired)
  if (requireNamespace("limma", quietly = TRUE)) {
    limma_result <- sn_test_programs(object, "s", "condition", sample_by = "donor",
      contrast = c("B", "A"), method = "limma", return_object = FALSE)
    condition <- factor(object$condition, levels = c("B", "A"))
    donor <- factor(object$donor)
    design <- stats::model.matrix(~0 + condition + donor)
    fit <- limma::lmFit(matrix(values, nrow = 1), design)
    fit <- limma::eBayes(limma::contrasts.fit(fit, matrix(c(1, -1, rep(0, 5)), ncol = 1)))
    expect_equal(limma_result$tables$primary$p_value, unname(fit$p.value[[1]]))
  }
  incomplete <- object[, -12]
  expect_error(sn_test_programs(incomplete, "s", "condition", sample_by = "donor"), "paired|missing|absent")
})

test_that("program metadata names do not overwrite programs or user columns (#14)", {
  object <- audit_de_object()
  object$score_T_cell <- 42
  object <- sn_score_programs(object, list(`T-cell` = "g1", T_cell = "g2"),
    method = "mean", layer = "counts", result_id = "score")
  expect_true(all(object$score_T_cell == 42))
  stored <- sn_get_result(object, "program_scoring", "score")
  mapping <- stored$tables$metadata_columns
  expect_setequal(mapping$program, c("T-cell", "T_cell"))
  expect_equal(length(unique(mapping$column)), 2L)
  for (i in seq_len(nrow(mapping))) {
    expected <- stored$tables$scores$score[stored$tables$scores$program == mapping$program[i]]
    expect_equal(unname(object[[mapping$column[i], drop = TRUE]]), expected)
  }
  again <- sn_score_programs(object, list(`T-cell` = "g1", T_cell = "g2"),
    method = "mean", layer = "counts", result_id = "score", overwrite = TRUE)
  expect_identical(sn_get_result(again, "program_scoring", "score")$tables$metadata_columns, mapping)
  other <- sn_score_programs(again, list(cell = "g3"), method = "mean", layer = "counts", result_id = "score-T")
  other <- sn_score_programs(other, list(cell = "g4"), method = "mean", layer = "counts", result_id = "score_T")
  first <- sn_get_result(other, "program_scoring", "score-T")$tables$metadata_columns$column
  second <- sn_get_result(other, "program_scoring", "score_T")$tables$metadata_columns$column
  expect_false(identical(first, second))
  expect_equal(unname(other[[first, drop = TRUE]]), as.numeric(SeuratObject::LayerData(other, layer = "counts")["g3", ]))
})

test_that("AUCell applies the recorded seed and preserves caller RNG (#15)", {
  skip_if_not_installed("AUCell")
  set.seed(19)
  m <- matrix(rbinom(500 * 30, 1, .15), 500,
    dimnames = list(paste0("g", 1:500), paste0("c", 1:30)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(m, sparse = TRUE))
  signatures <- list(P = paste0("g", 1:20), Q = paste0("g", 21:40))
  before <- .Random.seed
  a <- sn_score_programs(object, signatures, method = "aucell", layer = "counts",
    seed = 777, return_object = FALSE)
  expect_identical(.Random.seed, before)
  runif(5)
  b <- sn_score_programs(object, signatures, method = "aucell", layer = "counts",
    seed = 777, return_object = FALSE)
  expect_identical(a$tables$scores, b$tables$scores)
  set.seed(777)
  rankings <- AUCell::AUCell_buildRankings(SeuratObject::LayerData(object, layer = "counts"),
    plotStats = FALSE, verbose = FALSE)
  expected <- AUCell::getAUC(AUCell::AUCell_calcAUC(signatures, rankings, aucMaxRank = 25, verbose = FALSE))
  expect_equal(a$tables$scores$score, as.vector(t(expected[names(signatures), , drop = FALSE])))
})

test_that("marker backgrounds precede significance filtering and preserve upstream output (#9)", {
  object <- audit_de_object()
  object$stratum <- NULL
  actual <- sn_find_de(object, analysis = "markers", group_by = "group",
    only_pos = FALSE, min_pct = .25, logfc_threshold = 0, return.thresh = 1e-8, verbose = FALSE)
  stored <- sn_get_result(actual, "de", "default")
  direct <- Seurat::FindAllMarkers(object, group.by = "group", only.pos = FALSE,
    min.pct = .25, logfc.threshold = 0, return.thresh = 1e-8, verbose = FALSE)
  expect_equal(stored$tables$primary, tibble::as_tibble(direct))
  for (group in c("A", "B")) {
    all_tests <- Seurat::FindMarkers(object, ident.1 = group, group.by = "group",
      only.pos = FALSE, min.pct = .25, logfc.threshold = 0, verbose = FALSE)
    background <- stored$input$tested_features_by_comparison
    expect_setequal(background$gene[background$cluster == group], rownames(all_tests)[is.finite(all_tests$p_val)])
  }
})

test_that("stored-DE ORA uses a separate background for each comparison (#9)", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  object <- audit_de_object()
  de <- list(analysis = "markers", method = "wilcox", assay = "RNA", group_col = "cluster",
    p_col = "p_val_adj", rank_col = "avg_log2FC",
    table = data.frame(gene = c("CD3D", "MS4A1"), cluster = c("A", "B"),
      avg_log2FC = c(2, 3), p_val_adj = c(.001, .002)),
    input = list(tested_features = c("CD3D", "CD3E", "TRAC", "MS4A1", "CD79A"),
      tested_feature_groups = "cluster",
      tested_features_by_comparison = data.frame(
        gene = c("CD3D", "CD3E", "TRAC", "MS4A1", "CD79A"), cluster = c("A", "A", "A", "B", "B"))))
  object <- sn_store_result(object, "de", "markers", de)
  observed <- list()
  local_mocked_bindings(enrichGO = function(gene, universe, ...) {
    observed[[gene[[1]]]] <<- universe
    data.frame(ID = "GO:test", Description = "test", pvalue = .01, p.adjust = .02)
  }, .package = "clusterProfiler")
  output <- sn_run_enrichment(object, source_de_result_id = "markers", analysis = "ora",
    species = "human", database = "GOBP", result_id = "go", return_object = TRUE)
  expect_setequal(observed$CD3D, c("CD3D", "CD3E", "TRAC"))
  expect_setequal(observed$MS4A1, c("MS4A1", "CD79A"))
  result <- sn_get_result(output, "enrichment", "go")
  expect_setequal(result$tables$primary$Cluster, c("A", "B"))
  expect_length(result$parameters$comparison_backgrounds, 2L)
})
