.api_test_object <- function() {
  set.seed(91)
  counts <- matrix(rpois(40 * 40, 5), 40,
                   dimnames = list(paste0("gene", 1:40), paste0("cell", 1:40)))
  counts[1:5, 1:20] <- counts[1:5, 1:20] + 15
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$condition <- rep(c("A", "B"), each = 20)
  Seurat::NormalizeData(object, verbose = FALSE)
}

.api_test_results <- function() {
  object <- .api_test_object()
  de <- tibble::tibble(gene = paste0("gene", 1:6), cluster = rep(c("A", "B"), each = 3),
                      avg_log2FC = c(6, 5, -4, 3, 2, -1), p_val_adj = c(.01, .2, .01, .01, .2, .01))
  object <- sn_store_result(object, "de", "markers", list(
    table = de, analysis = "markers", rank_col = "avg_log2FC", group_col = "cluster", p_col = "p_val_adj"
  ))
  terms <- tibble::tibble(ID = paste0("term", 1:6), Description = paste0("term", 1:6),
                         Cluster = rep(c("A", "B"), each = 3), NES = 6:1,
                         p.adjust = c(.01, .2, .01, .01, .2, .01))
  sn_store_enrichment(object, terms, result_id = "markers.gsea.TEST", analysis = "gsea")
}

test_that("result selection applies direction independently of truncation", {
  object <- .api_test_results()
  up <- sn_get_de_result(object, "markers", direction = "up")
  down <- sn_get_de_result(object, "markers", direction = "down")
  expect_equal(up$avg_log2FC, c(6, 5, 3, 2))
  expect_equal(down$avg_log2FC, c(-4, -1))
  expect_equal(sn_get_de_result(object, "markers", direction = "up", top_n = 99)$gene, up$gene)
})

test_that("omitted result IDs resolve only unambiguous stored results", {
  object <- .api_test_results()
  expect_identical(sn_get_result(object, "de")$result_id, "markers")
  expect_equal(nrow(sn_get_enrichment_result(object)), 6L)
  expect_error(sn_get_result(object, "de", ""), "result_id")
  object <- sn_store_enrichment(object, data.frame(ID = "other", p.adjust = .1), result_id = "other")
  expect_error(sn_get_enrichment_result(object), "Multiple.*other.*markers|Multiple.*markers.*other")
  expect_equal(nrow(sn_get_enrichment_result(object, "other")), 1L)
})

test_that("DE and enrichment getters share explicit group or overall top-N scope", {
  object <- .api_test_results()
  expect_equal(nrow(sn_get_de_result(object, top_n = 2)), 4L)
  expect_equal(nrow(sn_get_enrichment_result(object, top_n = 2)), 4L)
  expect_equal(nrow(sn_get_de_result(object, top_n = 2, top_scope = "all")), 2L)
  expect_equal(nrow(sn_get_enrichment_result(object, top_n = 2, top_scope = "all")), 2L)
  expect_equal(nrow(sn_get_enrichment_result(object, top_n = 2, groups = "B")), 2L)
  expect_error(sn_get_de_result(object, top_n = -1), "top_n")
  expect_error(sn_get_enrichment_result(object, top_n = 1.5), "top_n")
  expect_error(sn_get_enrichment_result(object, with_metadata = TRUE), "unused argument")
})

test_that("significance and effect-size filtering belong to result retrieval", {
  object <- .api_test_results()
  result <- sn_get_de_result(object, p_adjusted_cutoff = .05, logfc_threshold = 2, direction = "up")
  expect_equal(result$gene, c("gene1", "gene4"))
  expect_equal(nrow(sn_get_enrichment_result(object, p_adjusted_cutoff = .05)), 4L)
  expect_equal(nrow(sn_get_result(object, "de")$tables$primary), 6L)
  expect_false(any(c("p_val_cutoff", "de_logfc") %in% names(formals(sn_find_de))))
})

test_that("enrichment selection diagnoses missing columns and ranks numeric ratios", {
  object <- .api_test_object()
  object <- sn_store_enrichment(object, data.frame(ID = c("t1", "t2"),
    GeneRatio = c("1/2", "2/3")), result_id = "ratio")
  expect_identical(sn_get_enrichment_result(object, top_n = 1)$ID, "t2")
  expect_error(sn_get_enrichment_result(object, groups = "A"), "group.*column|grouping")
  object <- sn_store_enrichment(object, data.frame(ID = c("t1", "t2")), result_id = "no_rank")
  expect_error(sn_get_enrichment_result(object, "no_rank", top_n = 1), "rank|ranking")
})

test_that("scoring preserves separate default runs and rejects accidental replacement", {
  object <- .api_test_object()
  first <- sn_score_programs(object, list(first = c("gene1", "gene2")), method = "mean")
  second <- sn_score_programs(first, list(second = c("gene3", "gene4")), method = "mean")
  ids <- sn_list_results(second, type = "program_scoring")$result_id
  expect_length(ids, 2L)
  expect_identical(sn_get_result(second, "program_scoring", "programs_mean")$tables$coverage$program, "first")
  expect_error(sn_score_programs(first, list(second = "gene3"), method = "mean", result_id = "programs_mean"),
               "already exists")
  replaced <- sn_score_programs(first, list(second = "gene3"), method = "mean",
                               result_id = "programs_mean", overwrite = TRUE)
  expect_false("programs_mean_first" %in% colnames(replaced[[]]))
  expect_true("programs_mean_second" %in% colnames(replaced[[]]))
  expect_length(sn_list_results(replaced)$result_id, 1L)
})

test_that("scoring aggregation explicitly distinguishes expression and scores", {
  skip_if_not_installed("UCell")
  object <- .api_test_object()
  signatures <- list(program = paste0("gene", 1:5))
  cell <- sn_score_programs(object, signatures, method = "ucell", return_object = FALSE)
  scores <- sn_score_programs(object, signatures, method = "ucell", group_by = "condition",
                             aggregate = "scores", return_object = FALSE)
  expression <- sn_score_programs(object, signatures, method = "ucell", group_by = "condition",
                                 aggregate = "expression", return_object = FALSE)
  expected <- tapply(cell$tables$primary$score, object$condition, mean)
  expect_equal(scores$tables$primary$score, as.numeric(expected))
  expect_false(isTRUE(all.equal(scores$tables$primary$score, expression$tables$primary$score)))
  expect_error(sn_score_programs(object, signatures, group_by = "condition"), "aggregate")
  expect_error(sn_score_programs(object, signatures, aggregate = "scores"), "group_by")
  expect_identical(scores$parameters$aggregate, "scores")
  expect_identical(expression$parameters$aggregate, "expression")
})

test_that("public clustering inputs and common controls are discoverable", {
  expected <- c("batch_by", "backend_control", "assay", "layer", "dims", "npcs", "hvg_group_by")
  expect_true(all(expected %in% names(formals(sn_run_cluster))))
  expect_false(any(c("batch", "integration_control", "cluster_random_seed") %in% names(formals(sn_run_cluster))))
  expect_true("seed" %in% names(formals(sn_score_programs)))
  for (fun in list(sn_find_doublets, sn_run_bulk_deconvolution, sn_compare_integrations, sn_run_infercnvpy)) {
    expect_true("n_workers" %in% names(formals(fun)))
  }
})

test_that("DE returns a unified result and tables are retrieved explicitly", {
  result <- sn_find_de(.api_test_object(), ident_1 = "A", ident_2 = "B", group_by = "condition",
                      min_pct = 0, logfc_threshold = 0, return_object = FALSE, verbose = FALSE)
  expect_true(sn_validate_result(result)$valid)
  expect_equal(nrow(result$tables$primary), 40L)
  expect_equal(nrow(sn_get_de_result(result, p_adjusted_cutoff = 0)), 0L)
})

test_that("integration shortcuts share the same clustering workflow and controls", {
  local_mocked_bindings(.sn_run_cluster_impl = function(args) args, .package = "Shennong")
  for (method in c("scvi", "scanvi", "scpoli")) {
    result <- get(paste0("sn_run_", method))(NULL, batch_by = "sample",
      backend_control = list(n_epochs = 3), assay = "RNA", layer = "counts", npcs = 12, dims = 1:12, seed = 29)
    expect_identical(result$integration_method, method)
    expect_identical(result$batch, "sample")
    expect_equal(result$backend_control$seed, 29)
    expect_equal(result$backend_control$n_epochs, 3)
    expect_equal(result$npcs, 12)
    expect_identical(result$dims, 1:12)
  }
})

test_that("Milo honors the return flag and seeds actual neighborhood sampling", {
  skip_if_not_installed("miloR")
  set.seed(4)
  counts <- matrix(rpois(60 * 60, 3), 60,
    dimnames = list(paste0("gene", 1:60), paste0("cell", 1:60)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$sample <- rep(paste0("S", 1:6), each = 10)
  object$condition <- rep(c("A", "B"), each = 30)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  embedding <- matrix(rnorm(60 * 5), 60, dimnames = list(colnames(object), paste0("PC_", 1:5)))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = embedding, key = "PC_", assay = "RNA")
  run <- function(seed, return_object = FALSE) suppressWarnings(sn_run_milo(
    object, sample_by = "sample", group_by = "condition", dims = 1:5, k = 10,
    prop = .5, refined = FALSE, seed = seed, return_object = return_object, keep_model = TRUE, verbose = FALSE))
  before <- .Random.seed
  a <- run(717)
  expect_identical(.Random.seed, before)
  runif(3)
  b <- run(717)
  c <- run(999)
  expect_identical(miloR::nhoodIndex(a$models$milo), miloR::nhoodIndex(b$models$milo))
  expect_false(identical(miloR::nhoodIndex(a$models$milo), miloR::nhoodIndex(c$models$milo)))
  expect_equal(a$provenance$random_seed, 717)
  expect_s4_class(run(717, TRUE), "Seurat")

  sce <- Seurat::as.SingleCellExperiment(object)
  SingleCellExperiment::reducedDim(sce, "shennong_milo") <- embedding
  direct <- miloR::buildGraph(miloR::Milo(sce), k = 10, d = 5, reduced.dim = "shennong_milo")
  set.seed(717)
  direct <- miloR::makeNhoods(direct, prop = .5, k = 10, d = 5, refined = FALSE, reduced_dims = "shennong_milo")
  expect_identical(miloR::nhoodIndex(a$models$milo), miloR::nhoodIndex(direct))
})


test_that("result storage requires explicit replacement", {
  object <- .api_test_results()
  result <- sn_get_result(object, "de")
  expect_error(sn_store_result(object, "de", "markers", result), "already exists")
  expect_s4_class(sn_store_result(object, "de", "markers", result, overwrite = TRUE), "Seurat")
  expect_error(sn_find_de(object, p_val_cutoff = .01), "Removed DE argument")
})

test_that("enrichment writer and reader roundtrip with one or multiple databases", {
  skip_if_not_installed("clusterProfiler")
  terms <- data.frame(term = rep(c("T1", "T2"), each = 5),
                      gene = paste0("gene", 1:10), description = rep(c("one", "two"), each = 5))
  local_mocked_bindings(.sn_enrich_get_msigdb_terms = function(...) terms, .package = "Shennong")
  object <- .api_test_results()
  run <- function(x, database, ...) sn_run_enrichment(x, analysis = "ora", species = "human",
    database = database, min_gs_size = 1, pvalue_cutoff = 1, qvalue_cutoff = 1, ...)
  stored <- run(object, "H")
  id <- setdiff(sn_list_results(stored, type = "enrichment")$result_id,
                sn_list_results(object, type = "enrichment")$result_id)
  expect_length(id, 1L)
  standalone <- run(paste0("gene", 1:3), "H")
  multiple <- run(paste0("gene", 1:3), c("H", "C2"))
  expect_true(sn_validate_result(standalone)$valid)
  expect_true(sn_validate_result(multiple)$valid)
  expect_setequal(unique(multiple$tables$primary$database), c("H", "C2"))
  expect_named(multiple$models$database_results, c("H", "C2"))
  fresh <- sn_find_de(.api_test_object(), group_by = "condition", min_pct = 0,
                      logfc_threshold = 0, verbose = FALSE)
  fresh <- run(fresh, "H", de_p_adjusted_cutoff = 1)
  expect_equal(sn_get_enrichment_result(fresh), sn_get_result(fresh, "enrichment")$tables$primary)
})
