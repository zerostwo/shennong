library(testthat)

make_de_test_object <- function() {
  set.seed(717)
  genes <- c(
    paste0("GENE", seq_len(80)),
    "CD3D", "CD3E", "TRAC", "LCK", "LTB",
    "MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1", "LYZ"
  )
  counts <- matrix(rpois(length(genes) * 80, lambda = 2), nrow = length(genes), ncol = 80)
  rownames(counts) <- genes
  colnames(counts) <- paste0("cell", seq_len(80))

  sample <- rep(c("s1", "s2", "s3", "s4"), each = 20)
  condition <- rep(rep(c("control", "treated"), each = 10), 4)
  cell_type <- rep(rep(c("Tcell", "Bcell"), each = 5), 8)

  tcell_treated <- cell_type == "Tcell" & condition == "treated"
  bcell_control <- cell_type == "Bcell" & condition == "control"

  counts[rownames(counts) %in% c("CD3D", "CD3E", "TRAC", "LCK", "LTB"), tcell_treated] <-
    counts[rownames(counts) %in% c("CD3D", "CD3E", "TRAC", "LCK", "LTB"), tcell_treated] + 8
  counts[rownames(counts) %in% c("MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"), bcell_control] <-
    counts[rownames(counts) %in% c("MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"), bcell_control] + 8
  counts[rownames(counts) %in% "LYZ", cell_type == "Bcell"] <-
    counts[rownames(counts) %in% "LYZ", cell_type == "Bcell"] + 4

  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    project = "de-test"
  )
  Seurat::Misc(object = object, slot = "species") <- "human"
  object <- Seurat::AddMetaData(
    object = object,
    metadata = data.frame(
      sample = sample,
      condition = condition,
      cell_type = cell_type,
      row.names = colnames(object)
    )
  )
  object <- Seurat::SetIdent(object = object, value = "cell_type")
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object
}

test_that("sn_find_de stores marker results on the Seurat object", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    result_id = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  stored <- sn_get_result(object, "de", "celltype_markers")
  expect_true("gene" %in% colnames(stored$tables$primary))
  expect_equal(stored$schema_version, "2.0.0")
  expect_equal(stored$analysis, "markers")
  expect_equal(stored$method, "wilcox")
  expect_identical(stored$result_id, "celltype_markers")
  expect_setequal(stored$input$tested_features, rownames(object[["RNA"]]))
  expect_equal(stored$input$tested_features_count, nrow(object[["RNA"]]))
  expect_identical(stored$input$tested_features_source, "assay_layer")
})

test_that("stored-DE ORA passes the exact tested feature subset to enrichment backends", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("msigdbr")

  object <- make_de_test_object()
  tested_features <- c("CD3D", "CD3E", "TRAC", "LCK")
  local_mocked_bindings(
    .sn_run_seurat_de = function(...) {
      data.frame(
        avg_log2FC = c(1.5, 1),
        p_val_adj = c(0.001, 0.002),
        row.names = c("CD3D", "CD3E")
      )
    },
    .sn_enrich_get_msigdb_terms = function(...) {
      tibble::tibble(
        term = c("TEST_TERM", "TEST_TERM"),
        description = c("test term", "test term"),
        gene = c("CD3D", "CD3E")
      )
    },
    .package = "Shennong"
  )
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    features = tested_features,
    result_id = "feature_subset",
    return_object = TRUE,
    verbose = FALSE
  )
  stored_de <- sn_get_result(object, "de", "feature_subset")
  expect_identical(stored_de$input$tested_features, tested_features)
  expect_equal(stored_de$input$tested_features_count, length(tested_features))
  expect_identical(stored_de$input$tested_features_source, "requested_features")

  go_universe <- NULL
  msigdb_universe <- NULL
  mock_enrichment <- data.frame(
    ID = "TEST_TERM",
    Description = "test term",
    pvalue = 0.01,
    p.adjust = 0.02,
    qvalue = 0.02
  )
  enrichment_outputs <- with_mocked_bindings(
    {
      go_object <- sn_run_enrichment(
        x = object,
        source_de_result_id = "feature_subset",
        analysis = "ora",
        species = "human",
        database = "GOBP",
        result_id = "go_subset",
        return_object = TRUE
      )
      hallmark_result <- sn_run_enrichment(
        x = object,
        source_de_result_id = "feature_subset",
        analysis = "ora",
        species = "human",
        database = "H",
        result_id = "hallmark_subset",
        return_object = FALSE
      )
      list(go_object = go_object, hallmark_result = hallmark_result)
    },
    enrichGO = function(..., universe) {
      go_universe <<- universe
      mock_enrichment
    },
    enricher = function(..., universe) {
      msigdb_universe <<- universe
      mock_enrichment
    },
    .package = "clusterProfiler"
  )

  expect_identical(go_universe, tested_features)
  expect_identical(msigdb_universe, tested_features)
  stored_go <- sn_get_result(enrichment_outputs$go_object, "enrichment", "go_subset")
  expect_identical(stored_go$parameters$universe_source, "stored_de_tested_features")
  expect_equal(stored_go$parameters$universe_size, length(tested_features))
})

test_that("stored-DE universe resolution only falls back for legacy results", {
  skip_if_not_installed("Seurat")
  object <- make_de_test_object()

  legacy <- Shennong:::.sn_enrich_stored_de_universe(
    de_result = list(assay = "RNA", input = list()),
    object = object
  )
  expect_identical(legacy$source, "stored_de_assay_fallback:RNA")
  expect_setequal(legacy$universe, rownames(object[["RNA"]]))

  expect_error(
    Shennong:::.sn_enrich_stored_de_universe(
      de_result = list(
        assay = "RNA",
        input = list(tested_features = c("CD3D", NA_character_))
      ),
      object = object
    ),
    "tested_features"
  )
})

test_that("sn_find_de accepts an explicit result_id", {
  skip_if_not_installed("Seurat")
  object <- make_de_test_object()
  testthat::local_mocked_bindings(
    .sn_run_seurat_de = function(...) {
      data.frame(avg_log2FC = 1, p_val_adj = 0.01, row.names = "GENE1")
    },
    .package = "Shennong"
  )

  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    result_id = "donor1_celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  expect_true("donor1_celltype_markers" %in% names(object@misc$shennong$results$de))
  stored <- sn_get_result(object, "de", "donor1_celltype_markers")
  expect_identical(stored$result_id, "donor1_celltype_markers")
  listing <- sn_list_results(object, type = "de")
  expect_identical(listing$result_id, "donor1_celltype_markers")
  expect_error(
    sn_find_de(object, result_id = "", return_object = TRUE),
    "result_id"
  )
})

test_that("stored-DE ORA uses a signed effect rather than an unsigned ranking score", {
  input <- data.frame(
    gene = c("UP", "DOWN"),
    cosg_score = c(0.9, 0.8),
    avg_log2FC = c(1, -1),
    p_val_adj = c(0.01, 0.01)
  )
  de_result <- list(rank_col = "cosg_score", p_col = "p_val_adj")

  up <- Shennong:::.sn_enrich_filter_stored_de_ora(
    input, de_result, direction = "up"
  )
  down <- Shennong:::.sn_enrich_filter_stored_de_ora(
    input, de_result, direction = "down"
  )

  expect_identical(up$gene, "UP")
  expect_identical(down$gene, "DOWN")
  expect_identical(attr(up, "sn_de_ora_selection")$effect_column, "avg_log2FC")

  expect_error(
    Shennong:::.sn_enrich_filter_stored_de_ora(
      input[, c("gene", "cosg_score", "p_val_adj")],
      de_result,
      direction = "up"
    ),
    "signed effect column"
  )
})

test_that("sn_find_de preserves scoped Seurat acceleration provenance", {
  skip_if_not_installed("Seurat")
  object <- make_de_test_object()
  testthat::local_mocked_bindings(
    .sn_run_seurat_de = function(...) {
      Shennong:::.sn_record_acceleration_usage("seurat")
      data.frame(avg_log2FC = 1, p_val_adj = 0.01, row.names = "GENE1")
    },
    .sn_acceleration_provenance = function() {
      context <- getOption("shennong.acceleration.provenance_context")
      active <- if (is.environment(context)) context$used_patches else character()
      if (length(active) == 0L) list() else list(active_patches = active)
    },
    .package = "Shennong"
  )

  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    result_id = "accelerated_markers",
    return_object = TRUE,
    verbose = FALSE
  )
  stored <- sn_get_de_result(
    object,
    result_id = "accelerated_markers",
    with_metadata = TRUE
  )

  expect_identical(
    stored$provenance$acceleration$active_patches,
    "seurat"
  )
})

test_that("sn_plot_dot can use stored top markers", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    result_id = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  plot <- suppressWarnings(
    sn_plot_dot(
      x = object,
      features = "top_markers",
      result_id = "celltype_markers",
      n = 2
    )
  )

  expect_s3_class(plot, "ggplot")
})

test_that("sn_find_de supports contrasts within each subset", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  result <- sn_find_de(
    object,
    analysis = "contrast",
    ident_1 = "treated",
    ident_2 = "control",
    group_by = "condition",
    subset_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    return_object = FALSE,
    verbose = FALSE
  )

  expect_true(all(c("gene", "cell_type", "comparison") %in% colnames(result)))
  expect_setequal(unique(result$cell_type), c("Tcell", "Bcell"))
})

test_that("sn_find_de supports pseudobulk contrasts", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")

  object <- make_de_test_object()
  result <- sn_find_de(
    object,
    analysis = "pseudobulk",
    ident_1 = "treated",
    ident_2 = "control",
    group_by = "condition",
    subset_by = "cell_type",
    sample_by = "sample",
    method = "edgeR",
    min_cells_per_sample = 5,
    return_object = FALSE,
    verbose = FALSE
  )

  expect_true(nrow(result) > 0)
  expect_true(all(c("gene", "comparison", "cell_type", "log2FoldChange") %in% colnames(result)))
})

test_that("sn_find_de supports limma pseudobulk contrasts", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")
  skip_if_not_installed("limma")

  object <- make_de_test_object()
  result <- sn_find_de(
    object,
    analysis = "pseudobulk",
    ident_1 = "treated",
    ident_2 = "control",
    group_by = "condition",
    subset_by = "cell_type",
    sample_by = "sample",
    method = "limma",
    min_cells_per_sample = 5,
    return_object = FALSE,
    verbose = FALSE
  )

  expect_true(nrow(result) > 0)
  expect_true(all(c("gene", "comparison", "cell_type", "log2FoldChange") %in% colnames(result)))
})

test_that("pseudobulk rejects normalized layers and applies method-specific integer rules", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")
  object <- make_de_test_object()

  expect_error(
    sn_find_de(
      object,
      analysis = "pseudobulk",
      ident_1 = "treated",
      ident_2 = "control",
      group_by = "condition",
      sample_by = "sample",
      layer = "data",
      method = "edgeR",
      min_cells_per_sample = 5,
      return_object = FALSE,
      verbose = FALSE
    ),
    "requires a raw or corrected count layer"
  )

  fractional <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  fractional@x <- fractional@x + 0.25
  SeuratObject::LayerData(object, assay = "RNA", layer = "corrected_counts") <- fractional
  expect_error(
    Shennong:::.sn_validate_pseudobulk_count_layer(
      fractional,
      layer = "corrected_counts"
    ),
    "DESeq2.*integer-valued"
  )
  expect_no_error(Shennong:::.sn_validate_pseudobulk_count_layer(
    fractional,
    layer = "corrected_counts",
    require_integer = FALSE
  ))
  result <- sn_find_de(
    object,
    analysis = "pseudobulk",
    ident_1 = "treated",
    ident_2 = "control",
    group_by = "condition",
    sample_by = "sample",
    layer = "corrected_counts",
    method = "edgeR",
    min_cells_per_sample = 5,
    return_object = FALSE,
    verbose = FALSE
  )
  expect_gt(nrow(result), 0L)

  expect_no_error(Shennong:::.sn_validate_pseudobulk_count_layer(
    SeuratObject::LayerData(object, assay = "RNA", layer = "counts"),
    layer = "raw"
  ))
  expect_no_error(Shennong:::.sn_validate_pseudobulk_count_layer(
    SeuratObject::LayerData(object, assay = "RNA", layer = "counts"),
    layer = "umi"
  ))
  expect_error(
    Shennong:::.sn_validate_pseudobulk_count_layer(
      SeuratObject::LayerData(object, assay = "RNA", layer = "counts"),
      layer = "normalized_counts"
    ),
    "requires a raw or corrected count layer"
  )
})

test_that("subset_levels is explicit, unique, and observed", {
  skip_if_not_installed("Seurat")
  object <- make_de_test_object()

  expect_error(
    sn_find_de(
      object,
      analysis = "markers",
      group_by = "cell_type",
      subset_levels = "T",
      return_object = FALSE,
      verbose = FALSE
    ),
    "only be supplied together"
  )
  expect_error(
    sn_find_de(
      object,
      analysis = "markers",
      group_by = "condition",
      subset_by = "cell_type",
      subset_levels = c("T", "T"),
      return_object = FALSE,
      verbose = FALSE
    ),
    "distinct"
  )
  expect_error(
    sn_find_de(
      object,
      analysis = "markers",
      group_by = "condition",
      subset_by = "cell_type",
      subset_levels = "not_observed",
      return_object = FALSE,
      verbose = FALSE
    ),
    "not observed"
  )
})

test_that("pseudobulk profile identifiers cannot collide on user labels", {
  keys <- Shennong:::.sn_pseudobulk_profile_keys(
    sample = c("a___b", "a", "a___b", "a"),
    group = c("c", "b___c", "c", "b___c")
  )
  expect_identical(keys, c("profile_1", "profile_2", "profile_1", "profile_2"))
  expect_length(unique(keys), 2L)
})

test_that("pseudobulk supports explicit paired designs and validates sample semantics", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("edgeR")
  object <- make_de_test_object()

  result <- sn_find_de(
    object,
    analysis = "pseudobulk",
    ident_1 = "treated",
    ident_2 = "control",
    group_by = "condition",
    sample_by = "sample",
    method = "edgeR",
    min_cells_per_sample = 5,
    design = ~sample + condition,
    contrast = c("condition", "treated", "control"),
    return_object = FALSE,
    verbose = FALSE
  )
  expect_gt(nrow(result), 0L)
  expect_true(all(result$comparison == "treated vs control"))

  expect_error(
    sn_find_de(
      object,
      analysis = "pseudobulk",
      ident_1 = "treated",
      ident_2 = "control",
      group_by = "condition",
      sample_by = "sample",
      method = "edgeR",
      min_cells_per_sample = 5,
      design = ~condition,
      contrast = c("condition", "treated", "control"),
      return_object = FALSE,
      verbose = FALSE
    ),
    "must include `sample_by`"
  )

  mixed <- object
  mixed$sample[as.character(mixed$sample) == "s4" & mixed$condition == "treated"] <- "s4_treated"
  expect_error(
    sn_find_de(
      mixed,
      analysis = "pseudobulk",
      ident_1 = "treated",
      ident_2 = "control",
      group_by = "condition",
      sample_by = "sample",
      method = "edgeR",
      min_cells_per_sample = 5,
      return_object = FALSE,
      verbose = FALSE
    ),
    "mixture of paired and unpaired"
  )
})

test_that("sn_find_de supports COSGR markers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("COSG")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    method = "COSGR",
    n_genes_user = 10,
    result_id = "cosgr_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  stored <- sn_get_result(object, "de", "cosgr_markers")
  result <- stored$tables$primary
  expect_true(all(c("gene", "cluster", "cosg_score", "rank") %in% colnames(result)))
  expect_equal(stored$method, "COSGR")
})

test_that("sn_annotate_de_features annotates marker tables with custom resources", {
  marker_tbl <- tibble::tibble(
    gene = c("TBX21", "CXCL10", "IL7R", "ACTB"),
    cluster = c("Tcell", "Myeloid", "Tcell", "Bcell"),
    avg_log2FC = c(2.1, 1.8, 1.2, 0.4)
  )
  custom_resource <- tibble::tibble(
    gene = c("TBX21", "CXCL10", "IL7R", "IL7R"),
    feature_class = c("transcription_factor", "chemokine", "surface_membrane", "cytokine"),
    feature_class_source = "test"
  )

  result <- sn_annotate_de_features(
    marker_tbl,
    species = "human",
    resource = "custom",
    custom_resource = custom_resource
  )

  expect_true(all(c(
    "feature_classes",
    "is_transcription_factor",
    "is_surface_membrane",
    "is_cytokine",
    "is_chemokine"
  ) %in% colnames(result)))
  expect_true(result$is_transcription_factor[result$gene == "TBX21"])
  expect_true(result$is_chemokine[result$gene == "CXCL10"])
  expect_true(result$is_surface_membrane[result$gene == "IL7R"])
  expect_true(result$is_cytokine[result$gene == "IL7R"])
  expect_false(result$is_chemokine[result$gene == "ACTB"])
})

test_that("sn_annotate_de_features stores annotated DE results on Seurat objects", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    result_id = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )
  custom_resource <- tibble::tibble(
    gene = c("CD3D", "MS4A1"),
    feature_class = c("surface_membrane", "surface_membrane"),
    feature_class_source = "test"
  )

  object <- sn_annotate_de_features(
    object,
    source_result_id = "celltype_markers",
    resource = "custom",
    custom_resource = custom_resource,
    return_object = TRUE
  )

  expect_true("celltype_markers_feature_classes" %in% names(object@misc$shennong$results$de))
  stored <- sn_get_de_result(object, result_id = "celltype_markers_feature_classes")
  expect_true("is_surface_membrane" %in% colnames(stored))
  expect_true(any(stored$is_surface_membrane, na.rm = TRUE))
  stored_result <- sn_get_result(object, "de", "celltype_markers_feature_classes")
  expect_null(stored_result[["table"]])
  expect_true("is_surface_membrane" %in% colnames(stored_result$tables$primary))
  expect_equal(
    stored_result$feature_annotation$source_result_id,
    "celltype_markers"
  )
})

test_that("sn_run_enrichment supports GSEA from ranked marker tables", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  ranked_markers <- tibble::tibble(
    gene = c("CD3D", "CD3E", "TRAC", "LCK", "LTB", "IL7R", "MALAT1", "ACTB"),
    avg_log2FC = c(3.2, 3.1, 2.9, 2.5, 2.2, 1.8, -0.5, -1)
  )

  expect_no_error({
    result <- suppressWarnings(sn_run_enrichment(
      ranked_markers,
      gene_clusters = gene ~ avg_log2FC,
      analysis = "gsea",
      species = "human",
      database = "GOBP",
      min_gs_size = 2,
      pvalue_cutoff = 1
    ))
  })

  expect_true(inherits(result, "gseaResult"))
})

test_that("sn_run_enrichment stores enrichment results on the Seurat object by default", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    result_id = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  object <- suppressWarnings(sn_run_enrichment(
    x = object,
    source_de_result_id = "celltype_markers",
    gene_clusters = gene ~ cluster,
    species = "human",
    database = "GOBP",
    result_id = "demo_gsea"
  ))

  expect_s4_class(object, "Seurat")
  expect_true("demo_gsea" %in% names(object@misc$shennong$results$enrichment))
  stored <- sn_get_enrichment_result(
    object,
    result_id = "demo_gsea",
    with_metadata = TRUE
  )
  expect_identical(stored$parameters$p_adjust_method, "BH")
  expect_identical(stored$parameters$min_gs_size, 10L)
  expect_identical(stored$parameters$max_gs_size, 500L)
  expect_identical(
    stored$parameters$backend_versions[["clusterProfiler"]],
    as.character(utils::packageVersion("clusterProfiler"))
  )
  expect_true(ncol(stored$tables$primary) > 0L)
})

test_that("sn_run_enrichment validates object type and GSEA input contracts", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  expect_error(
    sn_run_enrichment(
      x = make_de_test_object(),
      species = "human",
      database = "GOBP",
      return_object = TRUE
    ),
    "No stored results"
  )

  expect_error(
    sn_run_enrichment(
      x = tibble::tibble(symbol = c("CD3D", "CD3E"), avg_logFC = c(1, 2)),
      analysis = "gsea",
      species = "human",
      database = "GOBP"
    ),
    "gene_clusters"
  )

  expect_error(
    sn_run_enrichment(
      x = tibble::tibble(gene = c("CD3D", "CD3E")),
      analysis = "gsea",
      species = "human",
      database = "GOBP"
    ),
    "gene_clusters"
  )

  expect_error(
    sn_run_enrichment(
      x = c("CD3D", "CD3E"),
      analysis = "gsea",
      species = "human",
      database = "GOBP"
    ),
    "named numeric vector or a data frame"
  )
})

test_that("sn_run_enrichment can write GSEA results to disk and compare grouped GO sets", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  ranked_df <- tibble::tibble(
    gene = c("CD3D", "CD3E", "TRAC", "LCK", "LTB", "IL7R", "MALAT1", "ACTB"),
    avg_logFC = c(3.2, 3.1, 2.9, 2.5, 2.2, 1.8, -0.5, -1)
  )
  outdir <- tempfile("enrich-out-")

  gsea_result <- suppressWarnings(sn_run_enrichment(
    x = ranked_df,
    gene_clusters = gene ~ avg_logFC,
    analysis = "gsea",
    species = "human",
    database = "GOBP",
    prefix = "demo",
    outdir = outdir,
    min_gs_size = 2,
    pvalue_cutoff = 1
  ))

  grouped_genes <- tibble::tibble(
    gene = c("CD3D", "CD3E", "TRAC", "LCK", "MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"),
    cluster = c("Tcell", "Tcell", "Tcell", "Tcell", "Bcell", "Bcell", "Bcell", "Bcell")
  )
  compare_result <- sn_run_enrichment(
    x = grouped_genes,
    gene_clusters = gene ~ cluster,
    analysis = "ora",
    species = "human",
    database = "GOBP",
    min_gs_size = 2,
    pvalue_cutoff = 1
  )

  expect_true(inherits(gsea_result, "gseaResult"))
  expect_true(file.exists(file.path(outdir, "demo.enrichment.GOBP.rds")))
  expect_true(inherits(compare_result, "compareClusterResult"))
})

test_that("sn_run_enrichment auto-detects categorical ORA and requires explicit numeric GSEA", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  ora_input <- tibble::tibble(
    gene = c("CD3D", "CD3E", "TRAC", "LCK", "MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"),
    cell_type = c("Tcell", "Tcell", "Tcell", "Tcell", "Bcell", "Bcell", "Bcell", "Bcell")
  )
  gsea_input <- tibble::tibble(
    gene = c("CD3D", "CD3E", "TRAC", "LCK", "LTB", "IL7R", "MALAT1", "ACTB"),
    log2fc = c(3.2, 3.1, 2.9, 2.5, 2.2, 1.8, -0.5, -1)
  )

  ora_result <- sn_run_enrichment(
    ora_input,
    gene_clusters = gene ~ cell_type,
    species = "human",
    database = "GOBP",
    min_gs_size = 2,
    pvalue_cutoff = 1
  )
  gsea_result <- suppressWarnings(sn_run_enrichment(
    gsea_input,
    gene_clusters = gene ~ log2fc,
    analysis = "gsea",
    species = "human",
    database = "GOBP",
    min_gs_size = 2,
    pvalue_cutoff = 1
  ))

  expect_true(inherits(ora_result, "compareClusterResult"))
  expect_true(inherits(gsea_result, "gseaResult"))
})

test_that("sn_run_enrichment supports multi-database requests and database-specific storage names", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("msigdbr")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    result_id = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )
  outdir <- tempfile("multi-enrich-")

  object <- suppressWarnings(sn_run_enrichment(
    x = object,
    source_de_result_id = "celltype_markers",
    gene_clusters = gene ~ cluster,
    species = "human",
    database = c("GOBP", "H"),
    result_id = "combined",
    prefix = "bundle",
    outdir = outdir,
    pvalue_cutoff = 1
  ))

  expect_s4_class(object, "Seurat")
  expect_true(all(c("combined.GOBP", "combined.H") %in% names(object@misc$shennong$results$enrichment)))
  expect_true(file.exists(file.path(outdir, "bundle.enrichment.GOBP.rds")))
  expect_true(file.exists(file.path(outdir, "bundle.enrichment.H.rds")))
})

test_that("sn_run_enrichment supports Seurat x input via stored DE results", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  object <- make_de_test_object()
  object <- sn_find_de(
    object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    result_id = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  object <- suppressWarnings(sn_run_enrichment(
    x = object,
    source_de_result_id = "celltype_markers",
    species = "human",
    database = "GOBP",
    result_id = "from_de",
    pvalue_cutoff = 1
  ))

  expect_s4_class(object, "Seurat")
  expect_true("from_de" %in% names(object@misc$shennong$results$enrichment))
  stored <- sn_get_result(object, "enrichment", "from_de")
  expect_identical(stored$parameters$de_ora_selection$group_column, "cluster")
  expect_identical(stored$parameters$universe_source, "stored_de_tested_features")
  expect_equal(stored$parameters$universe_size, nrow(object[["RNA"]]))
  expect_identical(stored$parameters$de_ora_selection$universe_source, "stored_de_tested_features")
  expect_equal(stored$parameters$de_ora_selection$universe_size, nrow(object[["RNA"]]))
})

test_that("sn_run_enrichment helper parsers validate formulas and msigdb inputs", {
  expect_null(Shennong:::.sn_enrich_parse_msigdb_database("GOBP"))
  expect_equal(
    Shennong:::.sn_enrich_parse_msigdb_database("C2:CP:REACTOME"),
    list(collection = "C2", subcollection = "CP:REACTOME")
  )
  expect_equal(
    Shennong:::.sn_enrich_parse_msigdb_database("MSIGDB", collection = "h", subcollection = "all"),
    list(collection = "H", subcollection = "ALL")
  )
  expect_error(
    Shennong:::.sn_enrich_parse_msigdb_database("MSIGDB"),
    "collection"
  )

  expect_equal(
    Shennong:::.sn_enrich_parse_formula(gene ~ cluster),
    list(gene_col = "gene", value_col = "cluster")
  )
  expect_error(
    Shennong:::.sn_enrich_parse_formula(~ cluster),
    "two-sided formula"
  )
  expect_error(
    Shennong:::.sn_enrich_parse_formula("gene ~ cluster"),
    "`gene_clusters` must be a two-sided formula"
  )
})

test_that("sn_run_enrichment helper resolution covers store names and analysis inference", {
  expect_equal(
    Shennong:::.sn_enrich_result_ids("default", c("GOBP", "H")),
    stats::setNames(c("default.GOBP", "default.H"), c("GOBP", "H"))
  )
  expect_error(
    Shennong:::.sn_enrich_result_ids(c("only", "two"), c("GOBP", "H", "KEGG")),
    "length 1 or match the length"
  )

  expect_equal(
    Shennong:::.sn_enrich_resolve_analysis(
      input = tibble::tibble(gene = c("A", "B"), group = c("x", "y")),
      mapping = list(gene_col = "gene", value_col = "group")
    ),
    "ora"
  )
  expect_error(
    Shennong:::.sn_enrich_resolve_analysis(
      input = tibble::tibble(gene = c("A", "B"), score = c(1, -1)),
      mapping = list(gene_col = "gene", value_col = "score")
    ),
    "must be supplied"
  )
  expect_equal(
    Shennong:::.sn_enrich_resolve_analysis(stats::setNames(c(1, -1), c("A", "B"))),
    "gsea"
  )
})

test_that("stored DE ORA selects significant genes with an explicit direction", {
  input <- data.frame(
    gene = c("A", "B", "C", "D"),
    avg_log2FC = c(2, -2, 3, 0),
    p_val_adj = c(0.01, 0.02, 0.2, 0.001)
  )
  metadata <- list(p_col = "p_val_adj", rank_col = "avg_log2FC")

  up <- Shennong:::.sn_enrich_filter_stored_de_ora(input, metadata, direction = "up")
  down <- Shennong:::.sn_enrich_filter_stored_de_ora(input, metadata, direction = "down")
  both <- Shennong:::.sn_enrich_filter_stored_de_ora(input, metadata, direction = "both")

  expect_identical(up$gene, "A")
  expect_identical(down$gene, "B")
  expect_setequal(both$gene, c("A", "B"))
  expect_equal(attr(up, "sn_de_ora_selection")$retained_rows, 1L)
  expect_error(
    Shennong:::.sn_enrich_filter_stored_de_ora(
      input[, c("gene", "avg_log2FC")], metadata, direction = "up"
    ),
    "adjusted-p-value"
  )
  expect_error(
    Shennong:::.sn_enrich_filter_stored_de_ora(
      input[, c("gene", "p_val_adj")], metadata, direction = "up"
    ),
    "signed effect"
  )
})

test_that("enrichment backends do not change caller RNG state", {
  set.seed(103)
  before <- .Random.seed
  Shennong:::.sn_enrichment_with_rng_preserved(runif(10))
  expect_identical(.Random.seed, before)
})

test_that("sn_run_enrichment helper utilities normalize labels and resolve inputs consistently", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  object <- sn_store_result(
    object, "de", "markers",
    list(
      analysis = "markers", method = "test", backend = "test",
      tables = list(primary = data.frame(
        gene = c("CD3D", "MS4A1"), cluster = c("T", "B")
      )),
      created_at = "2026-03-29 00:00:00 UTC"
    )
  )

  resolved <- Shennong:::.sn_enrich_resolve_input(object, source_de_result_id = "markers")
  resolved_default <- Shennong:::.sn_enrich_resolve_input(object)
  passthrough <- Shennong:::.sn_enrich_resolve_input(c("CD3D", "MS4A1"))

  expect_true(is.data.frame(resolved$input))
  expect_identical(resolved$object, object)
  expect_equal(resolved_default$source_de_result_id, "markers")
  expect_equal(passthrough$input, c("CD3D", "MS4A1"))
  expect_null(passthrough$object)
  expect_error(
    Shennong:::.sn_enrich_resolve_input(object, source_de_result_id = "missing"),
    "No result"
  )

  expect_equal(
    Shennong:::.sn_enrich_normalize_database_labels(c("gobp", "H", "gobp")),
    c("GOBP", "H")
  )
  expect_equal(
    Shennong:::.sn_enrich_output_label("C2:CP:REACTOME"),
    "C2_CP_REACTOME"
  )
})

test_that("sn_run_enrichment prefers stored default DE results on Seurat objects", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  object <- sn_store_result(
    object, "de", "default",
    list(
      analysis = "markers", method = "test", backend = "test",
      tables = list(primary = data.frame(
        gene = c("CD3D", "MS4A1"), cluster = c("T", "B")
      )),
      created_at = "2026-03-29 00:00:00 UTC"
    )
  )
  object <- sn_store_result(
    object, "de", "older",
    list(
      analysis = "markers", method = "test", backend = "test",
      tables = list(primary = data.frame(
        gene = c("LYZ", "S100A8"), cluster = c("M", "M")
      )),
      created_at = "2026-03-28 00:00:00 UTC"
    )
  )

  resolved <- Shennong:::.sn_enrich_resolve_input(object)

  expect_equal(resolved$source_de_result_id, "default")
  expect_equal(resolved$input$gene, c("CD3D", "MS4A1"))
})

test_that("sn_run_enrichment caches msigdbr term tables within a session", {
  skip_if_not_installed("msigdbr")
  calls <- 0L
  rm(list = ls(envir = Shennong:::.sn_enrichment_cache_env$msigdb_terms), envir = Shennong:::.sn_enrichment_cache_env$msigdb_terms)

  first <- with_mocked_bindings(
    Shennong:::.sn_enrich_get_msigdb_terms("human", "H"),
    msigdbr = function(...) {
      calls <<- calls + 1L
      data.frame(
        gs_name = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
        gs_description = "demo",
        gene_symbol = c("CD3D", "CD3E"),
        stringsAsFactors = FALSE
      )
    },
    .package = "msigdbr"
  )
  second <- with_mocked_bindings(
    Shennong:::.sn_enrich_get_msigdb_terms("human", "H"),
    msigdbr = function(...) {
      calls <<- calls + 1L
      stop("cache miss")
    },
    .package = "msigdbr"
  )

  expect_equal(calls, 1L)
  expect_equal(first, second)
})

test_that("sn_run_enrichment gene resolvers require an explicit GSEA duplicate policy", {
  ora_df <- data.frame(gene = c("CD3D", "CD3D", "MS4A1"), stringsAsFactors = FALSE)
  gsea_df <- data.frame(
    gene = c("CD3D", "CD3D", "MS4A1"),
    score = c(1.5, 2.5, -1),
    stringsAsFactors = FALSE
  )

  expect_equal(
    Shennong:::.sn_enrich_resolve_gene_vector(ora_df, gene_col = "gene"),
    c("CD3D", "MS4A1")
  )
  expect_error(
    Shennong:::.sn_enrich_resolve_gene_vector(data.frame(symbol = "CD3D"), gene_col = "gene"),
    "was not found"
  )

  expect_error(
    Shennong:::.sn_enrich_resolve_gene_list(
      gsea_df,
      mapping = list(gene_col = "gene", value_col = "score")
    ),
    "must be unique"
  )
  gene_list <- Shennong:::.sn_enrich_resolve_gene_list(
    gsea_df,
    mapping = list(gene_col = "gene", value_col = "score"),
    duplicate_gene_method = "max_abs"
  )
  expect_equal(unname(gene_list), c(2.5, -1))
  expect_equal(names(gene_list), c("CD3D", "MS4A1"))
  expect_error(
    Shennong:::.sn_enrich_resolve_gene_list(
      data.frame(gene = c("A", "B"), score = c("up", "down")),
      mapping = list(gene_col = "gene", value_col = "score")
    ),
    "must be numeric"
  )
  expect_error(
    Shennong:::.sn_enrich_resolve_gene_list(
      stats::setNames(c(1, Inf), c("A", "B"))
    ),
    "finite"
  )
})

test_that("sn_run_enrichment caches symbol-to-ENTREZ mappings within a session", {
  skip_if_not_installed("clusterProfiler")

  calls <- 0L
  rm(list = ls(envir = Shennong:::.sn_enrichment_cache_env$symbol_to_entrez), envir = Shennong:::.sn_enrichment_cache_env$symbol_to_entrez)

  first <- with_mocked_bindings(
    Shennong:::.sn_enrich_symbol_to_entrez(c("CD3E", "CD3D"), "org.Hs.eg.db"),
    bitr = function(geneID, ...) {
      calls <<- calls + 1L
      data.frame(
        SYMBOL = geneID,
        ENTREZID = seq_along(geneID),
        stringsAsFactors = FALSE
      )
    },
    .package = "clusterProfiler"
  )
  second <- with_mocked_bindings(
    Shennong:::.sn_enrich_symbol_to_entrez(c("CD3D", "CD3E"), "org.Hs.eg.db"),
    bitr = function(...) {
      calls <<- calls + 1L
      stop("cache miss")
    },
    .package = "clusterProfiler"
  )

  expect_equal(calls, 1L)
  expect_equal(first$SYMBOL, c("CD3D", "CD3E"))
  expect_equal(second$SYMBOL, c("CD3D", "CD3E"))
})

test_that("sn_run_enrichment muffles empty-result warnings but preserves other warnings", {
  expect_warning(
    Shennong:::.sn_enrich_muffle_empty_warning(
      warning("Different warning")
    ),
    "Different warning"
  )
  expect_warning(
    Shennong:::.sn_enrich_muffle_empty_warning(
      warning("Invalid p-values detected")
    ),
    "Invalid p-values"
  )
  expect_warning(
    Shennong:::.sn_enrich_muffle_empty_warning(
      warning("qvalue::qvalue() failed, returning NA for qvalue.")
    ),
    "qvalue"
  )

  expect_no_warning(
    Shennong:::.sn_enrich_muffle_empty_warning(
      warning("No enrichment found for the supplied genes")
    )
  )
})

test_that("sn_run_enrichment supports Hallmark ORA with grouped marker tables", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  rm(list = ls(envir = Shennong:::.sn_enrichment_cache_env$msigdb_terms), envir = Shennong:::.sn_enrichment_cache_env$msigdb_terms)

  grouped_genes <- tibble::tibble(
    gene = c(
      "CD3D", "CD3E", "TRAC", "LCK", "IL7R",
      "MKI67", "TOP2A", "UBE2C", "BIRC5", "HMGB2"
    ),
    cluster = c(
      "Tcell", "Tcell", "Tcell", "Tcell", "Tcell",
      "Cycling", "Cycling", "Cycling", "Cycling", "Cycling"
    )
  )

  expect_no_error({
    result <- sn_run_enrichment(
      grouped_genes,
      gene_clusters = gene ~ cluster,
      analysis = "ora",
      species = "human",
      database = "H",
      min_gs_size = 2,
      pvalue_cutoff = 1
    )
  })

  expect_true(inherits(result, "compareClusterResult"))
})

test_that("sn_run_enrichment supports msigdbr subcollections via database strings", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  rm(list = ls(envir = Shennong:::.sn_enrichment_cache_env$msigdb_terms), envir = Shennong:::.sn_enrichment_cache_env$msigdb_terms)

  genes <- c("CD3D", "CD3E", "TRAC", "LCK", "LAT", "ZAP70")

  expect_no_error({
    result <- sn_run_enrichment(
      genes,
      analysis = "ora",
      species = "human",
      database = "C2:CP:REACTOME",
      min_gs_size = 2,
      pvalue_cutoff = 1
    )
  })

  expect_true(inherits(result, "enrichResult"))
})
