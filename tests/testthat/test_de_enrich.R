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
    store_name = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  expect_true("de_results" %in% names(methods::slot(object, "misc")))
  expect_true("celltype_markers" %in% names(object@misc$de_results))
  expect_true("gene" %in% colnames(object@misc$de_results$celltype_markers$table))
  expect_equal(object@misc$de_results$celltype_markers$schema_version, "1.0.0")
  expect_equal(object@misc$de_results$celltype_markers$analysis, "markers")
  expect_equal(object@misc$de_results$celltype_markers$method, "wilcox")
  expect_true("package_version" %in% names(object@misc$de_results$celltype_markers))
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
    store_name = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  plot <- suppressWarnings(
    sn_plot_dot(
      x = object,
      features = "top_markers",
      de_name = "celltype_markers",
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
    sample_col = "sample",
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
    sample_col = "sample",
    method = "limma",
    min_cells_per_sample = 5,
    return_object = FALSE,
    verbose = FALSE
  )

  expect_true(nrow(result) > 0)
  expect_true(all(c("gene", "comparison", "cell_type", "log2FoldChange") %in% colnames(result)))
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
    store_name = "cosgr_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  result <- object@misc$de_results$cosgr_markers$table
  expect_true(all(c("gene", "cluster", "cosg_score", "rank") %in% colnames(result)))
  expect_equal(object@misc$de_results$cosgr_markers$method, "COSGR")
})

test_that("sn_enrich supports GSEA from ranked marker tables", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  ranked_markers <- tibble::tibble(
    gene = c("CD3D", "CD3E", "TRAC", "LCK", "LTB", "IL7R", "MALAT1", "ACTB"),
    avg_log2FC = c(3.2, 3.1, 2.9, 2.5, 2.2, 1.8, -0.5, -1)
  )

  expect_no_error({
    result <- sn_enrich(
      ranked_markers,
      gene_clusters = gene ~ avg_log2FC,
      species = "human",
      database = "GOBP",
      pvalue_cutoff = 1
    )
  })

  expect_true(inherits(result, "gseaResult"))
})

test_that("sn_enrich stores enrichment results on the Seurat object by default", {
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
    store_name = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  object <- sn_enrich(
    x = object,
    source_de_name = "celltype_markers",
    gene_clusters = gene ~ cluster,
    species = "human",
    database = "GOBP",
    store_name = "demo_gsea"
  )

  expect_s4_class(object, "Seurat")
  expect_true("demo_gsea" %in% names(object@misc$enrichment_results))
})

test_that("sn_enrich validates object type and GSEA input contracts", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  expect_error(
    sn_enrich(
      x = make_de_test_object(),
      species = "human",
      database = "GOBP",
      return_object = TRUE
    ),
    "no stored DE results"
  )

  expect_error(
    sn_enrich(
      x = tibble::tibble(symbol = c("CD3D", "CD3E"), avg_logFC = c(1, 2)),
      analysis = "gsea",
      species = "human",
      database = "GOBP"
    ),
    "gene_clusters"
  )

  expect_error(
    sn_enrich(
      x = tibble::tibble(gene = c("CD3D", "CD3E")),
      analysis = "gsea",
      species = "human",
      database = "GOBP"
    ),
    "gene_clusters"
  )

  expect_error(
    sn_enrich(
      x = c("CD3D", "CD3E"),
      analysis = "gsea",
      species = "human",
      database = "GOBP"
    ),
    "named numeric vector or a data frame"
  )
})

test_that("sn_enrich can write GSEA results to disk and compare grouped GO sets", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  ranked_df <- tibble::tibble(
    gene = c("CD3D", "CD3E", "TRAC", "LCK", "LTB", "IL7R", "MALAT1", "ACTB"),
    avg_logFC = c(3.2, 3.1, 2.9, 2.5, 2.2, 1.8, -0.5, -1)
  )
  outdir <- tempfile("enrich-out-")

  gsea_result <- sn_enrich(
    x = ranked_df,
    gene_clusters = gene ~ avg_logFC,
    species = "human",
    database = "GOBP",
    prefix = "demo",
    outdir = outdir,
    pvalue_cutoff = 1
  )

  grouped_genes <- tibble::tibble(
    gene = c("CD3D", "CD3E", "TRAC", "LCK", "MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"),
    cluster = c("Tcell", "Tcell", "Tcell", "Tcell", "Bcell", "Bcell", "Bcell", "Bcell")
  )
  compare_result <- sn_enrich(
    x = grouped_genes,
    gene_clusters = gene ~ cluster,
    analysis = "ora",
    species = "human",
    database = "GOBP",
    pvalue_cutoff = 1
  )

  expect_true(inherits(gsea_result, "gseaResult"))
  expect_true(file.exists(file.path(outdir, "demo.enrichment.GOBP.rds")))
  expect_true(inherits(compare_result, "compareClusterResult"))
})

test_that("sn_enrich auto-detects ORA and GSEA from formula inputs", {
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

  ora_result <- sn_enrich(
    ora_input,
    gene_clusters = gene ~ cell_type,
    species = "human",
    database = "GOBP",
    pvalue_cutoff = 1
  )
  gsea_result <- sn_enrich(
    gsea_input,
    gene_clusters = gene ~ log2fc,
    species = "human",
    database = "GOBP",
    pvalue_cutoff = 1
  )

  expect_true(inherits(ora_result, "compareClusterResult"))
  expect_true(inherits(gsea_result, "gseaResult"))
})

test_that("sn_enrich supports multi-database requests and database-specific storage names", {
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
    store_name = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )
  outdir <- tempfile("multi-enrich-")

  object <- sn_enrich(
    x = object,
    source_de_name = "celltype_markers",
    gene_clusters = gene ~ cluster,
    species = "human",
    database = c("GOBP", "H"),
    store_name = "combined",
    prefix = "bundle",
    outdir = outdir,
    pvalue_cutoff = 1
  )

  expect_s4_class(object, "Seurat")
  expect_true(all(c("combined.GOBP", "combined.H") %in% names(object@misc$enrichment_results)))
  expect_true(file.exists(file.path(outdir, "bundle.enrichment.GOBP.rds")))
  expect_true(file.exists(file.path(outdir, "bundle.enrichment.H.rds")))
})

test_that("sn_enrich supports Seurat x input via stored DE results", {
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
    store_name = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )

  object <- sn_enrich(
    x = object,
    source_de_name = "celltype_markers",
    gene_clusters = gene ~ cluster,
    species = "human",
    database = "GOBP",
    store_name = "from_de",
    pvalue_cutoff = 1
  )

  expect_s4_class(object, "Seurat")
  expect_true("from_de" %in% names(object@misc$enrichment_results))
})

test_that("sn_enrich uses raw p-value filtering for enrichResult outputs", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  result <- methods::new(
    "enrichResult",
    result = data.frame(
      ID = c("A", "B"),
      Description = c("term A", "term B"),
      pvalue = c(0.01, 0.2),
      p.adjust = c(0.8, 0.01),
      qvalue = c(0.8, 0.01),
      geneID = c("CD3D/CD3E", "MS4A1/CD79A"),
      Count = c(2L, 2L)
    ),
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    organism = "human",
    ontology = "BP",
    gene = c("CD3D", "CD3E", "MS4A1", "CD79A"),
    keytype = "SYMBOL",
    universe = character(),
    gene2Symbol = character(),
    geneSets = list(A = c("CD3D", "CD3E"), B = c("MS4A1", "CD79A")),
    readable = FALSE,
    termsim = matrix(0, nrow = 0, ncol = 0),
    method = character(),
    dr = list()
  )

  filtered <- Shennong:::.sn_enrich_filter_by_pvalue(result, pvalue_cutoff = 0.05)

  expect_equal(nrow(filtered@result), 1)
  expect_equal(filtered@result$ID[[1]], "A")
})

test_that("sn_enrich helper parsers validate formulas and msigdb inputs", {
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

test_that("sn_enrich helper resolution covers store names and analysis inference", {
  expect_equal(
    Shennong:::.sn_enrich_store_names("default", c("GOBP", "H")),
    stats::setNames(c("default.GOBP", "default.H"), c("GOBP", "H"))
  )
  expect_error(
    Shennong:::.sn_enrich_store_names(c("only", "two"), c("GOBP", "H", "KEGG")),
    "length 1 or match the length"
  )

  expect_equal(
    Shennong:::.sn_enrich_resolve_analysis(
      input = tibble::tibble(gene = c("A", "B"), group = c("x", "y")),
      mapping = list(gene_col = "gene", value_col = "group")
    ),
    "ora"
  )
  expect_equal(
    Shennong:::.sn_enrich_resolve_analysis(
      input = tibble::tibble(gene = c("A", "B"), score = c(1, -1)),
      mapping = list(gene_col = "gene", value_col = "score")
    ),
    "gsea"
  )
  expect_equal(
    Shennong:::.sn_enrich_resolve_analysis(stats::setNames(c(1, -1), c("A", "B"))),
    "gsea"
  )
})

test_that("sn_enrich helper utilities normalize labels and resolve inputs consistently", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  object@misc$de_results <- list(
    markers = list(
      table = data.frame(gene = c("CD3D", "MS4A1"), cluster = c("T", "B")),
      analysis = "markers",
      created_at = "2026-03-29 00:00:00 UTC"
    )
  )

  resolved <- Shennong:::.sn_enrich_resolve_input(object, source_de_name = "markers")
  resolved_default <- Shennong:::.sn_enrich_resolve_input(object)
  passthrough <- Shennong:::.sn_enrich_resolve_input(c("CD3D", "MS4A1"))

  expect_true(is.data.frame(resolved$input))
  expect_identical(resolved$object, object)
  expect_equal(resolved_default$source_de_name, "markers")
  expect_equal(passthrough$input, c("CD3D", "MS4A1"))
  expect_null(passthrough$object)
  expect_error(
    Shennong:::.sn_enrich_resolve_input(object, source_de_name = "missing"),
    "was not found"
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

test_that("sn_enrich prefers stored default DE results on Seurat objects", {
  skip_if_not_installed("Seurat")

  object <- make_de_test_object()
  object@misc$de_results <- list(
    default = list(
      table = data.frame(gene = c("CD3D", "MS4A1"), cluster = c("T", "B")),
      analysis = "markers",
      created_at = "2026-03-29 00:00:00 UTC"
    ),
    older = list(
      table = data.frame(gene = c("LYZ", "S100A8"), cluster = c("M", "M")),
      analysis = "markers",
      created_at = "2026-03-28 00:00:00 UTC"
    )
  )

  resolved <- Shennong:::.sn_enrich_resolve_input(object)

  expect_equal(resolved$source_de_name, "default")
  expect_equal(resolved$input$gene, c("CD3D", "MS4A1"))
})

test_that("sn_enrich caches msigdbr term tables within a session", {
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

test_that("sn_enrich gene resolvers deduplicate ORA and GSEA inputs meaningfully", {
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

  gene_list <- Shennong:::.sn_enrich_resolve_gene_list(
    gsea_df,
    mapping = list(gene_col = "gene", value_col = "score")
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
})

test_that("sn_enrich caches symbol-to-ENTREZ mappings within a session", {
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

test_that("sn_enrich p-value filtering handles compareCluster and passthrough objects", {
  compare_result <- methods::new(
    "compareClusterResult",
    compareClusterResult = data.frame(
      Cluster = c("A", "A", "B"),
      ID = c("term1", "term2", "term3"),
      pvalue = c(0.01, 0.2, 0.03),
      stringsAsFactors = FALSE
    ),
    geneClusters = list(A = c("CD3D"), B = c("MS4A1")),
    fun = "enrichGO",
    gene2Symbol = character(),
    keytype = "SYMBOL",
    readable = FALSE,
    .call = call("compareCluster"),
    termsim = matrix(0, nrow = 0, ncol = 0),
    method = character(),
    dr = list()
  )

  filtered <- Shennong:::.sn_enrich_filter_by_pvalue(compare_result, pvalue_cutoff = 0.05)
  expect_equal(nrow(filtered@compareClusterResult), 2)
  expect_equal(filtered@compareClusterResult$ID, c("term1", "term3"))

  passthrough <- Shennong:::.sn_enrich_filter_by_pvalue(list(a = 1), pvalue_cutoff = 0.05)
  expect_equal(passthrough, list(a = 1))
})

test_that("sn_enrich supports Hallmark ORA with grouped marker tables", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")

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
    result <- sn_enrich(
      grouped_genes,
      gene_clusters = gene ~ cluster,
      analysis = "ora",
      species = "human",
      database = "H",
      pvalue_cutoff = 1
    )
  })

  expect_true(inherits(result, "compareClusterResult"))
})

test_that("sn_enrich supports msigdbr subcollections via database strings", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")

  genes <- c("CD3D", "CD3E", "TRAC", "LCK", "LAT", "ZAP70")

  expect_no_error({
    result <- sn_enrich(
      genes,
      analysis = "ora",
      species = "human",
      database = "C2:CP:REACTOME",
      pvalue_cutoff = 1
    )
  })

  expect_true(inherits(result, "enrichResult"))
})
