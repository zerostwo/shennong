.make_conformance_enrichment_terms <- function() {
  genes <- paste0("G", seq_len(40L))
  memberships <- list(
    UP = genes[1:10],
    DOWN = genes[31:40],
    MIXED = genes[c(1:3, 18:22, 38:40)],
    BACKGROUND_A = genes[8:19],
    BACKGROUND_B = genes[20:31],
    OUTSIDE_HEAVY = c("G1", "G2", paste0("X", 1:8))
  )
  dplyr::bind_rows(lapply(names(memberships), function(term) {
    data.frame(
      term = term,
      description = paste(term, "pathway"),
      gene = memberships[[term]],
      stringsAsFactors = FALSE
    )
  }))
}

.canonical_enrichment_table <- function(result) {
  table <- as.data.frame(result)
  if (nrow(table) == 0L) {
    return(table)
  }
  table <- table[order(table$ID), , drop = FALSE]
  rownames(table) <- NULL
  table
}

test_that("sn_enrich ORA matches clusterProfiler::enricher including its universe and cutoffs", {
  .conformance_require_package("clusterProfiler")
  terms <- .make_conformance_enrichment_terms()
  term2gene <- unique(terms[, c("term", "gene")])
  term2name <- unique(terms[, c("term", "description")])
  genes <- paste0("G", 1:8)
  universe <- paste0("G", 1:35)
  before <- .conformance_fingerprint(list(genes = genes, universe = universe))

  upstream <- .conformance_without_acceleration(clusterProfiler::enricher(
    gene = genes,
    universe = universe,
    pvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 3L,
    maxGSSize = 20L,
    TERM2GENE = term2gene,
    TERM2NAME = term2name
  ))

  local_mocked_bindings(
    .sn_enrich_get_msigdb_terms = function(...) terms,
    .package = "Shennong"
  )
  candidate <- .conformance_without_acceleration(sn_enrich(
    x = genes,
    analysis = "ora",
    species = "human",
    database = "H",
    universe = universe,
    pvalue_cutoff = 0.5,
    p_adjust_method = "BH",
    qvalue_cutoff = 1,
    min_gs_size = 3L,
    max_gs_size = 20L
  ))

  expect_equal(
    .canonical_enrichment_table(candidate),
    .canonical_enrichment_table(upstream),
    tolerance = 1e-12
  )
  expect_false("OUTSIDE_HEAVY" %in% as.data.frame(candidate)$ID)
  .conformance_expect_unchanged(
    list(genes = genes, universe = universe),
    before,
    "ORA inputs"
  )
})

test_that("sn_enrich GSEA matches clusterProfiler::GSEA with a controlled seed", {
  .conformance_require_package("clusterProfiler")
  terms <- .make_conformance_enrichment_terms()
  term2gene <- unique(terms[, c("term", "gene")])
  term2name <- unique(terms[, c("term", "description")])
  scores <- stats::setNames(
    c(seq(4, 0.1, length.out = 20L), seq(-0.1, -4, length.out = 20L)),
    paste0("G", seq_len(40L))
  )
  scores <- sort(scores, decreasing = TRUE)
  before <- .conformance_fingerprint(scores)

  set.seed(717L)
  upstream <- .conformance_without_acceleration(clusterProfiler::GSEA(
    geneList = scores,
    exponent = 1,
    minGSSize = 3L,
    maxGSSize = 20L,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    verbose = FALSE
  ))

  local_mocked_bindings(
    .sn_enrich_get_msigdb_terms = function(...) terms,
    .package = "Shennong"
  )
  set.seed(717L)
  candidate <- .conformance_without_acceleration(sn_enrich(
    x = scores,
    analysis = "gsea",
    species = "human",
    database = "H",
    pvalue_cutoff = 1,
    p_adjust_method = "BH",
    min_gs_size = 3L,
    max_gs_size = 20L,
    gsea_exponent = 1,
    duplicate_gene_method = "error"
  ))

  expect_equal(
    .canonical_enrichment_table(candidate),
    .canonical_enrichment_table(upstream),
    tolerance = 1e-12
  )
  expect_false("OUTSIDE_HEAVY" %in% as.data.frame(candidate)$ID)
  .conformance_expect_unchanged(scores, before, "GSEA ranked list")
})

test_that("sn_enrich retains parameters without false MSigDB patch usage", {
  .conformance_require_package("clusterProfiler")
  .conformance_require_package("SeuratObject")
  terms <- .make_conformance_enrichment_terms()
  counts <- matrix(
    1,
    nrow = 8L,
    ncol = 4L,
    dimnames = list(paste0("G", 1:8), paste0("stored_cell", 1:4))
  )
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object@misc$de_results <- list(markers = list(
    table = data.frame(gene = paste0("G", 1:8)),
    analysis = "markers",
    created_at = "2026-08-21 00:00:00 UTC"
  ))

  local_mocked_bindings(
    .sn_enrich_get_msigdb_terms = function(...) terms,
    .sn_with_default_autozyme = function(expr, patches, ...) {
      Shennong:::.sn_record_autozyme_usage(patches)
      force(expr)
    },
    .sn_autozyme_provenance = function() {
      context <- getOption("shennong.autozyme.provenance_context")
      used <- if (is.environment(context)) context$used_patches else character()
      if (length(used) == 0L) list() else list(active_patches = used)
    },
    .package = "Shennong"
  )
  object <- sn_enrich(
    x = object,
    source_de_name = "markers",
    analysis = "ora",
    species = "human",
    database = "H",
    universe = paste0("G", 1:35),
    min_gs_size = 3L,
    max_gs_size = 20L,
    pvalue_cutoff = 1,
    qvalue_cutoff = 1,
    store_name = "conformance_ora",
    return_object = TRUE
  )
  stored <- sn_get_enrichment_result(
    object,
    enrichment_name = "conformance_ora",
    with_metadata = TRUE
  )

  expect_null(stored$provenance$acceleration)
  expect_identical(stored$parameters$universe, paste0("G", 1:35))
  expect_identical(stored$parameters$min_gs_size, 3L)
  expect_identical(stored$parameters$backend_versions[["enrichit"]], "0.2.1")
})
