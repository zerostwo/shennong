library(testthat)

make_interpretation_base_object <- function() {
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

  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    project = "interpretation-test"
  )
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

make_interpretation_object <- function() {
  object <- make_interpretation_base_object()

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

  enrich_tbl <- tibble::tibble(
    Cluster = c("Tcell", "Bcell"),
    ID = c("GO:00001", "GO:00002"),
    Description = c("T cell activation", "B cell receptor signaling"),
    NES = c(2.1, 1.8),
    p.adjust = c(0.01, 0.03)
  )

  object <- sn_store_enrichment(
    object = object,
    result = enrich_tbl,
    store_name = "celltype_gsea",
    analysis = "gsea",
    database = "GOBP",
    species = "human",
    source_de_name = "celltype_markers",
    return_object = TRUE
  )

  object
}

test_that("sn_store_enrichment stores enrichment results on the Seurat object", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()

  expect_true("enrichment_results" %in% names(methods::slot(object, "misc")))
  expect_true("celltype_gsea" %in% names(object@misc$enrichment_results))
  expect_equal(object@misc$enrichment_results$celltype_gsea$database, "GOBP")
})

test_that("stored-result helpers list and retrieve DE, enrichment, and interpretation outputs", {
  skip_if_not_installed("Seurat")

  provider <- function(messages, model = NULL, ...) {
    list(text = "mock response")
  }

  object <- make_interpretation_object()
  object <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    provider = provider,
    store_name = "annotation_note",
    return_object = TRUE
  )

  result_index <- sn_list_results(object)
  expect_true(all(c("collection", "type", "name") %in% colnames(result_index)))
  expect_true(any(result_index$name == "celltype_markers"))
  expect_true(any(result_index$name == "celltype_gsea"))
  expect_true(any(result_index$name == "annotation_note"))

  marker_top <- sn_get_de_result(object, de_name = "celltype_markers", top_n = 2)
  expect_true(nrow(marker_top) > 0)

  enrich_top <- sn_get_enrichment_result(object, enrichment_name = "celltype_gsea", top_n = 1)
  expect_equal(nrow(enrich_top), 1)

  interpretation <- sn_get_interpretation_result(object, interpretation_name = "annotation_note")
  expect_match(interpretation$response$text, "mock response")
})

test_that("annotation and DE evidence helpers return structured outputs", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  object$percent.mt <- ifelse(object$cell_type == "Tcell", 4, 18)
  object$nCount_RNA_qc <- ifelse(object$cell_type == "Tcell", "Passed", "Failed")
  object$scDblFinder.class <- ifelse(object$cell_type == "Tcell", "singlet", "doublet")

  annotation_evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    n_markers = 3,
    enrichment_name = "celltype_gsea",
    n_terms = 2,
    include_qc = TRUE
  )
  de_evidence <- sn_prepare_de_evidence(
    object = object,
    de_name = "celltype_markers",
    n_genes = 3
  )

  expect_equal(annotation_evidence$task, "annotation")
  expect_true(all(c("cluster", "n_cells", "top_markers") %in% colnames(annotation_evidence$cluster_summary)))
  expect_true("top_functions" %in% colnames(annotation_evidence$cluster_summary))
  expect_true("doublet_fraction" %in% colnames(annotation_evidence$qc_summary))
  expect_equal(de_evidence$task, "de")
  expect_true("top_markers" %in% names(de_evidence))
})

test_that("annotation evidence and interpretation resolve a default marker result when de_name is omitted", {
  skip_if_not_installed("Seurat")

  provider <- function(messages, model = NULL, ...) {
    list(
      text = '{"cluster_annotations":[{"cluster":"Tcell","primary_label":"T cell","broad_label":"lymphocyte","confidence":"high","status":"confident","alternatives":[],"supporting_markers":["CD3D"],"supporting_functions":["T cell activation"],"risk_flags":[],"note":"marker support","recommended_checks":[]}],"narrative_summary":"demo"}'
    )
  }

  object <- make_interpretation_object()

  evidence <- sn_prepare_annotation_evidence(
    object = object,
    cluster_by = "cell_type",
    n_markers = 3
  )
  expect_equal(evidence$source_de_name, "celltype_markers")

  object <- sn_interpret_annotation(
    object = object,
    cluster_by = "cell_type",
    provider = provider,
    store_name = "annotation_default_de",
    return_object = TRUE
  )

  stored <- sn_get_interpretation_result(object, "annotation_default_de")
  expect_equal(stored$evidence$source_de_name, "celltype_markers")
})

test_that("enrichment and results evidence helpers work from stored results", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()

  enrichment_evidence <- sn_prepare_enrichment_evidence(
    object = object,
    enrichment_name = "celltype_gsea",
    n_terms = 1
  )
  results_evidence <- sn_prepare_results_evidence(
    object = object,
    cluster_de_name = "celltype_markers",
    enrichment_name = "celltype_gsea",
    cluster_by = "cell_type",
    n_markers = 2,
    n_terms = 1
  )

  expect_equal(enrichment_evidence$task, "enrichment")
  expect_equal(nrow(enrichment_evidence$top_terms), 1)
  expect_equal(results_evidence$task, "results")
  expect_true("cluster_markers" %in% names(results_evidence))
  expect_true("enrichment_summary" %in% names(results_evidence))
})

test_that("sn_build_prompt creates a prompt bundle from evidence", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    n_markers = 3
  )

  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "annotation",
    style = "concise annotation note",
    audience = "scientist"
  )

  expect_equal(prompt$task, "annotation")
  expect_true(all(c("system", "user", "messages") %in% names(prompt)))
  expect_length(prompt$messages, 2)
  expect_match(prompt$user, "^# Interpretation Request")
  expect_match(prompt$user, "## Task Instructions")
  expect_match(prompt$user, "## Evidence")
})

test_that("annotation prompts include all clusters instead of truncating at eight rows", {
  evidence <- list(
    task = "annotation",
    cluster_summary = tibble::tibble(
      cluster = as.character(0:11),
      top_markers = paste0("GENE", 0:11)
    ),
    top_marker_table = tibble::tibble(
      cluster = rep(as.character(0:11), each = 2),
      gene = paste0("GENE", seq_len(24)),
      avg_log2FC = seq_len(24)
    )
  )

  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "annotation",
    style = "cell type annotation",
    audience = "scientist"
  )

  expect_match(prompt$user, "\\b11\\b")
  expect_match(prompt$user, "GENE24")
})

test_that("sn_build_prompt supports human-readable output with background context", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    n_markers = 3
  )

  prompt <- sn_build_prompt(
    evidence = evidence,
    task = "annotation",
    background = "This dataset profiles PBMCs from a treatment study.",
    output_format = "human"
  )

  expect_equal(prompt$output_format, "human")
  expect_match(prompt$text, "Background")
  expect_match(prompt$text, "PBMCs")
})

test_that("high-level interpretation helpers can return prompts or store provider output", {
  skip_if_not_installed("Seurat")

  provider <- function(messages, model = NULL, ...) {
    list(text = paste("mock", model %||% "model", length(messages)))
  }

  object <- make_interpretation_object()

  prompt <- sn_write_results(
    object = object,
    cluster_de_name = "celltype_markers",
    enrichment_name = "celltype_gsea",
    cluster_by = "cell_type",
    background = "PBMC treatment study",
    return_prompt = TRUE
  )
  expect_equal(prompt$task, "results")

  human_prompt <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    background = "PBMC treatment study",
    output_format = "human"
  )
  expect_equal(human_prompt$output_format, "human")
  expect_match(human_prompt$text, "cell type annotation")

  object <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    provider = provider,
    model = "demo-model",
    store_name = "annotation_note",
    return_object = TRUE
  )

  expect_true("interpretation_results" %in% names(object@misc))
  expect_true("annotation_note" %in% names(object@misc$interpretation_results))
  expect_match(object@misc$interpretation_results$annotation_note$response$text, "demo-model")
})

test_that("structured annotation responses are normalized and written back to metadata", {
  skip_if_not_installed("Seurat")

  provider <- function(messages, model = NULL, ...) {
    list(text = paste(
      '{"cluster_annotations":[',
      '{"cluster":"Tcell","primary_label":"activated t cell","broad_label":"lymphocyte","confidence":"high","status":"confident","alternatives":["cytotoxic t cell"],',
      '"supporting_markers":["CD3D","TRAC"],"supporting_functions":["T cell activation"],"risk_flags":["transitional_state"],',
      '"note":"marker and pathway support a T-cell identity","recommended_checks":["confirm IL7R and LTB"]},',
      '{"cluster":"Bcell","primary_label":"b cell","broad_label":"lymphocyte","confidence":"medium","status":"possible_contamination","alternatives":["plasma cell"],',
      '"supporting_markers":["MS4A1","CD79A"],"supporting_functions":["B cell receptor signaling"],"risk_flags":["contamination"],',
      '"note":"mixed B-cell and stress-like signals","recommended_checks":["inspect percent.mt and doublet score"]}',
      '],"narrative_summary":"demo summary"}'
    ))
  }

  object <- make_interpretation_object()
  object$percent.mt <- ifelse(object$cell_type == "Tcell", 4, 22)
  object$nCount_RNA_qc <- ifelse(object$cell_type == "Tcell", "Passed", "Failed")
  object$scDblFinder.class <- ifelse(object$cell_type == "Tcell", "singlet", "doublet")

  object <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    enrichment_name = "celltype_gsea",
    cluster_by = "cell_type",
    provider = provider,
    model = "demo-model",
    store_name = "annotation_structured",
    metadata_prefix = "shennong_celltype",
    return_object = TRUE
  )

  stored <- sn_get_interpretation_result(object, "annotation_structured")
  expect_true(is.data.frame(stored$annotation_table))
  expect_true(all(c("cluster", "primary_label", "confidence", "status") %in% colnames(stored$annotation_table)))
  expect_equal(stored$narrative_summary, "demo summary")
  expect_true(all(c(
    "shennong_celltype_label",
    "shennong_celltype_broad_label",
    "shennong_celltype_confidence",
    "shennong_celltype_status",
    "shennong_celltype_risk_flags"
  ) %in% colnames(object[[]])))
  expect_false(any(c(
    "shennong_celltype_supporting_markers",
    "shennong_celltype_supporting_functions",
    "shennong_celltype_note",
    "shennong_celltype_recommended_checks"
  ) %in% colnames(object[[]])))
  expect_equal(unique(stats::na.omit(as.character(object$shennong_celltype_label))), c("Activated T Cell", "B Cell"))
  expect_true("possible_contamination" %in% unique(stats::na.omit(as.character(object$shennong_celltype_status))))
})

test_that("structured annotation parsing works from ellmer-native structured payloads", {
  response <- list(
    structured = list(
      cluster_annotations = tibble::tibble(
        cluster = c("0", "1"),
        primary_label = c("ILC2-like", "KIT+ ILC-like / ILCP-like"),
        broad_label = c("Innate lymphoid", "Innate lymphoid"),
        confidence = c("high", "medium"),
        status = c("confident", "ambiguous"),
        risk_flags = I(list(character(), c("transitional_state"))),
        supporting_markers = I(list(c("RORA", "GATA3"), c("KIT", "IL7R"))),
        note = c("type 2 program", "precursor-like state")
      ),
      narrative_summary = "demo structured summary"
    )
  )

  parsed <- .sn_parse_annotation_response(response)

  expect_true(is.data.frame(parsed$table))
  expect_equal(parsed$narrative, "demo structured summary")
  expect_equal(parsed$table$cluster, c("0", "1"))
  expect_equal(parsed$table$supporting_markers[[2]], "KIT; IL7R")
})

test_that("sn_make_ellmer_provider retries retryable proxy errors", {
  skip_if_not_installed("ellmer")

  attempts <- 0
  mock_chat <- list(
    chat = function(prompt_text) {
      attempts <<- attempts + 1
      if (attempts == 1) {
        stop("lexical error: invalid char in json text. <!DOCTYPE html>")
      }
      "OK"
    }
  )

  provider <- with_mocked_bindings(
    sn_make_ellmer_provider(
      api_key = "demo-key",
      base_url = "https://example.invalid/v1",
      model = "demo-model",
      retries = 2,
      retry_delay_sec = 0
    ),
    chat_openai_compatible = function(...) mock_chat,
    .package = "ellmer"
  )

  response <- with_mocked_bindings(
    provider(messages = list(list(role = "user", content = "hello"))),
    chat_openai_compatible = function(...) mock_chat,
    .package = "ellmer"
  )

  expect_equal(response$text, "OK")
  expect_equal(attempts, 2)
})

test_that("sn_make_ellmer_provider forwards reasoning_effort to ellmer api_args", {
  skip_if_not_installed("ellmer")

  captured <- NULL
  mock_chat <- list(
    chat = function(prompt_text) {
      "OK"
    }
  )

  provider <- with_mocked_bindings(
    sn_make_ellmer_provider(
      api_key = "demo-key",
      base_url = "https://example.invalid/v1",
      model = "demo-model",
      reasoning_effort = "high",
      retries = 1
    ),
    chat_openai_compatible = function(...) {
      captured <<- list(...)
      mock_chat
    },
    .package = "ellmer"
  )

  response <- with_mocked_bindings(
    provider(messages = list(list(role = "user", content = "hello"))),
    chat_openai_compatible = function(...) {
      captured <<- list(...)
      mock_chat
    },
    .package = "ellmer"
  )

  expect_equal(response$text, "OK")
  expect_equal(captured$api_args$reasoning_effort, "high")
})

test_that("sn_make_ellmer_provider uses ellmer structured output and registers tools", {
  skip_if_not_installed("ellmer")

  captured <- list()
  mock_chat <- list(
    register_tool = function(tool_def) {
      captured$tools <<- c(captured$tools %||% list(), list(tool_def))
      invisible(NULL)
    },
    chat = function(prompt_text) {
      captured$chat_prompt <<- prompt_text
      "comparison summary"
    },
    chat_structured = function(prompt_text, type, echo = "none", convert = TRUE) {
      captured$structured_prompt <<- prompt_text
      captured$structured_type <<- type
      list(
        cluster_annotations = tibble::tibble(
          cluster = "0",
          primary_label = "ILC2-like",
          broad_label = "Innate lymphoid",
          confidence = "high",
          status = "confident"
        ),
        narrative_summary = "summary"
      )
    }
  )

  provider <- with_mocked_bindings(
    sn_make_ellmer_provider(
      api_key = "demo-key",
      base_url = "https://example.invalid/v1",
      model = "demo-model",
      retries = 1
    ),
    chat_openai_compatible = function(...) mock_chat,
    .package = "ellmer"
  )

  tool_def <- ellmer::tool(
    function(cluster_ids) cluster_ids,
    name = "lookup_cluster_annotation_evidence",
    description = "demo tool",
    arguments = list(
      cluster_ids = ellmer::type_array(ellmer::type_string("cluster id"))
    )
  )

  tool_response <- with_mocked_bindings(
    provider(
      messages = list(list(role = "user", content = "hello")),
      tools = list(tool_def)
    ),
    chat_openai_compatible = function(...) mock_chat,
    .package = "ellmer"
  )
  structured_response <- with_mocked_bindings(
    provider(
      messages = list(list(role = "user", content = "hello")),
      structured_type = .sn_annotation_structured_type()
    ),
    chat_openai_compatible = function(...) mock_chat,
    .package = "ellmer"
  )

  expect_equal(tool_response$text, "comparison summary")
  expect_length(captured$tools, 1)
  expect_true(!is.null(structured_response$structured))
  expect_match(structured_response$text, "cluster_annotations")
})

test_that("sn_interpret_annotation falls back to the default ellmer provider path", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  object$percent.mt <- ifelse(object$cell_type == "Tcell", 4, 22)
  object$nCount_RNA_qc <- ifelse(object$cell_type == "Tcell", "Passed", "Failed")
  object$scDblFinder.class <- ifelse(object$cell_type == "Tcell", "singlet", "doublet")

  object <- with_mocked_bindings(
    sn_interpret_annotation(
      object = object,
      de_name = "celltype_markers",
      enrichment_name = "celltype_gsea",
      cluster_by = "cell_type",
      store_name = "annotation_structured_default",
      metadata_prefix = "shennong_celltype",
      return_object = TRUE
    ),
    .sn_get_default_ellmer_provider = function() function(messages, model = NULL, ...) {
      list(text = paste(
        '{"cluster_annotations":[',
        '{"cluster":"Tcell","primary_label":"activated t cell","broad_label":"lymphocyte","confidence":"high","status":"confident","alternatives":["cytotoxic t cell"],',
        '"supporting_markers":["CD3D","TRAC"],"supporting_functions":["T cell activation"],"risk_flags":["transitional_state"],',
        '"note":"marker and pathway support a T-cell identity","recommended_checks":["confirm IL7R and LTB"]},',
        '{"cluster":"Bcell","primary_label":"b cell","broad_label":"lymphocyte","confidence":"medium","status":"possible_contamination","alternatives":["plasma cell"],',
        '"supporting_markers":["MS4A1","CD79A"],"supporting_functions":["B cell receptor signaling"],"risk_flags":["contamination"],',
        '"note":"mixed B-cell and stress-like signals","recommended_checks":["inspect percent.mt and doublet score"]}',
        '],"narrative_summary":"demo summary"}'
      ))
    },
    .package = "Shennong"
  )

  expect_true("annotation_structured_default" %in% names(object@misc$interpretation_results))
  expect_true("shennong_celltype_label" %in% colnames(object[[]]))
})

test_that("annotation prompt can include candidate labels for sorted datasets", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()

  prompt <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    background = "Blood ILC-sorted dataset.",
    label_candidates = c("ILC1", "ILC2", "ILC3", "NK", "T cell"),
    return_prompt = TRUE
  )

  expect_match(prompt$user, "Annotation priors / candidate labels")
  expect_match(prompt$user, "ILC1")
  expect_match(prompt$user, "do not force unrelated T-cell or NK-cell labels")
})

test_that("annotation metadata_fields can opt into detailed metadata columns", {
  skip_if_not_installed("Seurat")

  provider <- function(messages, model = NULL, ...) {
    list(text = paste(
      '{"cluster_annotations":[',
      '{"cluster":"Tcell","primary_label":"activated t cell","broad_label":"lymphocyte","confidence":"high","status":"confident","alternatives":["cytotoxic t cell"],',
      '"supporting_markers":["CD3D","TRAC"],"supporting_functions":["T cell activation"],"risk_flags":["transitional_state"],',
      '"note":"marker and pathway support a T-cell identity","recommended_checks":["confirm IL7R and LTB"]},',
      '{"cluster":"Bcell","primary_label":"b cell","broad_label":"lymphocyte","confidence":"medium","status":"possible_contamination","alternatives":["plasma cell"],',
      '"supporting_markers":["MS4A1","CD79A"],"supporting_functions":["B cell receptor signaling"],"risk_flags":["contamination"],',
      '"note":"mixed B-cell and stress-like signals","recommended_checks":["inspect percent.mt and doublet score"]}',
      '],"narrative_summary":"demo summary"}'
    ))
  }

  object <- make_interpretation_object()
  object <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    provider = provider,
    metadata_prefix = "custom_ann",
    metadata_fields = c("primary_label", "confidence", "supporting_markers"),
    store_name = "annotation_custom_metadata",
    return_object = TRUE
  )

  expect_true(all(c(
    "custom_ann_label",
    "custom_ann_confidence",
    "custom_ann_supporting_markers"
  ) %in% colnames(object[[]])))
  expect_false("custom_ann_broad_label" %in% colnames(object[[]]))
})

test_that("misc-result helpers create collections and report missing stored results", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_base_object()
  stored_object <- Shennong:::.sn_store_misc_result(
    object = object,
    collection = "custom_results",
    store_name = "demo",
    result = list(value = 1, label = "stored")
  )

  stored_result <- Shennong:::.sn_get_misc_result(
    object = stored_object,
    collection = "custom_results",
    store_name = "demo"
  )

  expect_equal(stored_result$value, 1)
  expect_equal(stored_result$label, "stored")
  expect_error(
    Shennong:::.sn_get_misc_result(
      object = stored_object,
      collection = "custom_results",
      store_name = "missing"
    ),
    "No stored result named 'missing'"
  )
})

test_that("annotation default DE-name resolution prefers the default stored result", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  object@misc$de_results$default <- object@misc$de_results$celltype_markers

  resolved <- Shennong:::.sn_resolve_misc_result_name(
    object = object,
    collection = "de_results",
    preferred_analysis = "markers",
    arg_name = "de_name"
  )

  expect_equal(resolved, "default")
})

test_that("stored-result retrieval supports filtering, ranking, and metadata return", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  grouped_terms <- tibble::tibble(
    Cluster = c("Bcell", "Bcell", "Tcell"),
    Description = c("best B term", "second B term", "best T term"),
    p.adjust = c(0.01, 0.05, 0.02)
  )
  object <- sn_store_enrichment(
    object = object,
    result = grouped_terms,
    store_name = "cluster_ora",
    analysis = "ora",
    database = "TESTDB",
    return_object = TRUE
  )

  top_tcell_markers <- sn_get_de_result(
    object = object,
    de_name = "celltype_markers",
    top_n = 2,
    direction = "up",
    groups = "Tcell"
  )
  de_metadata <- sn_get_de_result(
    object = object,
    de_name = "celltype_markers",
    with_metadata = TRUE
  )
  top_b_terms <- sn_get_enrichment_result(
    object = object,
    enrichment_name = "cluster_ora",
    top_n = 1,
    groups = "Bcell"
  )
  enrichment_metadata <- sn_get_enrichment_result(
    object = object,
    enrichment_name = "cluster_ora",
    with_metadata = TRUE
  )

  expect_true(all(top_tcell_markers$cluster == "Tcell"))
  expect_lte(nrow(top_tcell_markers), 2)
  expect_true(all(c("table", "analysis", "group_col", "rank_col") %in% names(de_metadata)))
  expect_equal(nrow(top_b_terms), 1)
  expect_equal(top_b_terms$Cluster[[1]], "Bcell")
  expect_equal(top_b_terms$Description[[1]], "best B term")
  expect_equal(enrichment_metadata$analysis, "ora")
  expect_equal(enrichment_metadata$database, "TESTDB")
})

test_that("ranked-table subsetting handles grouped up and down selection", {
  ranked <- tibble::tibble(
    cluster = c("A", "A", "B", "B"),
    gene = c("g1", "g2", "g3", "g4"),
    score = c(3, -2, 5, -4)
  )

  top_up <- Shennong:::.sn_subset_ranked_table(
    table = ranked,
    rank_col = "score",
    group_col = "cluster",
    top_n = 1,
    direction = "up"
  )
  top_down <- Shennong:::.sn_subset_ranked_table(
    table = ranked,
    rank_col = "score",
    group_col = "cluster",
    top_n = 1,
    direction = "down"
  )

  expect_equal(top_up$gene, c("g1", "g3"))
  expect_equal(top_down$gene, c("g2", "g4"))
})

test_that("cluster, marker, and prediction summaries keep object-derived structure", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  object$celltypist_predicted_labels <- ifelse(object$cell_type == "Tcell", "T lineage", "B lineage")
  object$ref_majority_voting <- ifelse(object$cell_type == "Tcell", "T consensus", "B consensus")

  de_result <- sn_get_de_result(
    object = object,
    de_name = "celltype_markers",
    with_metadata = TRUE
  )
  cluster_summary <- Shennong:::.sn_prepare_cluster_summary(object, cluster_col = "cell_type")
  marker_summary <- Shennong:::.sn_prepare_marker_summary(de_result, n_markers = 3)
  marker_table <- Shennong:::.sn_prepare_marker_table(de_result, n_markers = 2)
  prediction_summary <- Shennong:::.sn_prepare_prediction_summary(object, cluster_col = "cell_type")

  expect_true(all(c("cluster", "n_cells", "fraction", "sample_distribution") %in% colnames(cluster_summary)))
  expect_equal(nrow(marker_summary), 2)
  expect_true(all(c("cluster", "top_markers") %in% colnames(marker_summary)))
  expect_equal(nrow(marker_table), 4)
  expect_true(all(c("cluster", "celltypist_predicted_labels", "ref_majority_voting") %in% colnames(prediction_summary)))
  expect_true(all(grepl("\\(40 cells\\)", prediction_summary$celltypist_predicted_labels)))
})

test_that("specific marker selection prefers cluster-restricted markers over shared top hits", {
  de_result <- list(
    table = tibble::tibble(
      cluster = c("A", "A", "B", "B"),
      gene = c("SHARED", "A_SPEC", "SHARED", "B_SPEC"),
      avg_log2FC = c(6, 5, 6, 5),
      p_val_adj = c(1e-6, 1e-6, 1e-6, 1e-6)
    ),
    group_col = "cluster",
    rank_col = "avg_log2FC",
    p_col = "p_val_adj",
    p_val_cutoff = 0.05,
    analysis = "markers"
  )

  top_summary <- Shennong:::.sn_prepare_marker_summary(de_result, n_markers = 1, selection = "top")
  specific_summary <- Shennong:::.sn_prepare_marker_summary(de_result, n_markers = 1, selection = "specific")

  expect_equal(top_summary$top_markers, c("SHARED", "SHARED"))
  expect_equal(specific_summary$top_markers, c("A_SPEC", "B_SPEC"))
})

test_that("annotation evidence can include cluster neighborhood geometry", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  embeddings <- matrix(
    c(
      rnorm(40, mean = 0, sd = 0.1),
      rnorm(40, mean = 0, sd = 0.1),
      rnorm(40, mean = 3, sd = 0.1),
      rnorm(40, mean = 0, sd = 0.1)
    ),
    ncol = 2
  )
  rownames(embeddings) <- colnames(object)
  colnames(embeddings) <- c("UMAP_1", "UMAP_2")
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = "UMAP_",
    assay = Seurat::DefaultAssay(object)
  )

  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    reduction = "umap",
    n_neighbor_clusters = 1
  )

  expect_true("geometry_summary" %in% names(evidence))
  expect_equal(nrow(evidence$geometry_summary), 2)
  expect_true("nearest_clusters" %in% colnames(evidence$cluster_summary))
  expect_true(any(grepl("Tcell|Bcell", evidence$geometry_summary$nearest_clusters)))
})

test_that("annotation evidence adds canonical lineage heuristic hints", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    marker_selection = "specific",
    reduction = NULL
  )

  expect_true("lineage_hints" %in% names(evidence))
  expect_true(all(c("heuristic_hint", "heuristic_top_signatures") %in% colnames(evidence$lineage_hints)))
  expect_true(all(c("heuristic_hint", "heuristic_rationale") %in% colnames(evidence$cluster_summary)))
  expect_true(any(grepl("T-cell-like", evidence$lineage_hints$heuristic_hint)))
  expect_true(any(grepl("B-cell-like", evidence$lineage_hints$heuristic_hint)))
})

test_that("annotation evidence exposes a canonical marker snapshot", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  evidence <- sn_prepare_annotation_evidence(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    marker_selection = "specific",
    reduction = NULL
  )

  expect_true("canonical_marker_snapshot" %in% names(evidence))
  expect_true(all(c("cluster", "CD3D", "MS4A1") %in% colnames(evidence$canonical_marker_snapshot)))
})

test_that("annotation heuristics identify erythroid contamination and mixed transitional states", {
  erythroid_hint <- .sn_select_annotation_hint(
    score_row = c(erythroid_contam = 3, nk_cytotoxic = 0.2),
    marker_support = c("HBB", "AHSP", "HBA1"),
    feature_values = c(HBB = 8, AHSP = 6, HBA1 = 7, HBA2 = 5, KLF1 = 4)
  )
  mixed_hint <- .sn_select_annotation_hint(
    score_row = c(t_cell = 2.5, nk_cytotoxic = 2.7, ilc3 = 1.4),
    marker_support = c("CD3E", "TRAC", "NKG7", "GNLY"),
    feature_values = c(CD3E = 3.2, TRAC = 2.9, NKG7 = 4.1, GNLY = 5.3, KLRD1 = 2.4)
  )
  transitional_hint <- .sn_select_annotation_hint(
    score_row = c(ilc3 = 2.8, nk_cytotoxic = 2.6, t_cell = 0.6),
    marker_support = c("IL23R", "RORC", "NCR2", "NKG7", "GNLY"),
    feature_values = c(IL23R = 2.4, RORC = 1.2, AHR = 1.8, NCR2 = 1.7, NKG7 = 3.5, GNLY = 3.8)
  )

  expect_equal(erythroid_hint$hint, "Erythroid contamination-like")
  expect_equal(mixed_hint$hint, "T/NK mixed lymphoid-like")
  expect_equal(transitional_hint$hint, "NK/ILC3 transitional-like")
})

test_that("annotation background injects blood ILC and contamination priors", {
  background <- .sn_compose_annotation_background(
    background = "This is blood ILC-sorted single-cell RNA-seq data.",
    label_candidates = c("ILC1", "ILC2", "ILC3", "NK", "T cell")
  )

  expect_match(background, "Blood ILC prior", fixed = TRUE)
  expect_match(background, "ILCP-like", fixed = TRUE)
  expect_match(background, "erythroid contamination", ignore.case = TRUE)
})

test_that("sn_interpret_annotation agentic mode performs a focused refinement pass", {
  skip_if_not_installed("Seurat")

  calls <- character()
  provider <- function(messages, model = NULL, ...) {
    prompt_text <- paste(vapply(messages, function(x) x$content, character(1)), collapse = "\n\n")
    if (grepl("Annotation stage: broad pass", prompt_text, fixed = TRUE)) {
      calls <<- c(calls, "broad")
      return(list(text = paste(
        '{"cluster_annotations":[',
        '{"cluster":"Tcell","primary_label":"lymphoid lineage","broad_label":"lymphoid","confidence":"low","status":"ambiguous","alternatives":[],"supporting_markers":["CD3D"],"supporting_functions":[],"risk_flags":[],"note":"broad pass","recommended_checks":[]},',
        '{"cluster":"Bcell","primary_label":"lymphoid lineage","broad_label":"lymphoid","confidence":"low","status":"ambiguous","alternatives":[],"supporting_markers":["MS4A1"],"supporting_functions":[],"risk_flags":[],"note":"broad pass","recommended_checks":[]}',
        '],"narrative_summary":"broad"}'
      )))
    }
    if (grepl("Annotation stage: focused refinement", prompt_text, fixed = TRUE)) {
      calls <<- c(calls, "focused")
      expect_match(prompt_text, "canonical_marker_snapshot")
      return(list(text = paste(
        '{"cluster_annotations":[',
        '{"cluster":"Tcell","primary_label":"activated t cell","broad_label":"lymphoid","confidence":"high","status":"confident","alternatives":[],"supporting_markers":["CD3D","TRAC"],"supporting_functions":["T cell activation"],"risk_flags":[],"note":"focused pass","recommended_checks":[]},',
        '{"cluster":"Bcell","primary_label":"b cell","broad_label":"lymphoid","confidence":"high","status":"confident","alternatives":[],"supporting_markers":["MS4A1","CD79A"],"supporting_functions":["B cell receptor signaling"],"risk_flags":[],"note":"focused pass","recommended_checks":[]}',
        '],"narrative_summary":"focused"}'
      )))
    }
    stop("Unexpected prompt")
  }

  object <- make_interpretation_object()
  object <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    annotation_mode = "agentic",
    provider = provider,
    store_name = "annotation_agentic",
    metadata_prefix = "agentic_ann",
    return_object = TRUE,
    show_progress = FALSE
  )

  stored <- sn_get_interpretation_result(object, "annotation_agentic")
  expect_equal(calls, c("broad", "focused"))
  expect_equal(stored$annotation_mode, "agentic")
  expect_equal(stored$narrative_summary, "focused")
  expect_equal(sort(stored$workflow$focus_clusters), c("Bcell", "Tcell"))
  expect_equal(unique(stats::na.omit(as.character(object$agentic_ann_label))), c("Activated T Cell", "B Cell"))
})

test_that("agentic annotation can return its staged prompt bundle", {
  skip_if_not_installed("Seurat")

  object <- make_interpretation_object()
  prompt <- sn_interpret_annotation(
    object = object,
    de_name = "celltype_markers",
    cluster_by = "cell_type",
    annotation_mode = "agentic",
    return_prompt = TRUE
  )

  expect_equal(prompt$annotation_mode, "agentic")
  expect_true(is.list(prompt$broad_prompt))
  expect_match(prompt$note, "Focused refinement prompt")
})
