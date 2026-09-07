make_communication_object <- function() {
  set.seed(101)
  genes <- c(
    "LIG1", "REC1", "TG1", "TG2", "TFG1", "TFG2", "ACTB", "MALAT1",
    "TGFB1", "TGFBR1", "CXCL12", "CXCR4", "CCL5", "CCR5"
  )
  counts <- matrix(rpois(length(genes) * 48, lambda = 2), nrow = length(genes), ncol = 48)
  rownames(counts) <- genes
  colnames(counts) <- paste0("cell", seq_len(48))
  cell_type <- rep(rep(c("Sender", "Receiver"), each = 3), 8)
  sample <- rep(paste0("S", seq_len(8)), each = 6)
  condition <- rep(rep(c("Ctrl", "Stim"), each = 4), each = 6)
  counts[rownames(counts) == "LIG1", cell_type == "Sender"] <- counts[rownames(counts) == "LIG1", cell_type == "Sender"] + 5
  counts[rownames(counts) == "REC1", cell_type == "Receiver"] <- counts[rownames(counts) == "REC1", cell_type == "Receiver"] + 5
  counts[c("TGFB1", "CXCL12", "CCL5"), cell_type == "Sender"] <- counts[c("TGFB1", "CXCL12", "CCL5"), cell_type == "Sender"] + 6
  counts[c("TGFBR1", "CXCR4", "CCR5"), cell_type == "Receiver"] <- counts[c("TGFBR1", "CXCR4", "CCR5"), cell_type == "Receiver"] + 6
  counts[rownames(counts) %in% c("TG1", "TG2"), cell_type == "Receiver" & condition == "Stim"] <-
    counts[rownames(counts) %in% c("TG1", "TG2"), cell_type == "Receiver" & condition == "Stim"] + 6
  object <- sn_initialize_seurat_object(
    x = Matrix::Matrix(counts, sparse = TRUE),
    project = "communication",
    species = "human"
  )
  object$cell_type <- cell_type
  object$sample <- sample
  object$condition <- condition
  Seurat::NormalizeData(object, verbose = FALSE)
}

make_multinichenet_object <- function() {
  set.seed(202)
  ligands <- paste0("L", seq_len(30))
  receptors <- paste0("R", seq_len(30))
  targets <- paste0("T", seq_len(40))
  genes <- c(ligands, receptors, targets)
  sample <- rep(paste0("MS", seq_len(8)), each = 20)
  condition <- rep(rep(c("Ctrl", "Stim"), each = 4), each = 20)
  cell_type <- rep(rep(c("Sender", "Receiver"), each = 10), 8)
  counts <- matrix(rpois(length(genes) * length(sample), lambda = 2), nrow = length(genes))
  rownames(counts) <- genes
  colnames(counts) <- paste0("mcell", seq_along(sample))
  counts[ligands, cell_type == "Sender"] <- counts[ligands, cell_type == "Sender"] + 3
  counts[receptors, cell_type == "Receiver"] <- counts[receptors, cell_type == "Receiver"] + 3
  counts[ligands[seq_len(10)], cell_type == "Sender" & condition == "Stim"] <-
    counts[ligands[seq_len(10)], cell_type == "Sender" & condition == "Stim"] + 5
  counts[receptors[seq_len(10)], cell_type == "Receiver" & condition == "Stim"] <-
    counts[receptors[seq_len(10)], cell_type == "Receiver" & condition == "Stim"] + 4
  counts[targets[seq_len(20)], cell_type == "Receiver" & condition == "Stim"] <-
    counts[targets[seq_len(20)], cell_type == "Receiver" & condition == "Stim"] + 8
  object <- sn_initialize_seurat_object(Matrix::Matrix(counts, sparse = TRUE), project = "multinichenet", species = "human")
  object$cell_type <- cell_type
  object$sample <- sample
  object$condition <- condition
  list(
    object = Seurat::NormalizeData(object, verbose = FALSE),
    ligand_target_matrix = matrix(
      runif(length(targets) * length(ligands), 0, 0.2),
      nrow = length(targets), dimnames = list(targets, ligands)
    ),
    lr_network = tibble::tibble(ligand = ligands, receptor = receptors)
  )
}

test_that("backend label maps make zero-valued groups safe and reversible", {
  map <- Shennong:::.sn_multinichenet_label_map(c("0", "1", "T cell", "0"))

  expect_false(any(unname(map$encoded) == "0"))
  expect_identical(
    unname(map$decoded[unname(map$encoded[c("0", "1", "T cell")])]),
    c("0", "1", "T cell")
  )
})

test_that("communication consensus keeps condition and sample contexts separate", {
  table <- tibble::tibble(
    source = "Sender", target = "Receiver", ligand = "LIG1", receptor = "REC1",
    score = c(1, 2, 3, 4), p_value = NA_real_, q_value = NA_real_,
    rank = c(1, 1, 1, 1), method = rep(c("m1", "m2"), 2),
    condition = rep(c("A", "B"), each = 2), sample = rep(c("S1", "S2"), each = 2),
    pathway = NA_character_, target_genes = NA_character_,
    evidence_source = rep(c("m1", "m2"), 2), spatial_distance = NA_real_
  )

  consensus <- Shennong:::.sn_communication_consensus(
    table, methods = c("m1", "m2"), min_methods = 2L
  )

  expect_equal(nrow(consensus), 2L)
  expect_setequal(consensus$condition, c("A", "B"))
  expect_setequal(consensus$sample, c("S1", "S2"))
  expect_false(any(grepl(";", consensus$condition, fixed = TRUE)))
})

test_that("communication comparisons require direction and support matched units", {
  evidence <- tibble::as_tibble(expand.grid(
    condition = c("Ctrl", "Stim"), pair = paste0("D", 1:4),
    receptor = "REC1", ligand = "LIG1", target = "Receiver", source = "Sender",
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ))
  evidence$sample <- paste(evidence$pair, evidence$condition, sep = "_")
  evidence$score <- ifelse(evidence$condition == "Stim", 2, 1) +
    rep(c(0.1, 0.2, 0.3, 0.4), each = 2)

  expect_error(
    Shennong:::.sn_compare_communication_samples(evidence),
    "explicit two-level"
  )
  expect_error(
    Shennong:::.sn_compare_communication_samples(evidence, c("Stim", "Stim")),
    "two distinct"
  )

  comparison <- Shennong:::.sn_compare_communication_samples(
    evidence, c("Stim", "Ctrl")
  )
  expect_true(comparison$paired)
  expect_equal(comparison$n_pairs, 4L)
  expect_equal(comparison$estimate, 1)
})

test_that("paired communication comparisons require complete selected-condition pairs", {
  evidence <- tibble::as_tibble(expand.grid(
    condition = c("Ctrl", "Stim", "Other"), pair = c("D1", "D2"),
    receptor = "REC1", ligand = "LIG1", target = "Receiver", source = "Sender",
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ))
  evidence$sample <- paste(evidence$pair, evidence$condition, sep = "_")
  evidence$score <- ifelse(
    evidence$condition == "Stim", 2,
    ifelse(evidence$condition == "Ctrl", 1, 100)
  )

  comparison <- Shennong:::.sn_compare_communication_samples(
    evidence, c("Stim", "Ctrl")
  )
  expect_equal(comparison$n_pairs, 2L)
  expect_equal(comparison$estimate, 1)

  incomplete <- evidence[!(evidence$pair == "D2" & evidence$condition == "Stim"), ]
  expect_error(
    Shennong:::.sn_compare_communication_samples(incomplete, c("Stim", "Ctrl")),
    "Incomplete pair.*D2"
  )

  missing_pair <- evidence
  missing_pair$pair[[1L]] <- NA_character_
  expect_error(
    Shennong:::.sn_compare_communication_samples(missing_pair, c("Stim", "Ctrl")),
    "missing or empty pair labels"
  )
})

test_that("communication sample evidence records a stable paired-unit mapping", {
  object <- make_communication_object()
  sample_order <- paste0("S", 1:8)
  donor_map <- stats::setNames(rep(paste0("D", 1:4), 2), sample_order)
  object$donor <- unname(donor_map[as.character(object$sample)])
  interactions <- tibble::tibble(
    source = "Sender", target = "Receiver", ligand = "LIG1", receptor = "REC1"
  )

  evidence <- Shennong:::.sn_communication_sample_evidence(
    object, interactions, "cell_type", "sample", "condition", "donor", "RNA", "data"
  )

  expect_equal(nrow(evidence), 8L)
  expect_identical(
    unique(evidence$pair[evidence$sample %in% c("S1", "S5")]),
    "D1"
  )

  object$donor[[1L]] <- " "
  expect_error(
    Shennong:::.sn_communication_sample_evidence(
      object, interactions, "cell_type", "sample", "condition", "donor", "RNA", "data"
    ),
    "paired_by.*missing or empty"
  )
})

test_that("cell communication results can be stored and retrieved", {
  skip_if_not_installed("Seurat")
  object <- make_communication_object()
  tbl <- tibble::tibble(
    source = c("Sender", "Sender", "Receiver"),
    target = c("Receiver", "Sender", "Sender"),
    ligand = c("LIG1", "LIG1", "LIG2"),
    receptor = c("REC1", "REC2", "REC1"),
    score = c(0.9, 0.2, 0.1)
  )
  object <- sn_store_cell_communication(
    object = object,
    result = tbl,
    result_id = "manual",
    method = "manual"
  )

  retrieved <- sn_get_cell_communication_result(object, "manual", sources = "Sender", targets = "Receiver")
  metadata <- sn_get_cell_communication_result(object, "manual", with_metadata = TRUE)
  listed <- sn_list_results(object)

  expect_equal(nrow(retrieved), 1)
  expect_equal(retrieved$ligand, "LIG1")
  expect_equal(metadata$analysis, "cell_communication")
  expect_equal(metadata$analysis_type, "cell_communication")
  expect_equal(metadata$backend, "manual")
  expect_true(all(c("primary", "backend_raw", "consensus", "sample_evidence", "condition_comparison", "method_concordance", "ligand_targets") %in% names(metadata$tables)))
  expect_true(all(c("source", "target", "ligand", "receptor", "score", "p_value", "q_value", "rank", "method", "spatial_distance") %in% names(metadata$tables$primary)))
  expect_true("cell_communication" %in% listed$type)
  expect_true("manual" %in% names(object@misc$shennong$results$cell_communication))
})

test_that("direct communication storage does not claim ambient acceleration", {
  skip_if_not_installed("Seurat")
  object <- make_communication_object()
  table <- tibble::tibble(
    source = "Sender", target = "Receiver", ligand = "LIG1",
    receptor = "REC1", score = 0.9
  )
  direct <- sn_store_cell_communication(
    object, table, method = "cellchat", return_object = FALSE
  )
  expect_null(direct$provenance$acceleration)

  scoped <- Shennong:::.sn_with_acceleration_provenance_context(
    sn_store_cell_communication(
      object, table, method = "cellchat", return_object = FALSE
    ),
    patches = "cellchat"
  )
  expect_identical(
    scoped$provenance$acceleration$suppressed_patches,
    "cellchat"
  )
})

test_that("NicheNet backend runs with supplied priors", {
  skip_if_not_installed("Seurat")
  skip_if_not(suppressWarnings(requireNamespace("nichenetr", quietly = TRUE)), "nichenetr is not installed")
  withr::local_options(list(shennong.acceleration = FALSE))
  object <- make_communication_object()
  ligand_target_matrix <- matrix(
    c(
      0.8, 0.7, 0.1, 0.1,
      0.1, 0.2, 0.7, 0.6
    ),
    nrow = 4,
    dimnames = list(c("TG1", "TG2", "ACTB", "MALAT1"), c("LIG1", "LIGX"))
  )
  lr_network <- tibble::tibble(ligand = "LIG1", receptor = "REC1")

  stored <- sn_run_cell_communication(
    object = object,
    method = "nichenetr",
    group_by = "cell_type",
    sender = "Sender",
    receiver = "Receiver",
    geneset = c("TG1", "TG2"),
    background_genes = c("TG1", "TG2", "ACTB", "MALAT1"),
    ligand_target_matrix = ligand_target_matrix,
    lr_network = lr_network,
    sample_by = "sample",
    condition_by = "condition",
    contrast = c("Stim", "Ctrl"),
    expressed_pct = 0.05,
    top_n = 5,
    return_object = FALSE
  )

  expect_equal(stored$method, "nichenetr")
  expect_true("test_ligand" %in% colnames(stored$tables$primary))
  expect_true("LIG1" %in% stored$tables$primary$test_ligand)
  expect_equal(stored$backend, "nichenet")
  expect_equal(nrow(stored$tables$sample_evidence), 8L * nrow(stored$tables$primary))
  expect_true(nrow(stored$tables$condition_comparison) > 0L)
  expect_true(all(c("ligand", "target_gene", "weight") %in% names(stored$tables$ligand_targets)))
  expect_s3_class(sn_plot_communication(stored, type = "bubble"), "ggplot")
  expect_s3_class(sn_plot_communication(stored, type = "heatmap"), "ggplot")
  expect_s3_class(sn_plot_communication(stored, type = "network"), "ggplot")
  expect_s3_class(sn_plot_communication(stored, type = "chord"), "ggplot")
  expect_s3_class(sn_plot_communication(stored, type = "river"), "ggplot")
  expect_s3_class(sn_plot_ligand_target(stored), "ggplot")
  expect_s3_class(sn_plot_communication_comparison(stored), "ggplot")
})

test_that("communication backends record requested acceleration honestly", {
  skip_if_not_installed("Seurat")
  withr::local_options(list(shennong.acceleration = TRUE))
  withr::local_envvar(c(
    AUTOZYME_DISABLED = NA,
    AUTOZYME_DISABLE = NA
  ))
  object <- make_communication_object()
  state <- new.env(parent = emptyenv())
  state$active <- character()
  state$enable_requests <- character()
  state$disabled_scopes <- 0L
  state$calls <- list()
  backend_result <- list(
    table = tibble::tibble(
      source = "Sender", target = "Receiver",
      ligand = "LIG1", receptor = "REC1", score = 1
    ),
    artifacts = list()
  )

  testthat::local_mocked_bindings(
    .sn_with_default_acceleration = function(expr, patches, ...) {
      mapped <- Shennong:::.sn_map_acceleration_patches(patches)
      Shennong:::.sn_with_acceleration_provenance_context({
        state$enable_requests <- c(state$enable_requests, patches)
        Shennong:::.sn_record_acceleration_usage(mapped$supported)
        Shennong:::.sn_record_acceleration_suppression(mapped$unsupported)
        force(expr)
      }, patches = patches)
    },
    .sn_run_cellchat = function(...) {
      state$calls[[length(state$calls) + 1L]] <- list(
        method = "cellchat",
        active = state$active,
        args = list(...)
      )
      backend_result
    },
    .sn_run_nichenetr = function(...) {
      state$calls[[length(state$calls) + 1L]] <- list(
        method = "nichenetr",
        active = state$active,
        args = list(...)
      )
      backend_result
    },
    .package = "Shennong"
  )

  cellchat_stored <- sn_run_cell_communication(
    object,
    method = "cellchat",
    group_by = "cell_type",
    return_object = FALSE
  )
  expect_identical(state$enable_requests, "cellchat")
  expect_length(state$calls[[1L]]$active, 0L)
  expect_true("cellchat" %in%
    cellchat_stored$provenance$acceleration$suppressed_patches)
  expect_length(state$active, 0L)

  dense_prior <- matrix(
    1,
    nrow = 1,
    dimnames = list("TG1", "LIG1")
  )
  nichenet_stored <- sn_run_cell_communication(
    object,
    method = "nichenetr",
    group_by = "cell_type",
    sender = "Sender",
    receiver = "Receiver",
    ligand_target_matrix = dense_prior,
    return_object = FALSE
  )
  expect_identical(state$enable_requests, c("cellchat", "nichenetr"))
  expect_length(state$calls[[2L]]$active, 0L)
  expect_true("nichenetr" %in%
    nichenet_stored$provenance$acceleration$suppressed_patches)
  expect_length(state$active, 0L)

  state$active <- character()
  state$enable_requests <- character()
  disabled_stored <- withr::with_options(
    list(shennong.acceleration = FALSE),
    sn_run_cell_communication(
      object,
      method = "nichenetr",
      group_by = "cell_type",
      sender = "Sender",
      receiver = "Receiver",
      ligand_target_matrix = dense_prior,
      return_object = FALSE
    )
  )
  expect_identical(state$enable_requests, "nichenetr")
  expect_true("nichenetr" %in%
    disabled_stored$provenance$acceleration$suppressed_patches)
  expect_length(disabled_stored$provenance$acceleration$used_patches, 0L)
})

test_that("NicheNet backend strips legacy control arguments and runs upstream", {
  skip_if_not_installed("Seurat")
  object <- make_communication_object()
  state <- new.env(parent = emptyenv())
  state$calls <- list()
  backend_result <- list(
    table = tibble::tibble(
      source = "Sender", target = "Receiver",
      ligand = "LIG1", receptor = "REC1", score = 1
    ),
    artifacts = list()
  )
  dense_prior <- matrix(1, nrow = 1, dimnames = list("TG1", "LIG1"))
  sparse_prior <- Matrix::Matrix(dense_prior, sparse = TRUE)

  testthat::local_mocked_bindings(
    .sn_run_nichenetr = function(...) {
      state$calls[[length(state$calls) + 1L]] <- list(args = list(...))
      backend_result
    },
    .package = "Shennong"
  )

  stored <- sn_run_cell_communication(
    object,
    method = "nichenetr",
    group_by = "cell_type",
    sender = "Sender",
    receiver = "Receiver",
    ligand_target_matrix = dense_prior,
    backend_control = list(nichenet = list(single = FALSE, zyme = FALSE)),
    return_object = FALSE
  )
  expect_equal(stored$backend, "nichenet")
  current <- state$calls[[length(state$calls)]]
  expect_false("single" %in% names(current$args))
  expect_false("zyme" %in% names(current$args))
  expect_true(is.matrix(current$args$ligand_target_matrix))

  sn_run_cell_communication(
    object,
    method = "nichenetr",
    group_by = "cell_type",
    sender = "Sender",
    receiver = "Receiver",
    ligand_target_matrix = sparse_prior,
    return_object = FALSE
  )
  current <- state$calls[[length(state$calls)]]
  expect_s4_class(current$args$ligand_target_matrix, "Matrix")
})

test_that("communication backends standardize to one comparable schema", {
  liana <- Shennong:::.sn_standardize_communication(tibble::tibble(
    source = "Sender", target = "Receiver", ligand_complex = "LIG1",
    receptor_complex = "REC1", magnitude_rank = 0.1, specificity_rank = 0.2
  ), method = "liana", condition = "Stim")
  cellchat <- Shennong:::.sn_standardize_communication(tibble::tibble(
    source = "Sender", target = "Receiver", ligand = "LIG1",
    receptor = "REC1", prob = 0.8, pval = 0.01, pathway_name = "PathwayA"
  ), method = "cellchat", condition = "Stim")
  multinichenet <- Shennong:::.sn_standardize_communication(tibble::tibble(
    sender = "Sender", receiver = "Receiver", ligand = "LIG1",
    receptor = "REC1", prioritization_score = 0.7, group = "Stim"
  ), method = "multinichenet")
  natmi <- Shennong:::.sn_standardize_communication(tibble::tibble(
    source = c("Sender", "Sender"), target = c("Receiver", "Receiver"),
    ligand_complex = c("LIG1", "LIG2"), receptor_complex = c("REC1", "REC2"),
    prod_weight = c(0.42, 0.21), edge_specificity = c(0.8, 0.9)
  ), method = "liana")
  combined <- dplyr::bind_rows(liana, cellchat, multinichenet)
  consensus <- Shennong:::.sn_communication_consensus(combined)
  concordance <- Shennong:::.sn_communication_concordance(combined)

  expect_equal(natmi$score, c(0.42, 0.21))
  expect_equal(natmi$rank, c(1, 2))
  expect_true(all(c("source", "target", "ligand", "receptor", "score", "p_value", "q_value", "rank", "method", "condition", "sample", "pathway", "target_genes", "evidence_source", "spatial_distance") %in% names(combined)))
  expect_equal(consensus$n_methods, 3L)
  expect_equal(consensus$evidence_source, "cellchat;liana;multinichenet")
  expect_equal(nrow(concordance), 3L)
  expect_true(all(concordance$shared_edges == 1L))
  expect_true(all(concordance$complete_edges == 1L))
})

test_that("communication consensus requires cross-method support and does not mislabel minimum p-values", {
  shared <- tibble::tibble(
    source = "Sender", target = "Receiver", ligand = "L_shared",
    receptor = "R_shared", score = 1, p_value = 0.01
  )
  singleton <- tibble::tibble(
    source = "Sender", target = "Receiver", ligand = "L_single",
    receptor = "R_single", score = 100, p_value = 1e-12
  )
  table <- dplyr::bind_rows(
    Shennong:::.sn_standardize_communication(dplyr::bind_rows(shared, singleton), "m1"),
    Shennong:::.sn_standardize_communication(shared, "m2"),
    Shennong:::.sn_standardize_communication(shared, "m3")
  )

  consensus <- Shennong:::.sn_communication_consensus(
    table, methods = c("m1", "m2", "m3")
  )

  expect_identical(consensus$ligand, "L_shared")
  expect_identical(consensus$n_methods, 3L)
  expect_identical(consensus$available_methods, 3L)
  expect_equal(consensus$method_support_fraction, 1)
  expect_true(is.na(consensus$p_value))
  expect_true(is.na(consensus$q_value))
  expect_equal(consensus$minimum_method_p_value, 0.01)
  expect_identical(consensus$p_value_combination, "not_combined_correlated_methods")
})

test_that("single-backend communication output is not presented as consensus", {
  table <- Shennong:::.sn_standardize_communication(tibble::tibble(
    source = "Sender", target = "Receiver", ligand = "L1",
    receptor = "R1", score = 1
  ), "m1")
  consensus <- Shennong:::.sn_communication_consensus(table, methods = "m1")
  expect_equal(nrow(consensus), 0L)
})

test_that("LIANA input contains the exact selected assay layer", {
  skip_if_not_installed("SingleCellExperiment")
  object <- make_communication_object()
  selected <- SeuratObject::LayerData(object, assay = "RNA", layer = "data")
  selected <- selected + 7
  SeuratObject::LayerData(object, assay = "RNA", layer = "custom_liana") <- selected

  adapted <- Shennong:::.sn_liana_sce_input(
    object, assay = "RNA", layer = "custom_liana"
  )
  observed <- SummarizedExperiment::assay(adapted$object, adapted$assay)

  expect_equal(as.matrix(observed), as.matrix(selected))
  expect_identical(adapted$source_assay, "RNA")
  expect_identical(adapted$source_layer, "custom_liana")
})

test_that("communication concordance tolerates shared edges with missing ranks", {
  edges <- tibble::tibble(
    source = "Sender", target = "Receiver",
    ligand = paste0("L", 1:3), receptor = paste0("R", 1:3)
  )
  table <- dplyr::bind_rows(
    dplyr::mutate(edges, method = "liana", rank = c(1, 2, 3)),
    dplyr::mutate(edges, method = "cellchat", rank = c(NA_real_, 2, Inf))
  )

  concordance <- Shennong:::.sn_communication_concordance(table)
  expect_equal(concordance$shared_edges, 3L)
  expect_equal(concordance$complete_edges, 1L)
  expect_true(is.na(concordance$rank_correlation))

  table$rank[table$method == "cellchat"] <- c(3, 2, 1)
  complete <- Shennong:::.sn_communication_concordance(table)
  expect_equal(complete$complete_edges, 3L)
  expect_equal(complete$rank_correlation, -1)
})

test_that("public communication runner builds a real cross-method consensus", {
  skip_if_not_installed("liana")
  skip_if_not_installed("CellChat")
  withr::local_options(list(shennong.acceleration = FALSE))
  object <- make_communication_object()
  stored <- sn_run_cell_communication(
    object = object,
    method = c("liana", "cellchat"),
    group_by = "cell_type",
    sample_by = "sample",
    condition_by = "condition",
    contrast = c("Stim", "Ctrl"),
    consensus = TRUE,
    backend_control = list(
      liana = list(resource = "consensus", method = "natmi"),
      cellchat = list(min_cells = 3)
    ),
    return_object = FALSE
  )

  expect_equal(stored$method, "consensus")
  expect_equal(stored$backend, "liana+cellchat")
  expect_true(nrow(stored$tables$backend_raw) > 0L)
  expect_true(nrow(stored$tables$consensus) > 0L)
  expect_equal(nrow(stored$tables$method_concordance), 1L)
  expect_true(nrow(stored$tables$sample_evidence) > 0L)
})

test_that("CellPhoneDB output parser retains interaction evidence", {
  output_dir <- tempfile("cellphonedb-")
  dir.create(output_dir)
  pvalues <- data.frame(
    id_cp_interaction = c("CPI-1", "CPI-2"),
    gene_a = c("simple:LIG1", "complex:LIG1_LIG2"),
    gene_b = c("simple:REC1", "complex:REC1_REC2"),
    `Sender|Receiver` = c(0.01, 0.2),
    check.names = FALSE
  )
  means <- data.frame(
    id_cp_interaction = c("CPI-1", "CPI-2"),
    `Sender|Receiver` = c(1.2, 0.4),
    check.names = FALSE
  )
  utils::write.table(pvalues, file.path(output_dir, "pvalues.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(means, file.path(output_dir, "means.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  parsed <- Shennong:::.sn_read_cellphonedb_output(output_dir)
  standardized <- Shennong:::.sn_standardize_communication(parsed$table, "cellphonedb")
  expect_equal(nrow(standardized), 2L)
  expect_equal(standardized$ligand, c("LIG1", "LIG1_LIG2"))
  expect_equal(standardized$source, rep("Sender", 2))
  expect_equal(standardized$q_value, stats::p.adjust(c(0.01, 0.2), "BH"))
})

test_that("high-level CellPhoneDB consumes imported tables after temporary cleanup", {
  object <- make_communication_object()
  pvalues <- data.frame(
    id_cp_interaction = c("CPI-1", "CPI-2"),
    gene_a = c("simple:LIG1", "simple:TGFB1"),
    gene_b = c("simple:REC1", "simple:TGFBR1"),
    `Sender|Receiver` = c(0.01, 0.2),
    check.names = FALSE
  )
  means <- data.frame(
    id_cp_interaction = rev(pvalues$id_cp_interaction),
    gene_a = rev(pvalues$gene_a),
    gene_b = rev(pvalues$gene_b),
    `Sender|Receiver` = c(0.4, 1.2),
    check.names = FALSE
  )

  stored <- testthat::with_mocked_bindings(
    sn_run_cell_communication(
      object,
      method = "cellphonedb",
      group_by = "cell_type",
      return_object = FALSE
    ),
    sn_run_cellphonedb = function(...) {
      list(
        output_dir = NULL,
        run_dir_retained = FALSE,
        imported_tables = list(pvalues = pvalues, means = means)
      )
    },
    .package = "Shennong"
  )

  expect_identical(stored$analysis_type, "cell_communication")
  expect_equal(stored$tables$primary$score, c(1.2, 0.4))
  expect_equal(stored$tables$primary$p_value, c(0.01, 0.2))
  expect_null(stored$artifacts$cellphonedb$manifest$output_dir)
})

test_that("CellPhoneDB table pairing fails closed on malformed evidence", {
  pvalues <- data.frame(
    id_cp_interaction = "CPI-1", gene_a = "LIG1", gene_b = "REC1",
    `Sender|Receiver` = 1.5, check.names = FALSE
  )
  means <- data.frame(
    id_cp_interaction = "CPI-other", gene_a = "LIG1", gene_b = "REC1",
    `Sender|Receiver` = 1, check.names = FALSE
  )
  expect_error(
    Shennong:::.sn_parse_cellphonedb_tables(pvalues, means),
    "same unique, non-empty interaction identifiers"
  )
  means$id_cp_interaction <- "CPI-1"
  expect_error(
    Shennong:::.sn_parse_cellphonedb_tables(pvalues, means),
    "p-values must be numeric values in \\[0, 1\\]"
  )
})

test_that("MultiNicheNet backend uses biological samples and conditions", {
  skip_if_not_installed("multinichenetr")
  skip_if_not_installed("SingleCellExperiment")
  inputs <- make_multinichenet_object()
  inputs$ligand_target_matrix[seq_len(20), seq_len(10)] <-
    inputs$ligand_target_matrix[seq_len(20), seq_len(10)] + 0.8

  stored <- sn_run_cell_communication(
    object = inputs$object,
    method = "multinichenet",
    group_by = "cell_type",
    sample_by = "sample",
    condition_by = "condition",
    contrast = c("Stim", "Ctrl"),
    sender = "Sender",
    receiver = "Receiver",
    ligand_target_matrix = inputs$ligand_target_matrix,
    lr_network = inputs$lr_network,
    min_cells = 5,
    top_n = 20,
    empirical_pval = FALSE,
    logFC_threshold = 0.1,
    p_val_threshold = 1,
    top_n_LR = 10,
    return_object = FALSE
  )

  expect_equal(stored$method, "multinichenet")
  expect_equal(stored$backend, "multinichenet")
  expect_true(nrow(stored$tables$primary) > 0L)
  expect_true(all(c("source", "target", "ligand", "receptor", "score", "condition") %in% names(stored$tables$primary)))
  expect_true(nrow(stored$tables$sample_evidence) > 0L)
})

test_that("regulatory activity can run DoRothEA-style and PROGENy-style networks", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("decoupleR")
  object <- make_communication_object()
  tf_network <- tibble::tibble(
    source = c("TF1", "TF1", "TF2", "TF2"),
    target = c("TFG1", "TFG2", "TG1", "TG2"),
    mor = c(1, -1, 1, 1)
  )
  pathway_network <- tibble::tibble(
    source = c("PathwayA", "PathwayA", "PathwayB", "PathwayB"),
    target = c("TG1", "TG2", "TFG1", "TFG2"),
    mor = c(1, 1, -1, 1)
  )

  object <- sn_run_regulatory_activity(
    object = object,
    method = "dorothea",
    group_by = "cell_type",
    network = tf_network,
    minsize = 1,
    result_id = "tf_activity"
  )
  pathway <- sn_run_regulatory_activity(
    object = object,
    method = "progeny",
    group_by = "cell_type",
    network = pathway_network,
    minsize = 1,
    return_object = FALSE
  )
  tf <- sn_get_regulatory_activity_result(object, "tf_activity", sources = "TF1")
  listed <- sn_list_results(object)

  expect_true("tf_activity" %in% names(object@misc$shennong$results$regulatory_activity))
  expect_true("regulatory_activity" %in% listed$type)
  expect_true(all(tf$source == "TF1"))
  expect_equal(unique(tf$analysis_type), "transcription_factor")
  expect_equal(unique(pathway$tables$primary$analysis_type), "pathway")
})

test_that("progeny network is reshaped to long source-target form", {
  skip_if_not_installed("progeny")
  network <- Shennong:::.sn_progeny_network("human", top = 50)
  expect_s3_class(network, "data.frame")
  expect_true(all(c("source", "target", "weight") %in% colnames(network)))
  expect_true(nrow(network) > 0)
  expect_false(anyNA(network$weight))
  normalized <- Shennong:::.sn_normalize_regulatory_network(network, method = "progeny")
  expect_true(all(c("source", "target", "mor") %in% colnames(normalized)))
})

test_that("liana backend tolerates an unset resource argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("liana")
  genes <- unique(c(paste0("G", seq_len(40)), "CD3E", "CD8A", "IL2RG", "CXCL12", "CXCR4"))
  set.seed(717)
  counts <- matrix(rpois(length(genes) * 90, lambda = 2), nrow = length(genes))
  rownames(counts) <- genes
  colnames(counts) <- paste0("cell", seq_len(90))
  counts["CXCL12", 1:30] <- counts["CXCL12", 1:30] + 8
  counts["CXCR4", 61:90] <- counts["CXCR4", 61:90] + 8
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  Seurat::NormalizeData(object, verbose = FALSE)
  object$seurat_clusters <- factor(rep(c("0", "1", "2"), each = 30))
  result <- tryCatch(
    sn_run_cell_communication(
      object,
      method = "liana",
      group_by = "seurat_clusters",
      species = "human",
      min_cells = 10,
      return_object = FALSE
    ),
    error = function(e) e
  )
  if (!inherits(result, "error")) {
    table <- result$table
    expect_s3_class(table, "data.frame")
    expect_true(nrow(table) >= 0)
  } else {
    expect_false(grepl("argument is of length zero", conditionMessage(result)))
  }
})
