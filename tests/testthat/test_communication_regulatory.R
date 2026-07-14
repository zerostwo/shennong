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
    store_name = "manual",
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
  expect_true(all(c("source", "target", "ligand", "receptor", "score", "p_value", "q_value", "rank", "method", "spatial_distance") %in% names(metadata$table)))
  expect_true("cell_communication" %in% listed$type)
  expect_true("cell_communication_results" %in% names(object@misc))
})

test_that("NicheNet backend runs with supplied priors", {
  skip_if_not_installed("Seurat")
  skip_if_not(suppressWarnings(requireNamespace("nichenetr", quietly = TRUE)), "nichenetr is not installed")
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
  expect_true("test_ligand" %in% colnames(stored$table))
  expect_true("LIG1" %in% stored$table$test_ligand)
  expect_equal(stored$backend, "nichenet")
  expect_equal(nrow(stored$tables$sample_evidence), 8L * nrow(stored$table))
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

test_that("communication backends standardize to one comparable schema", {
  liana <- Shennong:::.sn_standardize_communication(tibble::tibble(
    source = "Sender", target = "Receiver", ligand_complex = "LIG1",
    receptor_complex = "REC1", magnitude_rank = 0.1, specificity_rank = 0.2
  ), method = "liana")
  cellchat <- Shennong:::.sn_standardize_communication(tibble::tibble(
    source = "Sender", target = "Receiver", ligand = "LIG1",
    receptor = "REC1", prob = 0.8, pval = 0.01, pathway_name = "PathwayA"
  ), method = "cellchat")
  multinichenet <- Shennong:::.sn_standardize_communication(tibble::tibble(
    sender = "Sender", receiver = "Receiver", ligand = "LIG1",
    receptor = "REC1", prioritization_score = 0.7, group = "Stim"
  ), method = "multinichenet")
  combined <- dplyr::bind_rows(liana, cellchat, multinichenet)
  consensus <- Shennong:::.sn_communication_consensus(combined)
  concordance <- Shennong:::.sn_communication_concordance(combined)

  expect_true(all(c("source", "target", "ligand", "receptor", "score", "p_value", "q_value", "rank", "method", "condition", "sample", "pathway", "target_genes", "evidence_source", "spatial_distance") %in% names(combined)))
  expect_equal(consensus$n_methods, 3L)
  expect_equal(consensus$evidence_source, "cellchat;liana;multinichenet")
  expect_equal(nrow(concordance), 3L)
  expect_true(all(concordance$shared_edges == 1L))
})

test_that("public communication runner builds a real cross-method consensus", {
  skip_if_not_installed("liana")
  skip_if_not_installed("CellChat")
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
  expect_true(nrow(stored$table) > 0L)
  expect_true(all(c("source", "target", "ligand", "receptor", "score", "condition") %in% names(stored$table)))
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
    store_name = "tf_activity"
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

  expect_true("regulatory_activity_results" %in% names(object@misc))
  expect_true("regulatory_activity" %in% listed$type)
  expect_true(all(tf$source == "TF1"))
  expect_equal(unique(tf$analysis_type), "transcription_factor")
  expect_equal(unique(pathway$table$analysis_type), "pathway")
})
