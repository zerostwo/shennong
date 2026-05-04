make_communication_object <- function() {
  genes <- c("LIG1", "REC1", "TG1", "TG2", "TFG1", "TFG2", "ACTB", "MALAT1")
  counts <- matrix(rpois(length(genes) * 24, lambda = 2), nrow = length(genes), ncol = 24)
  rownames(counts) <- genes
  colnames(counts) <- paste0("cell", seq_len(24))
  cell_type <- rep(c("Sender", "Receiver"), each = 12)
  condition <- rep(rep(c("Ctrl", "Stim"), each = 6), 2)
  counts[rownames(counts) == "LIG1", cell_type == "Sender"] <- counts[rownames(counts) == "LIG1", cell_type == "Sender"] + 5
  counts[rownames(counts) == "REC1", cell_type == "Receiver"] <- counts[rownames(counts) == "REC1", cell_type == "Receiver"] + 5
  counts[rownames(counts) %in% c("TG1", "TG2"), cell_type == "Receiver" & condition == "Stim"] <-
    counts[rownames(counts) %in% c("TG1", "TG2"), cell_type == "Receiver" & condition == "Stim"] + 6
  object <- sn_initialize_seurat_object(
    x = Matrix::Matrix(counts, sparse = TRUE),
    project = "communication",
    species = "human"
  )
  object$cell_type <- cell_type
  object$condition <- condition
  Seurat::NormalizeData(object, verbose = FALSE)
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
    expressed_pct = 0.05,
    top_n = 5,
    return_object = FALSE
  )

  expect_equal(stored$method, "nichenetr")
  expect_true("test_ligand" %in% colnames(stored$table))
  expect_true("LIG1" %in% stored$table$test_ligand)
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
