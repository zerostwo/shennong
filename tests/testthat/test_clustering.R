library(testthat)

make_test_object <- function(seed, prefix, n_genes = 200, n_cells = 40) {
  set.seed(seed)
  counts <- matrix(rpois(n_genes * n_cells, lambda = 2), nrow = n_genes, ncol = n_cells)
  counts[1:20, 1:(n_cells / 2)] <- counts[1:20, 1:(n_cells / 2)] + 5
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0(prefix, "_cell", seq_len(n_cells))

  sn_initialize_seurat_object(
    x = counts,
    project = prefix
  )
}

make_zero_layer <- function(object) {
  Matrix::Matrix(
    0,
    nrow = nrow(object),
    ncol = ncol(object),
    sparse = TRUE,
    dimnames = dimnames(SeuratObject::LayerData(object, layer = "counts"))
  )
}

make_rare_test_object <- function() {
  set.seed(88)
  counts <- matrix(rpois(220 * 60, lambda = 2), nrow = 220, ncol = 60)
  counts[1:20, 1:30] <- counts[1:20, 1:30] + 4
  counts[201:220, ] <- 0
  counts[201:210, 55:60] <- 25
  counts[211:220, 57:60] <- 18
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  rownames(counts) <- c(
    paste0("gene", seq_len(200)),
    paste0("raregene", seq_len(20))
  )
  colnames(counts) <- paste0("rare_cell", seq_len(60))

  object <- sn_initialize_seurat_object(
    x = counts,
    project = "rare"
  )
  object$sample <- rep(c("A", "B"), each = 30)
  object$seed_cluster <- c(rep("major", 54), rep("rare", 6))
  object
}

make_fake_python_runner <- function(payload_expr) {
  script_path <- tempfile("fake-python-")
  lines <- c(
    "#!/bin/sh",
    "runner=$(printf \"%s\" \"$1\" | sed \"s/^'//; s/'$//\")",
    "out=$(awk -F\"'\" '/with open\\(/ {print $2}' \"$runner\")",
    "if [ -z \"$out\" ]; then",
    "  echo 'missing output path' >&2",
    "  exit 1",
    "fi",
    payload_expr
  )
  writeLines(lines, script_path)
  Sys.chmod(script_path, mode = "755")
  script_path
}

test_that("sn_run_cluster clusters a single dataset with the standard workflow", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 1, prefix = "single")
  clustered <- sn_run_cluster(
    object = object,
    normalization_method = "seurat",
    nfeatures = 50,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
  expect_true("umap" %in% names(clustered@reductions))
})

test_that("sn_run_cluster can return cluster assignments directly and supports scran batch workflows", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 101, prefix = "cluster-return")
  object$batch <- rep(c("A", "B"), each = ncol(object) / 2)

  clusters <- sn_run_cluster(
    object = object,
    normalization_method = "seurat",
    nfeatures = 50,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    species = "human",
    return_cluster = TRUE,
    verbose = FALSE
  )

  expect_length(clusters, ncol(object))
  expect_true(all(!is.na(clusters)))

  skip_if_not_installed("scran")
  scran_object <- make_test_object(seed = 102, prefix = "cluster-scran", n_cells = 120)
  scran_object$batch <- rep(c("A", "B"), each = ncol(scran_object) / 2)
  integration_called <- FALSE
  clustered <- suppressWarnings(
    with_mocked_bindings(
      sn_run_cluster(
        object = scran_object,
        batch = "batch",
        normalization_method = "scran",
        integration_method = "harmony",
        nfeatures = 50,
        npcs = 10,
        dims = 1:10,
        species = "human",
        verbose = FALSE
      ),
      .sn_run_batch_integration = function(object, method, batch, reduction, features, assay, dims, npcs, theta, group_by_vars, integration_control = list(), verbose = TRUE) {
        integration_called <<- TRUE
        object@misc$integration <- list(
          method = method,
          batch_by = batch,
          reduction = reduction,
          input_features = features
        )
        list(object = object, reduction = reduction)
      },
      .package = "Shennong"
    ),
    classes = "warning"
  )

  expect_s4_class(clustered, "Seurat")
  expect_true(integration_called)
  expect_equal(clustered@misc$integration$method, "harmony")
  expect_equal(clustered@misc$integration$batch_by, "batch")
})

test_that("clustering helper functions normalize dims and select HVGs by group", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 21, prefix = "helpers")
  object$sample <- rep(c("A", "B"), each = ncol(object) / 2)
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  global_hvg <- Shennong:::.sn_select_variable_features(
    object = object,
    nfeatures = 25,
    split_by = NULL,
    verbose = FALSE
  )
  grouped_hvg <- Shennong:::.sn_select_variable_features(
    object = object,
    nfeatures = 25,
    split_by = "sample",
    verbose = FALSE
  )

  expect_equal(Shennong:::.sn_resolve_cluster_dims(dims = NULL, npcs = 8), 1:8)
  expect_equal(Shennong:::.sn_resolve_cluster_dims(dims = c(0, -1, 2, 3, 3), npcs = 10), c(2L, 3L, 3L))
  expect_lte(length(global_hvg$features), 25)
  expect_lte(length(grouped_hvg$features), 25)
  expect_identical(Seurat::VariableFeatures(grouped_hvg$object), grouped_hvg$features)
  expect_error(
    Shennong:::.sn_select_variable_features(
      object = object,
      nfeatures = 25,
      split_by = "missing_group",
      verbose = FALSE
    ),
    "hvg_group_by"
  )
})

test_that("rare-feature helper functions cover gini, local markers, and grouped selection", {
  skip_if_not_installed("Seurat")

  object <- make_rare_test_object()
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  expect_equal(Shennong:::.sn_gini_coefficient(c(0, 0, 0)), 0)
  expect_gt(Shennong:::.sn_gini_coefficient(c(0, 0, 10, 20)), 0)

  empty_gini <- Shennong:::.sn_detect_rare_features_gini(
    expr = SeuratObject::LayerData(object, layer = "counts"),
    nfeatures = 10,
    min_cells = 100,
    max_fraction = 0.01
  )
  expect_equal(nrow(empty_gini), 0)

  group_info <- Shennong:::.sn_resolve_rare_groups(
    object = object,
    group_by = "seed_cluster",
    max_fraction = 0.2,
    max_cells = 10
  )
  expect_true("rare" %in% group_info$rare_groups)
  expect_equal(unname(group_info$group_sizes["rare"]), 6)

  marker_features <- Shennong:::.sn_detect_rare_features_local_markers(
    object = object,
    groups = group_info$groups,
    rare_groups = group_info$rare_groups,
    assay = "RNA",
    layer = "data",
    nfeatures = 10,
    min_cells = 3
  )
  expect_true(any(grepl("^raregene", marker_features)))

  selected <- Shennong:::.sn_select_rare_features(
    object = object,
    base_features = rownames(object),
    method = c("gini", "local_markers"),
    assay = "RNA",
    layer = "data",
    nfeatures = 10,
    group_by = "seed_cluster",
    control = list(
      group_max_fraction = 0.2,
      group_max_cells = 10,
      gene_max_fraction = 0.2,
      min_cells = 3
    ),
    verbose = FALSE
  )

  expect_true(length(selected$features) > 0)
  expect_true(any(selected$metadata$method == "gini"))
  expect_true(any(selected$metadata$method == "local_markers"))
  expect_true(any(grepl("^raregene", selected$features)))
})

test_that("rare-feature grouping helpers cover temporary clustering and controls", {
  skip_if_not_installed("Seurat")

  object <- make_rare_test_object()
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  base_features <- c(rownames(object)[1:40], paste0("raregene", 1:10))
  grouped <- Shennong:::.sn_make_temporary_grouping(
    object = object,
    features = base_features,
    npcs = 8,
    dims = 1:6,
    resolution = 0.3
  )
  auto_group_info <- Shennong:::.sn_resolve_rare_groups(
    object = object,
    features = base_features,
    npcs = 8,
    dims = 1:6,
    resolution = 0.3,
    max_fraction = 0.2,
    max_cells = 10
  )

  expect_length(grouped, ncol(object))
  expect_true(all(grouped %in% auto_group_info$groups))
  expect_length(auto_group_info$groups, ncol(object))
  expect_true(all(names(auto_group_info$group_sizes) %in% unique(auto_group_info$groups)))

  selected <- Shennong:::.sn_select_rare_features(
    object = object,
    base_features = base_features,
    method = "local_markers",
    assay = "RNA",
    layer = "data",
    nfeatures = 10,
    group_by = "seed_cluster",
    control = list(group_max_fraction = 0.2, group_max_cells = 10),
    npcs = 8,
    dims = 1:6,
    resolution = 0.3,
    verbose = FALSE
  )

  expect_true(all(selected$metadata$method == "local_markers"))
  expect_true(any(grepl("^raregene", selected$features)))
  expect_error(
    Shennong:::.sn_select_rare_features(
      object = object,
      base_features = base_features,
      method = "local_hvg",
      verbose = FALSE
    ),
    "'arg' should be one of"
  )
})

test_that("rare-feature selection handles none and missing grouping errors explicitly", {
  skip_if_not_installed("Seurat")

  object <- make_rare_test_object()
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  none_selected <- Shennong:::.sn_select_rare_features(
    object = object,
    base_features = rownames(object),
    method = "none",
    verbose = FALSE
  )
  expect_equal(none_selected$features, character(0))
  expect_s3_class(none_selected$metadata, "data.frame")
  expect_null(none_selected$groups)

  expect_error(
    Shennong:::.sn_resolve_rare_groups(
      object = object,
      group_by = "missing_group"
    ),
    "rare_feature_group_by"
  )
})

test_that("rare-cell helper functions support challenging-groups and mocked backends", {
  skip_if_not_installed("Seurat")

  object <- make_rare_test_object()
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- Seurat::ScaleData(object, verbose = FALSE)
  object <- Seurat::FindVariableFeatures(object, nfeatures = 60, verbose = FALSE)
  object <- Seurat::RunPCA(object, npcs = 10, verbose = FALSE, seed.use = 717)

  challenging <- sn_detect_rare_cells(
    object = object,
    method = "challenging_groups",
    group_by = "seed_cluster",
    reduction = "pca",
    dims = 1:10,
    k = 5
  )
  expect_warning(
    legacy_challenging <- sn_detect_rare_cells(
      object = object,
      method = "challenging_groups",
      group = "seed_cluster",
      reduction = "pca",
      dims = 1:10,
      k = 5
    ),
    "`group` is deprecated"
  )

  expect_s3_class(challenging, "data.frame")
  expect_equal(colnames(legacy_challenging), colnames(challenging))
  expect_true(all(c("cell_id", "method", "rare_score", "rare_cell") %in% colnames(challenging)))
  expect_equal(unique(challenging$method), "challenging_groups")

  testthat::with_mocked_bindings(
    {
      sccad_tbl <- sn_detect_rare_cells(
        object = object,
        method = "sccad",
        layer = "data"
      )
      sca_tbl <- sn_detect_rare_cells(
        object = object,
        method = "sca",
        layer = "data",
        k = 5
      )

      expect_true(any(sccad_tbl$rare_cell))
      expect_true("subcluster" %in% colnames(sccad_tbl))
      expect_equal(unique(sccad_tbl$method), "sccad")
      expect_equal(unique(sca_tbl$method), "sca")
      expect_true(all(is.finite(sca_tbl$rare_score)))
    },
    .sn_run_sccad = function(...) {
      list(
        rare_sets = list(colnames(object)[55:60]),
        scores = 0.95,
        sub_clusters = rep("rare_subcluster", ncol(object))
      )
    },
    .sn_run_sca = function(...) {
      embedding <- cbind(seq_len(ncol(object)), rep(1, ncol(object)))
      rownames(embedding) <- colnames(object)
      embedding
    }
  )
})

test_that("rare-cell helper functions validate required inputs", {
  skip_if_not_installed("Seurat")

  object <- make_rare_test_object()
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  expect_error(
    sn_detect_rare_cells(
      object = object,
      method = "challenging_groups"
    ),
    "group"
  )
})

test_that("python rare-cell backends serialize inputs and parse JSON outputs", {
  skip_if_not_installed("Seurat")

  expr <- Matrix::Matrix(
    matrix(
      c(5, 0, 3,
        2, 4, 1,
        0, 6, 2,
        7, 1, 5),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(
        paste0("gene", 1:4),
        paste0("cell", 1:3)
      )
    ),
    sparse = TRUE
  )

  fake_sccad_script <- tempfile("scCAD-", fileext = ".py")
  writeLines("# placeholder", fake_sccad_script)
  fake_python_sccad <- make_fake_python_runner(
    "printf '%s' '{\"rare_sets\":[[\"cell1\",\"cell2\"]],\"scores\":[0.91],\"sub_clusters\":[\"rare_a\",\"rare_a\",\"major_b\"],\"degs_list\":[[\"gene1\",\"gene2\"]]}' > \"$out\""
  )
  sccad_result <- Shennong:::.sn_run_sccad(
    expr = expr,
    cell_ids = colnames(expr),
    gene_ids = rownames(expr),
    python = fake_python_sccad,
    script = fake_sccad_script,
    normalization = TRUE,
    save_full = TRUE
  )

  expect_equal(as.vector(unlist(sccad_result$rare_sets, use.names = FALSE)), c("cell1", "cell2"))
  expect_equal(sccad_result$sub_clusters, c("rare_a", "rare_a", "major_b"))
  expect_equal(as.numeric(sccad_result$scores), 0.91)

  fake_python_sca <- make_fake_python_runner(
    "printf '%s' '{\"reduction\":[[1.0,0.0],[0.5,0.5],[0.0,1.0]]}' > \"$out\""
  )
  sca_embedding <- Shennong:::.sn_run_sca(
    expr = expr,
    python = fake_python_sca,
    n_comps = 5,
    iters = 2,
    nbhd_size = 3,
    model = "wilcoxon",
    seed = 11
  )

  expect_equal(dim(sca_embedding), c(ncol(expr), 2))
  expect_equal(rownames(sca_embedding), colnames(expr))
  expect_equal(unname(sca_embedding[, 1]), c(1, 0.5, 0))
})

test_that("embedding rarity helpers score neighbors and reject degenerate inputs", {
  embeddings <- rbind(
    cell1 = c(0, 0),
    cell2 = c(1, 0),
    cell3 = c(0, 1)
  )

  knn <- Shennong:::.sn_get_embedding_knn(embeddings = embeddings, k = 2)
  rarity <- Shennong:::.sn_score_embedding_rarity(embeddings = embeddings, k = 2)

  expect_equal(dim(knn$dist), c(3, 2))
  expect_equal(sort(names(rarity)), rownames(embeddings))
  expect_true(all(is.finite(rarity)))
  expect_true(all(rarity >= 0))

  expect_error(
    Shennong:::.sn_get_embedding_knn(embeddings = matrix(c(0, 0), nrow = 1), k = 2),
    "At least two cells are required"
  )
})

test_that("python rare-cell backends surface runner failures clearly", {
  expr <- Matrix::Matrix(
    matrix(
      c(1, 2,
        3, 4),
      nrow = 2,
      dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
    ),
    sparse = TRUE
  )

  failing_python <- tempfile("fake-python-fail-")
  writeLines(c("#!/bin/sh", "echo 'backend failure' >&2", "exit 1"), failing_python)
  Sys.chmod(failing_python, mode = "755")

  fake_sccad_script <- tempfile("scCAD-", fileext = ".py")
  writeLines("# placeholder", fake_sccad_script)

  expect_error(
    Shennong:::.sn_run_sccad(
      expr = expr,
      cell_ids = colnames(expr),
      gene_ids = rownames(expr),
      python = failing_python,
      script = fake_sccad_script
    ),
    "scCAD execution failed"
  )

  expect_error(
    Shennong:::.sn_run_sca(
      expr = expr,
      python = failing_python
    ),
    "SCA execution failed"
  )
})

test_that("sn_run_cluster supports the SCTransform workflow for a single dataset", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glmGamPoi")

  object <- make_test_object(seed = 2, prefix = "sct")
  clusters <- sn_run_cluster(
    object = object,
    normalization_method = "sctransform",
    nfeatures = 50,
    npcs = 10,
    dims = 1:10,
    return_cluster = TRUE,
    verbose = FALSE
  )

  expect_length(clusters, ncol(object))
})

test_that("sn_run_cluster supports SCTransform followed by Harmony integration", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glmGamPoi")
  skip_if_not_installed("harmony")

  object1 <- make_test_object(seed = 21, prefix = "sctbatch1")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 22, prefix = "sctbatch2")
  object2$sample <- "pbmc3k"
  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))

  clustered <- sn_run_cluster(
    object = merged,
    batch_by = "sample",
    normalization_method = "sctransform",
    nfeatures = 50,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("harmony" %in% names(clustered@reductions))
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
})

test_that("sn_run_cluster integrates batches with harmony", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("harmony")

  object1 <- make_test_object(seed = 3, prefix = "batch1")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 4, prefix = "batch2")
  object2$sample <- "pbmc3k"

  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))

  clustered <- sn_run_cluster(
    object = merged,
    batch_by = "sample",
    normalization_method = "seurat",
    nfeatures = 50,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("harmony" %in% names(clustered@reductions))
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
})

test_that("sn_run_cluster dispatches the selected integration backend", {
  skip_if_not_installed("Seurat")

  object1 <- make_test_object(seed = 13, prefix = "dispatcha")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 14, prefix = "dispatchb")
  object2$sample <- "pbmc3k"
  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))

  for (backend in c("coralysis", "seurat_cca", "seurat_rpca")) {
    captured <- NULL
    clustered <- testthat::with_mocked_bindings(
      sn_run_cluster(
        object = merged,
        batch_by = "sample",
        normalization_method = "seurat",
        integration_method = backend,
        integration_control = list(example = backend),
        nfeatures = 50,
        block_genes = NULL,
        npcs = 10,
        dims = 1:10,
        verbose = FALSE
      ),
      .sn_run_batch_integration = function(object,
                                           method,
                                           batch,
                                           reduction,
                                           features,
                                           assay,
                                           dims,
                                           npcs,
                                           theta,
                                           group_by_vars,
                                           integration_control = list(),
                                           verbose = TRUE) {
        captured <<- list(
          method = method,
          batch_by = batch,
          reduction = reduction,
          features = features,
          integration_control = integration_control
        )
        embeddings <- Seurat::Embeddings(object = object[[reduction]])
        new_reduction <- switch(
          method,
          coralysis = "coralysis",
          seurat_cca = "integrated.cca",
          seurat_rpca = "integrated.rpca"
        )
        key <- switch(
          method,
          coralysis = "CORALYSIS_",
          seurat_cca = "INTCCA_",
          seurat_rpca = "INTRPCA_"
        )
        object@reductions[[new_reduction]] <- Seurat::CreateDimReducObject(
          embeddings = embeddings,
          key = key,
          assay = assay
        )
        list(object = object, reduction = new_reduction)
      },
      .package = "Shennong"
    )

    expect_equal(captured$method, backend)
    expect_equal(captured$batch, "sample")
    expect_equal(captured$integration_control, list(example = backend))
    expect_true(switch(
      backend,
      coralysis = "coralysis",
      seurat_cca = "integrated.cca",
      seurat_rpca = "integrated.rpca"
    ) %in% names(clustered@reductions))
    expect_true("seurat_clusters" %in% colnames(clustered[[]]))
  }
})

test_that("sn_run_cluster imports scVI and scANVI pixi backend outputs", {
  skip_if_not_installed("Seurat")

  object1 <- make_test_object(seed = 17, prefix = "scvia")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 18, prefix = "scvib")
  object2$sample <- "pbmc3k"
  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))
  merged$cell_label <- rep(c("T", "B"), length.out = ncol(merged))

  for (backend in c("scvi", "scanvi")) {
    captured <- NULL
    runtime_dir <- tempfile("shennong-runtime-")
    control <- list(
      runtime_dir = runtime_dir,
      pixi_project = file.path(runtime_dir, "pixi", "scvi"),
      n_latent = 6,
      max_epochs = 1,
      write_h5ad = FALSE
    )
    if (identical(backend, "scanvi")) {
      control$labels_key <- "cell_label"
    }

    clustered <- testthat::with_mocked_bindings(
      sn_run_cluster(
        object = merged,
        batch_by = "sample",
        normalization_method = "seurat",
        integration_method = backend,
        integration_control = control,
        nfeatures = 50,
        block_genes = NULL,
        npcs = 10,
        dims = 1:6,
        resolution = 0.3,
        verbose = FALSE
      ),
      .sn_execute_scvi_pixi = function(pixi,
                                       manifest_path,
                                       script,
                                       input_dir,
                                       output_dir,
                                       config_path,
                                       environment = NULL,
                                       pixi_home = NULL,
                                       install_pixi = TRUE,
                                       pixi_version = "latest",
                                       pixi_download_url = NULL,
                                       verbose = TRUE) {
        captured <<- list(
          manifest_path = manifest_path,
          script = script,
          input_dir = input_dir,
          output_dir = output_dir,
          environment = environment,
          pixi_home = pixi_home,
          install_pixi = install_pixi,
          config = jsonlite::read_json(config_path, simplifyVector = TRUE)
        )
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        cells <- utils::read.csv(file.path(input_dir, "cells.csv"), stringsAsFactors = FALSE)$cell_id
        latent <- as.data.frame(matrix(seq_along(cells), nrow = length(cells), ncol = 6))
        rownames(latent) <- cells
        utils::write.csv(latent, file.path(output_dir, "latent.csv"))
        metadata <- data.frame(scvi_qc = seq_along(cells), row.names = cells)
        if (identical(captured$config$method, "scanvi")) {
          metadata$scanvi_prediction <- "predicted"
        }
        utils::write.csv(metadata, file.path(output_dir, "obs.csv"))
        jsonlite::write_json(
          list(output_h5ad = file.path(output_dir, "integrated.h5ad")),
          file.path(output_dir, "manifest.json"),
          auto_unbox = TRUE
        )
      },
      .package = "Shennong"
    )

    expect_s4_class(clustered, "Seurat")
    expect_true(backend %in% names(clustered@reductions))
    expect_true("seurat_clusters" %in% colnames(clustered[[]]))
    expect_equal(clustered@misc$integration$method, backend)
    expect_equal(clustered@misc$integration$batch, "sample")
    expect_true(file.exists(file.path(control$pixi_project, "pixi.toml")))
    expect_equal(captured$config$method, backend)
    expect_equal(captured$config$batch_key, "sample")
    expect_equal(captured$config$n_latent, 6)
    expect_true(captured$environment %in% c("cpu", "gpu"))
    expect_true(grepl("pixi/home$", captured$pixi_home))
    expect_true("scvi_qc" %in% colnames(clustered[[]]))
    if (identical(backend, "scanvi")) {
      expect_true("scanvi_prediction" %in% colnames(clustered[[]]))
    }
  }
})

test_that("sn_run_infercnvpy accepts a Seurat object and imports outputs", {
  skip_if_not_installed("Seurat")

  genes <- c("ACTB", "GAPDH", "TP53", "EGFR", "MYC", "PTEN", "CD3D", "MS4A1")
  counts <- Matrix::Matrix(
    matrix(rpois(length(genes) * 6, lambda = 5), nrow = length(genes)),
    sparse = TRUE
  )
  rownames(counts) <- genes
  colnames(counts) <- paste0("cnv_cell", seq_len(6))
  object <- sn_initialize_seurat_object(x = counts, project = "infercnvpy")
  object$cell_type <- rep(c("normal", "tumor"), each = 3)

  captured <- NULL
  runtime_dir <- tempfile("shennong-infercnvpy-")
  updated <- testthat::with_mocked_bindings(
    sn_run_infercnvpy(
      object = object,
      assay = "RNA",
      layer = "counts",
      species = "human",
      reference_by = "cell_type",
      reference_cat = "normal",
      runtime_dir = runtime_dir,
      window_size = 5,
      cnv_score_group_by = "cell_type",
      step = 1,
      run_umap = FALSE,
      install_pixi = FALSE
    ),
    .sn_execute_infercnvpy_pixi = function(script,
                                           input_dir,
                                           output_dir,
                                           config_path,
                                           ...) {
      captured <<- list(
        script = script,
        input_dir = input_dir,
        output_dir = output_dir,
        config = jsonlite::read_json(config_path, simplifyVector = TRUE),
        dots = list(...)
      )
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      obs <- utils::read.csv(file.path(input_dir, "obs.csv"), stringsAsFactors = FALSE)
      metadata <- data.frame(
        cnv_score = seq_len(nrow(obs)) / 10,
        cnv_leiden = rep(c("0", "1"), length.out = nrow(obs)),
        row.names = obs$cell_id
      )
      utils::write.csv(metadata, file.path(output_dir, "obs.csv"))
      pca <- data.frame(
        PC1 = seq_len(nrow(obs)),
        PC2 = rev(seq_len(nrow(obs))),
        row.names = obs$cell_id
      )
      utils::write.csv(pca, file.path(output_dir, "cnv_pca.csv"))
      jsonlite::write_json(
        list(output_h5ad = file.path(output_dir, "infercnvpy.h5ad")),
        file.path(output_dir, "manifest.json"),
        auto_unbox = TRUE
      )
    },
    .package = "Shennong"
  )

  expect_s4_class(updated, "Seurat")
  expect_true("infercnvpy_cnv_score" %in% colnames(updated[[]]))
  expect_true("infercnvpy_cnv_leiden" %in% colnames(updated[[]]))
  expect_true("infercnvpy_cnv_pca" %in% names(updated@reductions))
  expect_equal(captured$config$reference_key, "cell_type")
  expect_equal(captured$config$reference_cat, "normal")
  expect_equal(captured$config$window_size, 5)
  expect_equal(captured$config$cnv_score_groupby, "cell_type")
  expect_true(file.exists(file.path(captured$input_dir, "matrix.mtx")))
  expect_true(file.exists(file.path(captured$input_dir, "var.csv")))
  expect_true("infercnvpy" %in% names(updated@misc))
  expect_true("infercnvpy" %in% names(updated@misc$infercnvpy))
})

test_that("python method wrappers accept Seurat object inputs", {
  skip_if_not_installed("Seurat")

  genes <- c("ACTB", "GAPDH", "TP53", "EGFR", "MYC", "PTEN")
  counts <- Matrix::Matrix(
    matrix(rpois(length(genes) * 6, lambda = 5), nrow = length(genes)),
    sparse = TRUE
  )
  rownames(counts) <- genes
  colnames(counts) <- paste0("py_cell", seq_len(6))
  object <- sn_initialize_seurat_object(x = counts, project = "python-methods")
  object$cell_type <- rep(c("T", "B"), each = 3)
  object$sample <- rep(c("A", "B"), length.out = ncol(object))
  object$x <- seq_len(ncol(object))
  object$y <- rev(seq_len(ncol(object)))
  signatures <- data.frame(T = runif(length(genes)), B = runif(length(genes)), row.names = genes)

  calls <- list(
    scarches = function() sn_run_scarches(object = object, batch_by = "sample", install_pixi = FALSE),
    scpoli = function() sn_run_scpoli(object = object, batch_by = "sample", label_by = "cell_type", install_pixi = FALSE),
    cellphonedb = function() sn_run_cellphonedb(object = object, group_by = "cell_type", install_pixi = FALSE),
    cell2location = function() sn_run_cell2location(object = object, reference_signatures = signatures, spatial_cols = c("x", "y"), install_pixi = FALSE),
    tangram = function() sn_run_tangram(object = object, reference_object = object, cell_type_by = "cell_type", spatial_cols = c("x", "y"), install_pixi = FALSE),
    squidpy = function() sn_run_squidpy(object = object, spatial_cols = c("x", "y"), cluster_by = "cell_type", install_pixi = FALSE),
    spatialdata = function() sn_run_spatialdata(object = object, spatial_cols = c("x", "y"), install_pixi = FALSE),
    stlearn = function() sn_run_stlearn(object = object, spatial_cols = c("x", "y"), install_pixi = FALSE)
  )

  for (method in names(calls)) {
    captured <- NULL
    updated <- testthat::with_mocked_bindings(
      calls[[method]](),
      .sn_execute_python_object_pixi = function(environment,
                                                script,
                                                input_dir,
                                                output_dir,
                                                config_path,
                                                ...) {
        captured <<- list(
          environment = environment,
          script = script,
          input_dir = input_dir,
          output_dir = output_dir,
          config = jsonlite::read_json(config_path, simplifyVector = TRUE)
        )
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        obs <- utils::read.csv(file.path(input_dir, "query", "obs.csv"), stringsAsFactors = FALSE)
        metadata <- data.frame(method_score = seq_len(nrow(obs)), row.names = obs$cell_id)
        utils::write.csv(metadata, file.path(output_dir, "obs.csv"))
        latent <- data.frame(AXIS1 = seq_len(nrow(obs)), AXIS2 = rev(seq_len(nrow(obs))), row.names = obs$cell_id)
        utils::write.csv(latent, file.path(output_dir, "latent.csv"))
        jsonlite::write_json(
          list(output_h5ad = file.path(output_dir, paste0(method, ".h5ad"))),
          file.path(output_dir, "manifest.json"),
          auto_unbox = TRUE
        )
      },
      .package = "Shennong"
    )
    expect_s4_class(updated, "Seurat")
    expect_true(paste0(method, "_method_score") %in% colnames(updated[[]]))
    expect_true(paste0(method, "_latent") %in% names(updated@reductions))
    expect_true(method %in% names(updated@misc))
    expect_true(file.exists(file.path(captured$input_dir, "query", "matrix.mtx")))
    if (identical(method, "cellphonedb")) {
      expect_equal(captured$config$groupby, "cell_type")
    }
    if (method %in% c("tangram")) {
      expect_true(file.exists(file.path(captured$input_dir, "reference", "matrix.mtx")))
    }
    if (method %in% c("cell2location", "tangram", "squidpy", "spatialdata", "stlearn")) {
      expect_true(file.exists(file.path(captured$input_dir, "query", "spatial.csv")))
    }
  }

  legacy_captured <- NULL
  expect_warning(
    testthat::with_mocked_bindings(
      sn_run_cellphonedb(object = object, groupby = "cell_type", install_pixi = FALSE),
      .sn_execute_python_object_pixi = function(environment,
                                                script,
                                                input_dir,
                                                output_dir,
                                                config_path,
                                                ...) {
        legacy_captured <<- jsonlite::read_json(config_path, simplifyVector = TRUE)
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        obs <- utils::read.csv(file.path(input_dir, "query", "obs.csv"), stringsAsFactors = FALSE)
        metadata <- data.frame(method_score = seq_len(nrow(obs)), row.names = obs$cell_id)
        utils::write.csv(metadata, file.path(output_dir, "obs.csv"))
        jsonlite::write_json(
          list(output_h5ad = file.path(output_dir, "cellphonedb.h5ad")),
          file.path(output_dir, "manifest.json"),
          auto_unbox = TRUE
        )
      },
      .package = "Shennong"
    ),
    "`groupby` is deprecated"
  )
  expect_equal(legacy_captured$groupby, "cell_type")
})

test_that("sn_run_cluster runs Coralysis integration end to end", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Coralysis")
  skip_if_not_installed("SingleCellExperiment")

  object <- make_test_object(seed = 16, prefix = "coralysis", n_genes = 120, n_cells = 30)
  object$sample <- rep(c("A", "B", "C"), each = 10)

  clustered <- suppressWarnings(sn_run_cluster(
    object = object,
    batch_by = "sample",
    normalization_method = "seurat",
    integration_method = "coralysis",
    nfeatures = 40,
    block_genes = NULL,
    npcs = 5,
    dims = 1:5,
    resolution = 0.2,
    integration_control = list(
      icp_args = list(
        k = 2,
        L = 3,
        C = 1,
        threads = 1,
        train.with.bnn = FALSE,
        train.k.nn.prop = NULL,
        build.train.set = FALSE,
        use.cluster.seed = FALSE,
        ari.cutoff = 0.1
      ),
      pca_args = list(
        p = 5,
        pca.method = "stats",
        return.model = TRUE
      )
    ),
    verbose = FALSE
  ))

  expect_s4_class(clustered, "Seurat")
  expect_true("coralysis" %in% names(clustered@reductions))
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
  expect_true(inherits(clustered@misc$coralysis, "SingleCellExperiment"))
  expect_equal(clustered@misc$integration$method, "coralysis")
})

test_that("sn_run_cluster restricts SCTransform integration to Harmony", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 15, prefix = "sctmethod")
  object$sample <- rep(c("A", "B"), each = ncol(object) / 2)

  expect_error(
    sn_run_cluster(
      object = object,
      batch_by = "sample",
      normalization_method = "sctransform",
      integration_method = "coralysis",
      nfeatures = 50,
      npcs = 10,
      dims = 1:10,
      verbose = FALSE
    ),
    "SCTransform integration is currently supported only with `integration_method = \"harmony\"`"
  )
})

test_that("sn_run_cluster can select HVGs by metadata group during integration", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("harmony")

  object1 <- make_test_object(seed = 31, prefix = "batcha")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 32, prefix = "batchb")
  object2$sample <- "pbmc3k"
  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))

  clustered <- sn_run_cluster(
    object = merged,
    batch_by = "sample",
    normalization_method = "seurat",
    hvg_group_by = "sample",
    nfeatures = 50,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("harmony" %in% names(clustered@reductions))
  expect_lte(length(Seurat::VariableFeatures(clustered)), 50)
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
})

test_that("sn_run_cluster defaults hvg_group_by to batch unless overridden", {
  skip_if_not_installed("Seurat")

  object1 <- make_test_object(seed = 61, prefix = "autoa")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 62, prefix = "autob")
  object2$sample <- "pbmc3k"
  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))

  captured <- list()
  local_mocked_bindings(
    .sn_select_variable_features = function(object, nfeatures, split_by = NULL, verbose = TRUE) {
      captured <<- c(captured, list(split_by))
      list(object = object, features = rownames(object)[seq_len(min(nfeatures, nrow(object)))])
    },
    .env = asNamespace("Shennong")
  )

  sn_run_cluster(
    object = merged,
    batch_by = "sample",
    normalization_method = "seurat",
    nfeatures = 30,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_equal(captured[[1]], "sample")

  captured <- list()
  local_mocked_bindings(
    .sn_select_variable_features = function(object, nfeatures, split_by = NULL, verbose = TRUE) {
      captured <<- c(captured, list(split_by))
      list(object = object, features = rownames(object)[seq_len(min(nfeatures, nrow(object)))])
    },
    .env = asNamespace("Shennong")
  )

  sn_run_cluster(
    object = merged,
    batch_by = "sample",
    hvg_group_by = "orig.ident",
    normalization_method = "seurat",
    nfeatures = 30,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_equal(captured[[1]], "orig.ident")
})

test_that("sn_run_cluster handles merged split layers with differing feature sets", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("harmony")

  object1 <- make_test_object(seed = 33, prefix = "mergea")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 34, prefix = "mergeb")
  object2$sample <- "pbmc3k"

  object1 <- object1[1:180, ]
  object2 <- object2[c(1:150, 181:200), ]
  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))

  clustered <- sn_run_cluster(
    object = merged,
    batch_by = "sample",
    hvg_group_by = "sample",
    normalization_method = "seurat",
    nfeatures = 50,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("harmony" %in% names(clustered@reductions))
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
})

test_that("sn_run_cluster errors for missing HVG grouping columns", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 41, prefix = "missing-hvg")

  expect_error(
    sn_run_cluster(
      object = object,
      normalization_method = "seurat",
      hvg_group_by = "missing_column",
      nfeatures = 50,
      block_genes = NULL,
      npcs = 10,
      dims = 1:10,
      verbose = FALSE
    ),
    "metadata column name"
  )
})

test_that("sn_detect_rare_cells detects rare-score outliers with the gini method", {
  skip_if_not_installed("Seurat")

  object <- make_rare_test_object()
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  rare_tbl <- sn_detect_rare_cells(
    object,
    method = "gini",
    nfeatures = 20,
    max_fraction = 0.15
  )

  expect_s3_class(rare_tbl, "data.frame")
  expect_true(all(c("cell_id", "method", "rare_score", "rare_cell") %in% colnames(rare_tbl)))
  expect_true(any(rare_tbl$rare_cell))
  expect_true(any(rare_tbl$cell_id[rare_tbl$rare_cell] %in% colnames(object)[55:60]))
})

test_that("sn_run_cluster can append rare-aware features before PCA", {
  skip_if_not_installed("Seurat")

  object <- make_rare_test_object()
  clustered <- sn_run_cluster(
    object = object,
    normalization_method = "seurat",
    nfeatures = 40,
    rare_feature_method = "gini",
    rare_feature_n = 15,
    rare_feature_control = list(gene_max_fraction = 0.15),
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("rare_feature_selection" %in% names(clustered@misc))
  expect_true(length(clustered@misc$rare_feature_selection$rare_features) > 0)
  expect_true(any(grepl("^raregene", clustered@misc$rare_feature_selection$rare_features)))
  expect_true(length(Seurat::VariableFeatures(clustered)) >= 40)
})

test_that("sn_run_cluster merges user-supplied HVGs into the PCA feature set", {
  skip_if_not_installed("Seurat")

  object <- make_rare_test_object()
  user_features <- c("raregene18", "raregene19", "missing_feature")

  clustered <- sn_run_cluster(
    object = object,
    normalization_method = "seurat",
    nfeatures = 40,
    hvg_features = user_features,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_true(all(c("raregene18", "raregene19") %in% Seurat::VariableFeatures(clustered)))
  expect_true(all(c("raregene18", "raregene19") %in% rownames(clustered@reductions$pca@feature.loadings)))
  expect_equal(clustered@misc$hvg_selection$user_features, c("raregene18", "raregene19"))
  expect_equal(clustered@misc$hvg_selection$missing_user_features, "missing_feature")
})

test_that("sn_run_cluster applies rare_feature_n per selected method", {
  skip_if_not_installed("Seurat")

  object <- make_rare_test_object()
  clustered <- sn_run_cluster(
    object = object,
    normalization_method = "seurat",
    nfeatures = 40,
    rare_feature_method = c("gini", "local_markers"),
    rare_feature_group_by = "seed_cluster",
    rare_feature_n = 5,
    rare_feature_control = list(
      group_max_fraction = 0.2,
      group_max_cells = 10,
      gene_max_fraction = 0.2
    ),
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  rare_store <- clustered@misc$rare_feature_selection
  method_counts <- table(rare_store$rare_feature_table$method)

  expect_equal(rare_store$rare_feature_n, 5)
  expect_equal(as.integer(method_counts["gini"]), 5)
  expect_equal(as.integer(method_counts["local_markers"]), 5)
  expect_lte(length(rare_store$rare_features), 10)
  expect_gte(length(rare_store$rare_features), 5)
  expect_equal(rare_store$selected_rare_feature_n, length(rare_store$rare_features))
  expect_true(any(grepl("^raregene", rare_store$rare_features)))
})

test_that("sn_transfer_labels adds compact Seurat label-transfer metadata", {
  skip_if_not_installed("Seurat")

  reference <- make_test_object(seed = 200, prefix = "reference", n_genes = 80, n_cells = 20)
  query <- make_test_object(seed = 201, prefix = "query", n_genes = 80, n_cells = 12)
  reference$cell_type <- rep(c("T cell", "B cell"), each = 10)
  predictions <- data.frame(
    predicted.id = rep(c("T cell", "B cell"), each = 6),
    prediction.score.max = rep(c(0.92, 0.88), each = 6),
    prediction.score.T.cell = rep(c(0.92, 0.12), each = 6),
    prediction.score.B.cell = rep(c(0.08, 0.88), each = 6),
    row.names = colnames(query),
    check.names = FALSE
  )

  mapped <- testthat::with_mocked_bindings(
    sn_transfer_labels(
      object = query,
      reference = reference,
      label_by = "cell_type",
      prediction_prefix = "pbmc_ref",
      dims = 1:5,
      store_prediction_scores = TRUE,
      verbose = FALSE
    ),
    .sn_find_transfer_anchors_backend = function(...) {
      list(anchor_count = 5)
    },
    .sn_transfer_data_backend = function(...) {
      predictions
    }
  )

  expect_s4_class(mapped, "Seurat")
  expect_true(all(c(
    "pbmc_ref_label",
    "pbmc_ref_score",
    "pbmc_ref_score_T_cell",
    "pbmc_ref_score_B_cell"
  ) %in% colnames(mapped[[]])))
  expect_equal(unname(mapped$pbmc_ref_label), predictions$predicted.id)
  expect_equal(mapped@misc$label_transfer$pbmc_ref$label_col, "cell_type")
})

test_that("sn_transfer_labels supports Coralysis reference mapping", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Coralysis")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")

  reference <- make_test_object(seed = 210, prefix = "coralref", n_genes = 40, n_cells = 10)
  query <- make_test_object(seed = 211, prefix = "coralquery", n_genes = 40, n_cells = 6)
  reference$cell_type <- rep(c("T cell", "B cell"), each = 5)
  reference <- Seurat::NormalizeData(reference, verbose = FALSE)
  query <- Seurat::NormalizeData(query, verbose = FALSE)
  ref_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = SeuratObject::LayerData(reference, assay = "RNA", layer = "data")),
    colData = reference[[]]
  )
  reference@misc$coralysis <- ref_sce

  mapped_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = SeuratObject::LayerData(query, assay = "RNA", layer = "data")),
    colData = data.frame(
      coral_labels = rep(c("T cell", "B cell"), each = 3),
      pruned_coral_labels = rep(c("T cell", "B cell"), each = 3),
      coral_probability = rep(c(0.91, 0.87), each = 3),
      row.names = colnames(query)
    )
  )

  mapped <- testthat::with_mocked_bindings(
    sn_transfer_labels(
      object = query,
      reference = reference,
      label_by = "cell_type",
      method = "coralysis",
      prediction_prefix = "coral_ref",
      transfer_control = list(k.nn = 5),
      verbose = FALSE
    ),
    .sn_coralysis_reference_mapping_backend = function(...) {
      mapped_sce
    },
    .package = "Shennong"
  )

  expect_s4_class(mapped, "Seurat")
  expect_true(all(c("coral_ref_label", "coral_ref_score", "coral_ref_raw_label") %in% colnames(mapped[[]])))
  expect_equal(unname(mapped$coral_ref_label), mapped_sce$pruned_coral_labels)
  expect_equal(mapped@misc$label_transfer$coral_ref$method, "coralysis")
  expect_equal(mapped@misc$label_transfer$coral_ref$transfer_control, list(k.nn = 5))
})

test_that("sn_transfer_labels supports scANVI and scArches-style label transfer", {
  skip_if_not_installed("Seurat")

  reference <- make_test_object(seed = 220, prefix = "scanviref", n_genes = 60, n_cells = 14)
  query <- make_test_object(seed = 221, prefix = "scanviquery", n_genes = 60, n_cells = 8)
  reference$cell_type <- rep(c("T cell", "B cell"), each = 7)

  for (backend in c("scanvi", "scarches")) {
    captured <- NULL
    mapped <- testthat::with_mocked_bindings(
      sn_transfer_labels(
        object = query,
        reference = reference,
        label_by = "cell_type",
        method = backend,
        prediction_prefix = paste0(backend, "_ref"),
        features = rownames(reference)[1:30],
        transfer_control = list(
          runtime_dir = tempfile("shennong-transfer-"),
          pixi_project = tempfile("shennong-transfer-pixi-"),
          max_epochs = 1,
          scanvi_max_epochs = 1,
          write_h5ad = FALSE,
          install_pixi = FALSE
        ),
        verbose = FALSE
      ),
      .sn_run_scvi_integration = function(object,
                                          method,
                                          batch,
                                          features,
                                          assay,
                                          integration_control = list(),
                                          verbose = TRUE) {
        captured <<- list(
          method = method,
          batch = batch,
          features = features,
          labels_key = integration_control$labels_key,
          unlabeled_category = integration_control$unlabeled_category
        )
        transfer_role <- object[[".sn_transfer_role"]][, 1]
        transfer_label <- object[[integration_control$labels_key]][, 1]
        prediction <- ifelse(transfer_role == "query", "T cell", transfer_label)
        object$scanvi_prediction <- prediction
        object@misc$integration <- list(run_dir = tempfile("scanvi-run-"), output_h5ad = tempfile("scanvi-.h5ad"))
        list(object = object, reduction = "scanvi")
      },
      .package = "Shennong"
    )

    expect_s4_class(mapped, "Seurat")
    expect_true(paste0(backend, "_ref_label") %in% colnames(mapped[[]]))
    expect_equal(unique(unname(mapped[[paste0(backend, "_ref_label")]][, 1])), "T cell")
    expect_equal(mapped@misc$label_transfer[[paste0(backend, "_ref")]]$method, backend)
    expect_equal(captured$method, "scanvi")
    expect_equal(captured$labels_key, ".sn_transfer_label")
    expect_equal(captured$unlabeled_category, "Unknown")
  }
})

test_that("sn_simulate dispatches scDesign3 output into Seurat objects", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("scDesign3")

  object <- make_test_object(seed = 202, prefix = "simulate", n_genes = 30, n_cells = 10)
  object$cell_type <- rep(c("A", "B"), each = 5)
  sim_counts <- Matrix::Matrix(
    matrix(rpois(30 * 6, lambda = 2), nrow = 30, ncol = 6),
    sparse = TRUE
  )
  rownames(sim_counts) <- rownames(object)
  colnames(sim_counts) <- paste0("sim_cell_", seq_len(6))

  simulated <- testthat::with_mocked_bindings(
    sn_simulate_scdesign3(
      object = object,
      celltype = "cell_type",
      ncell = 6,
      n_cores = 1,
      return = "seurat",
      project = "mock_scdesign3"
    ),
    .sn_run_scdesign3_backend = function(...) {
      args <- list(...)
      expect_equal(args$celltype, "cell_type")
      expect_equal(args$mu_formula, "cell_type")
      expect_equal(args$corr_formula, "cell_type")
      expect_equal(args$ncell, 6)
      list(
        new_count = sim_counts,
        new_covariate = data.frame(
          cell_type = rep(c("A", "B"), each = 3),
          row.names = colnames(sim_counts)
        )
      )
    }
  )

  expect_s4_class(simulated, "Seurat")
  expect_equal(ncol(simulated), 6)
  expect_equal(rownames(simulated), rownames(object))
  expect_equal(simulated@misc$scdesign3$mu_formula, "cell_type")

  generic <- testthat::with_mocked_bindings(
    sn_simulate(
      object = object,
      method = "scdesign3",
      celltype = "cell_type",
      ncell = 6,
      n_cores = 1,
      return = "counts"
    ),
    .sn_run_scdesign3_backend = function(...) {
      list(
        new_count = sim_counts,
        new_covariate = data.frame(
          cell_type = rep(c("A", "B"), each = 3),
          row.names = colnames(sim_counts)
        )
      )
    }
  )
  expect_s4_class(generic, "dgCMatrix")

  combined <- testthat::with_mocked_bindings(
    sn_simulate(
      object = object,
      method = "scdesign3",
      celltype = "cell_type",
      ncell = 6,
      n_cores = 1,
      combine_original = TRUE
    ),
    .sn_run_scdesign3_backend = function(...) {
      list(
        new_count = sim_counts,
        new_covariate = data.frame(
          cell_type = rep(c("A", "B"), each = 3),
          row.names = colnames(sim_counts)
        )
      )
    }
  )
  expect_s4_class(combined, "Seurat")
  expect_true(all(c("original", "simulated") %in% combined$simulation_source))
})

test_that("sn_run_cluster skips cell-cycle scoring cleanly when markers do not overlap", {
  skip_if_not_installed("Seurat")

  object1 <- make_test_object(seed = 13, prefix = "batchx")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 14, prefix = "batchy")
  object2$sample <- "pbmc3k"
  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))

  expect_no_error(
    clustered <- suppressWarnings(sn_run_cluster(
      object = merged,
      batch_by = "sample",
      species = "human",
      normalization_method = "seurat",
      nfeatures = 50,
      npcs = 10,
      dims = 1:10,
      verbose = FALSE
    ))
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
  expect_false(all(c("S.Score", "G2M.Score", "Phase", "CC.Difference") %in% colnames(clustered[[]])))
})

test_that("sn_run_cluster can use a non-default layer without overwriting counts", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 11, prefix = "layered")
  original_counts <- SeuratObject::LayerData(object, layer = "counts")
  zero_counts <- make_zero_layer(object)

  SeuratObject::LayerData(object, layer = "counts") <- zero_counts
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- original_counts

  clustered <- sn_run_cluster(
    object = object,
    normalization_method = "seurat",
    assay = "RNA",
    layer = "decontaminated_counts",
    nfeatures = 50,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
  expect_equal(
    as.matrix(SeuratObject::LayerData(clustered, layer = "counts")),
    as.matrix(zero_counts)
  )
})

test_that("sn_remove_ambient_contamination requires raw counts for SoupX", {
  counts <- matrix(rpois(100 * 20, lambda = 3), nrow = 100, ncol = 20)
  rownames(counts) <- paste0("gene", seq_len(100))
  colnames(counts) <- paste0("cell", seq_len(20))

  expect_error(
    sn_remove_ambient_contamination(
      x = counts,
      method = "soupx",
      return_object = FALSE
    ),
    "`raw` is required"
  )
})

test_that("sn_remove_ambient_contamination reuses stored raw paths from initialization metadata", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(100 * 20, lambda = 3), nrow = 100, ncol = 20), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(100))
  colnames(counts) <- paste0("cell", seq_len(20))
  raw_counts <- counts + 1
  object <- SeuratObject::CreateSeuratObject(counts = counts, project = "stored-raw")
  raw_dir <- tempfile("raw-feature-")
  dir.create(raw_dir)
  Seurat::Misc(object, "input_source") <- list(
    type = "10x_outs",
    raw_path = raw_dir
  )

  captured <- NULL
  corrected <- with_mocked_bindings(
    sn_remove_ambient_contamination(
      x = object,
      method = "soupx",
      return_object = FALSE,
      verbose = FALSE
    ),
    .sn_resolve_counts_input = function(x, arg = "x") {
      if (inherits(x, "Seurat")) {
        return(list(object = x, counts = counts))
      }
      if (identical(x, raw_dir)) {
        return(list(object = NULL, counts = raw_counts))
      }
      stop("unexpected input")
    },
    .sn_remove_ambient_soupx = function(x_info, raw_info, ...) {
      captured <<- list(
        x_counts = x_info$counts,
        raw_counts = raw_info$counts
      )
      list(
        counts = x_info$counts,
        metadata = NULL,
        zero_cells = character(0),
        removed_cells = character(0)
      )
    },
    .package = "Shennong"
  )

  expect_equal(as.matrix(corrected), as.matrix(counts))
  expect_equal(as.matrix(captured$x_counts), as.matrix(counts))
  expect_equal(as.matrix(captured$raw_counts), as.matrix(raw_counts))
})

test_that("sn_remove_ambient_contamination supports decontX on matrices", {
  skip_if_not_installed("celda")
  skip_if_not_installed("SingleCellExperiment")

  counts <- matrix(rpois(200 * 20, lambda = 3), nrow = 200, ncol = 20)
  rownames(counts) <- paste0("gene", seq_len(200))
  colnames(counts) <- paste0("cell", seq_len(20))
  cluster <- rep(c("a", "b"), each = 10)

  corrected <- sn_remove_ambient_contamination(
    x = counts,
    method = "decontx",
    cluster = cluster,
    return_object = FALSE,
    verbose = FALSE,
    estimateDelta = FALSE,
    maxIter = 5
  )

  expect_equal(dim(corrected), dim(counts))
  expect_equal(rownames(corrected), rownames(counts))
  expect_equal(colnames(corrected), colnames(counts))
  expect_true(all(corrected == round(corrected)))
})

test_that("sn_remove_ambient_contamination defaults to decontX and writes a new layer for Seurat", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("celda")
  skip_if_not_installed("SingleCellExperiment")

  object <- make_test_object(seed = 5, prefix = "decontx", n_genes = 120, n_cells = 20)
  cluster <- rep(c("a", "b"), each = 10)

  updated <- sn_remove_ambient_contamination(
    x = object,
    cluster = cluster,
    verbose = FALSE,
    estimateDelta = FALSE,
    maxIter = 5
  )

  expect_s4_class(updated, "Seurat")
  expect_true(all(c(
    "decontX_contamination",
    "decontX_clusters",
    "nCount_RNA_corrected",
    "nFeature_RNA_corrected"
  ) %in% colnames(updated[[]])))
  corrected <- SeuratObject::LayerData(updated, layer = "decontaminated_counts")
  expect_equal(dim(corrected), dim(SeuratObject::LayerData(object, layer = "counts")))
  expect_true(all(corrected == round(corrected)))
})

test_that("sn_remove_ambient_contamination marks zero-count corrected cells in metadata", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 45, prefix = "ambient-flag", n_genes = 60, n_cells = 6)
  zero_flag <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
  counts <- SeuratObject::LayerData(object, layer = "counts")
  metadata <- data.frame(
    decontX_contamination = rep(0.1, ncol(counts)),
    decontX_clusters = rep(c("a", "b"), each = 3),
    nCount_corrected = Matrix::colSums(counts),
    nFeature_corrected = Matrix::colSums(counts > 0),
    row.names = colnames(counts)
  )
  out <- list(
    counts = counts,
    metadata = metadata,
    zero_cells = colnames(counts)[zero_flag],
    removed_cells = character(0)
  )

  updated <- Shennong:::.sn_apply_ambient_result_to_object(
    object = object,
    out = out,
    assay = "RNA",
    layer = "decontaminated_counts"
  )

  expect_true("decontaminated_counts_zero_count" %in% colnames(updated[[]]))
  expect_identical(
    as.logical(updated$decontaminated_counts_zero_count),
    zero_flag
  )
})

test_that("sn_remove_ambient_contamination supports BPCells-backed Seurat layers", {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("celda")
  skip_if_not_installed("SingleCellExperiment")

  counts <- Matrix::Matrix(matrix(rpois(120 * 20, lambda = 3), nrow = 120), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(120))
  colnames(counts) <- paste0("cell", seq_len(20))

  dir_path <- tempfile("bpc-")
  BPCells::write_matrix_dir(round(counts), dir_path)
  bp_counts <- BPCells::open_matrix_dir(dir_path)
  object <- SeuratObject::CreateSeuratObject(counts = bp_counts, project = "bpcells")
  cluster <- rep(c("a", "b"), each = 10)

  updated <- sn_remove_ambient_contamination(
    x = object,
    method = "decontx",
    cluster = cluster,
    verbose = FALSE,
    estimateDelta = FALSE,
    maxIter = 5
  )

  expect_s4_class(updated, "Seurat")
  expect_true("decontaminated_counts" %in% SeuratObject::Layers(updated[["RNA"]]))
  corrected <- SeuratObject::LayerData(updated, layer = "decontaminated_counts")
  expect_equal(dim(corrected), c(120, 20))
})

test_that("decontX zero-count handling restores original cells by default", {
  original_counts <- Matrix::Matrix(
    matrix(
      c(
        5, 0,
        1, 2,
        0, 3
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(original_counts) <- paste0("gene", 1:3)
  colnames(original_counts) <- paste0("cell", 1:2)
  corrected_counts <- Matrix::Matrix(
    matrix(
      c(
        0.4, 0,
        0.3, 0,
        0.2, 0
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(corrected_counts) <- rownames(original_counts)
  colnames(corrected_counts) <- colnames(original_counts)

  handled <- Shennong:::.sn_handle_zero_count_cells(
    original_counts = original_counts,
    corrected_counts = round(corrected_counts),
    remove_zero_count_cells = FALSE
  )

  expect_equal(as.matrix(handled$counts), as.matrix(original_counts))
  expect_equal(handled$zero_cells, colnames(original_counts))
  expect_length(handled$removed_cells, 0)
})

test_that("decontX zero-count handling can remove affected cells", {
  original_counts <- Matrix::Matrix(
    matrix(
      c(
        5, 0,
        1, 2,
        0, 3
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(original_counts) <- paste0("gene", 1:3)
  colnames(original_counts) <- paste0("cell", 1:2)
  corrected_counts <- Matrix::Matrix(0, nrow = 3, ncol = 2, sparse = TRUE)
  colnames(corrected_counts) <- colnames(original_counts)
  rownames(corrected_counts) <- rownames(original_counts)

  handled <- Shennong:::.sn_handle_zero_count_cells(
    original_counts = original_counts,
    corrected_counts = corrected_counts,
    remove_zero_count_cells = TRUE
  )

  expect_equal(ncol(handled$counts), 0)
  expect_equal(handled$removed_cells, colnames(original_counts))
})

test_that("sn_find_doublets can analyze a non-default layer", {
  skip_if_not_installed("Seurat")
  skip_if_not(suppressWarnings(requireNamespace("scDblFinder", quietly = TRUE)))
  skip_if_not_installed("SingleCellExperiment")

  object <- make_test_object(seed = 12, prefix = "doublets", n_genes = 250, n_cells = 40)
  original_counts <- SeuratObject::LayerData(object, layer = "counts")
  zero_counts <- make_zero_layer(object)

  object$precluster <- rep(c("a", "b"), each = 20)
  SeuratObject::LayerData(object, layer = "counts") <- zero_counts
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- original_counts

  updated <- suppressWarnings(sn_find_doublets(
    object = object,
    clusters = "precluster",
    assay = "RNA",
    layer = "decontaminated_counts",
    ncores = 1
  ))

  expect_true(all(c("scDblFinder.class_corrected", "scDblFinder.score_corrected") %in% colnames(updated[[]])))
  expect_equal(
    as.matrix(SeuratObject::LayerData(updated, layer = "counts")),
    as.matrix(zero_counts)
  )
})

test_that("sn_find_doublets skips zero-count cells in corrected layers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")

  object <- make_test_object(seed = 15, prefix = "doublets-corrected", n_genes = 40, n_cells = 6)
  object$precluster <- rep(c("a", "b"), each = 3)

  corrected_counts <- SeuratObject::LayerData(object, layer = "counts")
  corrected_counts[, 1] <- 0
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- corrected_counts

  local_mocked_bindings(
    .sn_run_scDblFinder = function(sce, ...) {
      n <- ncol(sce)
      colData <- SummarizedExperiment::colData(sce)
      colData$scDblFinder.class <- rep("singlet", n)
      colData$scDblFinder.score <- seq_len(n) / 10
      SummarizedExperiment::colData(sce) <- colData
      sce
    },
    .env = asNamespace("Shennong")
  )

  updated <- sn_find_doublets(
    object = object,
    clusters = "precluster",
    assay = "RNA",
    layer = "decontaminated_counts",
    min_features = 1,
    ncores = 1
  )

  expect_true(all(c("scDblFinder.class_corrected", "scDblFinder.score_corrected") %in% colnames(updated[[]])))
  expect_identical(as.character(updated$scDblFinder.class_corrected[[1]]), "unresolved")
  expect_true(is.na(updated$scDblFinder.score_corrected[[1]]))
  expect_true(all(!is.na(updated$scDblFinder.class_corrected[-1])))
  expect_true(all(!is.na(updated$scDblFinder.score_corrected[-1])))
  expect_identical(
    levels(updated$scDblFinder.class_corrected),
    c("singlet", "doublet", "unresolved")
  )
})

test_that("sn_find_doublets skips low-feature cells before running scDblFinder", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")

  object <- make_test_object(seed = 18, prefix = "doublets-lowfeat", n_genes = 300, n_cells = 6)
  object$precluster <- rep(c("a", "b"), each = 3)

  corrected_counts <- Matrix::Matrix(0, nrow = 300, ncol = 6, sparse = TRUE)
  rownames(corrected_counts) <- rownames(object)
  colnames(corrected_counts) <- colnames(object)
  corrected_counts[1:50, 1] <- 1
  corrected_counts[1:220, 2:6] <- 1
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- corrected_counts

  local_mocked_bindings(
    .sn_run_scDblFinder = function(sce, ...) {
      n <- ncol(sce)
      colData <- SummarizedExperiment::colData(sce)
      colData$scDblFinder.class <- rep("singlet", n)
      colData$scDblFinder.score <- seq_len(n) / 10
      SummarizedExperiment::colData(sce) <- colData
      sce
    },
    .env = asNamespace("Shennong")
  )

  updated <- sn_find_doublets(
    object = object,
    clusters = "precluster",
    assay = "RNA",
    layer = "decontaminated_counts",
    min_features = 200,
    ncores = 1
  )

  expect_identical(as.character(updated$scDblFinder.class_corrected[[1]]), "unresolved")
  expect_true(is.na(updated$scDblFinder.score_corrected[[1]]))
  expect_true(all(!is.na(updated$scDblFinder.class_corrected[-1])))
  expect_identical(
    levels(updated$scDblFinder.class_corrected),
    c("singlet", "doublet", "unresolved")
  )
})

test_that("sn_run_celltypist adds predicted labels back onto the Seurat object", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 52, prefix = "celltypist", n_genes = 80, n_cells = 8)
  object$precluster <- rep(c("c1", "c2"), each = 4)
  original_counts <- SeuratObject::LayerData(object, layer = "counts")

  fake_celltypist <- tempfile("fake-celltypist-")
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "outdir=''",
      "prefix=''",
      "indata=''",
      "majority='false'",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    --indata) indata=$2; shift 2 ;;",
      "    --outdir) outdir=$2; shift 2 ;;",
      "    --prefix) prefix=$2; shift 2 ;;",
      "    --majority-voting) majority='true'; shift ;;",
      "    *) shift ;;",
      "  esac",
      "done",
      "python3 - <<'PY' \"$indata\" \"$outdir\" \"$prefix\" \"$majority\"",
      "import csv, os, sys",
      "indata, outdir, prefix, majority = sys.argv[1:]",
      "with open(indata, newline='') as handle:",
      "    reader = csv.reader(handle)",
      "    header = next(reader)",
      "cells = header[1:]",
      "output_path = os.path.join(outdir, prefix + 'predicted_labels.csv')",
      "with open(output_path, 'w', newline='') as handle:",
      "    writer = csv.writer(handle)",
      "    if majority == 'true':",
      "        writer.writerow(['', 'predicted_labels', 'over_clustering', 'majority_voting'])",
      "        for i, cell in enumerate(cells):",
      "            writer.writerow([cell, 'Tcell' if i % 2 == 0 else 'Bcell', f'cluster_{(i % 2) + 1}', 'T lineage' if i % 2 == 0 else 'B lineage'])",
      "    else:",
      "        writer.writerow(['', 'predicted_labels'])",
      "        for i, cell in enumerate(cells):",
      "            writer.writerow([cell, 'Tcell' if i % 2 == 0 else 'Bcell'])",
      "PY"
    ),
    fake_celltypist
  )
  Sys.chmod(fake_celltypist, mode = "0755")

  outdir <- tempfile("celltypist-out-")
  dir.create(outdir)

  updated <- sn_run_celltypist(
    x = object,
    celltypist = fake_celltypist,
    model = "Immune_All_Low.pkl",
    outdir = outdir,
    over_clustering = "precluster",
    majority_voting = TRUE,
    quiet = TRUE
  )

  expect_true(file.exists(file.path(outdir, "over_clustering.txt")))
  expect_true(all(c(
    "Immune_All_Low_predicted_labels",
    "Immune_All_Low_over_clustering",
    "Immune_All_Low_majority_voting"
  ) %in% colnames(updated[[]])))
  expect_true("sn_run_celltypist" %in% names(updated@commands))
  expect_equal(
    as.matrix(SeuratObject::LayerData(updated, layer = "counts")),
    as.matrix(original_counts)
  )
})

test_that("sn_run_celltypist returns prediction tables for path inputs", {
  input_data <- tempfile(fileext = ".csv")
  utils::write.csv(
    matrix(
      1:6,
      nrow = 2,
      dimnames = list(c("gene1", "gene2"), c("cell1", "cell2", "cell3"))
    ),
    input_data
  )

  fake_celltypist <- tempfile("fake-celltypist-path-")
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "outdir=''",
      "prefix=''",
      "indata=''",
      "while [[ $# -gt 0 ]]; do",
      "  case \"$1\" in",
      "    --indata) indata=$2; shift 2 ;;",
      "    --outdir) outdir=$2; shift 2 ;;",
      "    --prefix) prefix=$2; shift 2 ;;",
      "    *) shift ;;",
      "  esac",
      "done",
      "python3 - <<'PY' \"$indata\" \"$outdir\" \"$prefix\"",
      "import csv, os, sys",
      "indata, outdir, prefix = sys.argv[1:]",
      "with open(indata, newline='') as handle:",
      "    reader = csv.reader(handle)",
      "    header = next(reader)",
      "cells = header[1:]",
      "output_path = os.path.join(outdir, prefix + 'predicted_labels.csv')",
      "with open(output_path, 'w', newline='') as handle:",
      "    writer = csv.writer(handle)",
      "    writer.writerow(['', 'predicted_labels'])",
      "    for i, cell in enumerate(cells):",
      "        writer.writerow([cell, 'Tcell' if i % 2 == 0 else 'Bcell'])",
      "PY"
    ),
    fake_celltypist
  )
  Sys.chmod(fake_celltypist, mode = "0755")

  outdir <- tempfile("celltypist-path-out-")
  dir.create(outdir)

  predicted <- sn_run_celltypist(
    x = input_data,
    celltypist = fake_celltypist,
    model = "Immune_All_Low.pkl",
    outdir = outdir,
    majority_voting = FALSE,
    over_clustering = "obs_cluster",
    quiet = TRUE
  )

  expect_s3_class(predicted, "tbl_df")
  expect_equal(predicted$cell, c("cell1", "cell2", "cell3"))
  expect_true("Immune_All_Low_predicted_labels" %in% colnames(predicted))
})

test_that("sn_run_celltypist errors when expected outputs are missing", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 53, prefix = "celltypist-missing", n_genes = 40, n_cells = 6)
  fake_celltypist <- tempfile("fake-celltypist-missing-")
  writeLines(
    c(
      "#!/usr/bin/env bash",
      "set -euo pipefail",
      "exit 0"
    ),
    fake_celltypist
  )
  Sys.chmod(fake_celltypist, mode = "0755")

  outdir <- tempfile("celltypist-missing-out-")
  dir.create(outdir)

  expect_error(
    sn_run_celltypist(
      x = object,
      celltypist = fake_celltypist,
      model = "Immune_All_Low.pkl",
      outdir = outdir,
      majority_voting = FALSE,
      quiet = TRUE
    ),
    "did not generate expected output files"
  )
})
