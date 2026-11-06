make_reproducibility_trajectory_object <- function() {
  skip_if_not_installed("SeuratObject")
  Shennong:::.sn_with_seed(29L, {
    n_each <- 30L
    progress <- seq(0, 1, length.out = n_each)
    embedding <- rbind(
      cbind(progress, rnorm(n_each, 0, .025)),
      cbind(1 + progress, progress + rnorm(n_each, 0, .025)),
      cbind(1 + progress, -progress + rnorm(n_each, 0, .025))
    )
    cells <- paste0("cell", seq_len(nrow(embedding)))
    dimnames(embedding) <- list(cells, c("PC_1", "PC_2"))
    means <- rbind(
      1 + 5 * embedding[, 1], 11 - 5 * embedding[, 1],
      2 + 8 * abs(embedding[, 2]), 2 + 5 * (embedding[, 2] > 0),
      2 + 5 * (embedding[, 2] < 0), matrix(3, nrow = 3L, ncol = length(cells))
    )
    counts <- matrix(rpois(length(means), as.numeric(means)), nrow = nrow(means),
                     dimnames = list(paste0("G", seq_len(nrow(means))), cells))
    object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
    object[["pca"]] <- SeuratObject::CreateDimReducObject(
      embeddings = embedding, key = "PC_", assay = "RNA"
    )
    object$seurat_clusters <- rep(c("root", "effector", "memory"), each = n_each)
    object
  })
}

test_that("trajectory seeds control backend draws and restore caller RNG", {
  object <- make_reproducibility_trajectory_object()
  runner <- function(embedding, clusters, ...) {
    pseudotime <- stats::setNames(runif(nrow(embedding)), rownames(embedding))
    list(pseudotime = pseudotime, lineages = list(Lineage1 = unique(clusters)))
  }
  run <- function(ambient, seed = NULL, control_seed = NULL) {
    set.seed(ambient)
    before <- .Random.seed
    result <- sn_run_trajectory(
      object, method = "palantir", reduction = "pca", test_dynamic = FALSE,
      seed = seed, backend_control = list(seed = control_seed, runner = runner),
      return_object = FALSE
    )
    expect_identical(.Random.seed, before)
    result
  }
  first <- run(101L, 0L, 42L)
  second <- run(102L, 0L, 99L)
  expected <- Shennong:::.sn_with_seed(0L, runif(ncol(object)))
  expect_equal(first$tables$cells$primary_pseudotime, expected, tolerance = 1e-14)
  expect_identical(first$tables$cells, second$tables$cells)
  expect_identical(first$provenance$random_seed, 0L)
  control_only <- run(103L, control_seed = 0L)
  expect_identical(control_only$tables$cells, first$tables$cells)
  default <- run(104L)
  explicit_default <- run(105L, seed = 717L)
  expect_identical(default$tables$cells, explicit_default$tables$cells)
  expect_identical(default$provenance$random_seed, 717L)
  expect_false(identical(default$tables$cells, first$tables$cells))

  set.seed(106L)
  before <- .Random.seed
  expect_error(sn_run_trajectory(
    object, method = "palantir", reduction = "pca", test_dynamic = FALSE,
    seed = 5L, backend_control = list(runner = function(...) {
      runif(1)
      stop("backend failure")
    })
  ), "backend failure")
  expect_identical(.Random.seed, before)
  withr::local_preserve_seed()
  rm(".Random.seed", envir = .GlobalEnv)
  sn_run_trajectory(
    object, method = "palantir", reduction = "pca", test_dynamic = FALSE,
    seed = 0L, backend_control = list(runner = runner), return_object = FALSE
  )
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("trajectory rejects invalid seeds before executing a backend", {
  object <- make_reproducibility_trajectory_object()
  runner <- function(...) stop("backend should not run")
  for (seed in list(NA_real_, NaN, Inf, -1, .5, c(1, 2), "1", TRUE, numeric(), 2^31)) {
    for (control in c(FALSE, TRUE)) {
      expect_error(sn_run_trajectory(
        object, method = "palantir", reduction = "pca", test_dynamic = FALSE,
        seed = if (control) NULL else seed,
        backend_control = list(runner = runner, seed = if (control) seed else NULL)
      ), "seed.*non-negative.*integer")
    }
  }
})

test_that("real Slingshot and tradeSeq honor the declared seed across ambient RNG states", {
  skip_if_not_installed("slingshot")
  skip_if_not_installed("tradeSeq")
  object <- make_reproducibility_trajectory_object()
  run <- function(ambient) {
    set.seed(ambient)
    before <- .Random.seed
    result <- sn_run_trajectory(
      object, reduction = "pca", start = "root", end = c("effector", "memory"),
      test_dynamic = TRUE, dynamic_features = paste0("G", 1:8),
      max_dynamic_features = 8L, nknots = 3L, trend_features = c("G1", "G2"),
      trend_points = 10L, seed = 777L, return_object = FALSE, verbose = FALSE
    )
    expect_identical(.Random.seed, before)
    result
  }
  first <- run(101L)
  second <- run(102L)
  expect_equal(first$tables$cells, second$tables$cells, tolerance = 1e-12)
  expect_equal(first$tables$dynamic_genes, second$tables$dynamic_genes, tolerance = 1e-12)
  expect_equal(first$tables$fitted_trends, second$tables$fitted_trends, tolerance = 1e-12)
  expect_identical(first$provenance$random_seed, 777L)
})

test_that("direct Monocle uses effective graph partition and seed settings", {
  skip_if_not_installed("monocle3")
  skip_if_not_installed("Seurat")
  object <- make_reproducibility_trajectory_object()
  object <- suppressWarnings(suppressMessages(Seurat::RunUMAP(
    object, reduction = "pca", dims = 1:2, seed.use = 5L, verbose = FALSE
  )))
  cluster_cells <- monocle3::cluster_cells
  observed_seeds <- integer()
  testthat::local_mocked_bindings(
    cluster_cells = function(...) {
      args <- list(...)
      observed_seeds <<- c(observed_seeds, args$random_seed)
      do.call(cluster_cells, args)
    },
    .package = "monocle3"
  )
  run <- function(ambient, use_partition = FALSE) {
    set.seed(ambient)
    before <- .Random.seed
    on.exit(expect_identical(.Random.seed, before), add = TRUE)
    sn_run_trajectory(
      object, method = "monocle3", start = "root", test_dynamic = FALSE,
      seed = 777L, return_object = FALSE, verbose = FALSE,
      backend_control = list(monocle3 = list(
        use_partition = TRUE,
        cluster = list(k = 10L, cluster_method = "louvain", random_seed = 1L),
        learn_graph = list(use_partition = use_partition)
      ))
    )
  }
  first <- run(101L)
  second <- run(102L)
  expect_equal(first$tables$cells, second$tables$cells, tolerance = 1e-12)
  expect_gt(first$diagnostics$backend$partitions, 1L)
  expect_identical(first$diagnostics$backend$use_partition, FALSE)
  expect_identical(first$provenance$random_seed, 777L)
  expect_identical(observed_seeds, c(777L, 777L))
  expect_error(run(103L, use_partition = TRUE), "multiple disconnected partitions")
})

test_that("Python integration configuration preserves a zero seed", {
  python <- Sys.which("python3")
  skip_if(!nzchar(python), "python3 is required for the configuration probe")
  pixi_root <- system.file("pixi", package = "Shennong")
  script <- tempfile(fileext = ".py")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "import ast, pathlib, sys",
    "root = pathlib.Path(sys.argv[1])",
    "paths = ('scvi/scripts/scvi_integration.py', 'scarches/scripts/scpoli_integration.py', 'bbknn/scripts/bbknn_integration.py')",
    "for path in paths:",
    "    tree = ast.parse((root / path).read_text())",
    "    run = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'run')",
    "    seed_statement = next(n for n in run.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'seed' for t in n.targets))",
    "    code = compile(ast.Module(body=[seed_statement], type_ignores=[]), path, 'exec')",
    "    for config, expected in (({'seed': 0}, 0), ({'seed': 42}, 42), ({'seed': None}, 717), ({}, 717)):",
    "        scope = {'config': config}",
    "        exec(code, scope)",
    "        assert scope['seed'] == expected, (path, config, scope['seed'], expected)",
    "print('seed configuration probes passed')"
  ), script)
  output <- suppressWarnings(system2(python, c(shQuote(script), shQuote(pixi_root)),
                                     stdout = TRUE, stderr = TRUE))
  expect_equal(attr(output, "status") %||% 0L, 0L, info = paste(output, collapse = "\n"))
})
