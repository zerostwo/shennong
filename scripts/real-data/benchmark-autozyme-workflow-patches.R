#!/usr/bin/env Rscript

# Fresh-process, three-arm benchmark for AutoZyme patches reached through
# Shennong workflows. This runner deliberately distinguishes patch eligibility
# and activation from evidence that an input-specific fast path actually ran.

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  file.path(getwd(), "scripts", "real-data", "benchmark-autozyme-workflow-patches.R")
}
script_file <- normalizePath(script_file, winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(
  file.path(dirname(script_file), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

.option <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  inline <- args[startsWith(args, prefix)]
  if (length(inline)) {
    return(sub(prefix, "", inline[[length(inline)]], fixed = TRUE))
  }
  index <- match(paste0("--", name), args)
  if (!is.na(index) && index < length(args)) return(args[[index + 1L]])
  default
}

.flag <- function(name) paste0("--", name) %in% args

.absolute_path <- function(path, must_work = FALSE) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(repo_root, path)
  normalizePath(path, winslash = "/", mustWork = must_work)
}

.integer_option <- function(name, default, minimum = 0L) {
  value <- suppressWarnings(as.integer(.option(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop("`--", name, "` must be one integer >= ", minimum, ".", call. = FALSE)
  }
  value
}

.csv_option <- function(name, default) {
  value <- .option(name, default)
  value <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  unique(value[nzchar(value)])
}

.default_input <- function() {
  from_env <- Sys.getenv("SHENNONG_AUTOZYME_FIXTURE", unset = "")
  candidates <- c(
    from_env,
    file.path(repo_root, "data-local", "pkgdown-real", "single-cell", "kotliarov_pbmc.qs2"),
    file.path(repo_root, "data-local", "single-cell", "capacity_30k.qs2"),
    file.path(repo_root, "data-local", "capacity", "seurat_30k.qs2"),
    file.path(repo_root, "data-local", "capacity", "30k.qs2")
  )
  candidates <- path.expand(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) existing[[1L]] else candidates[[1L]]
}

.operation_specs <- list(
  lisi = list(
    patch = "lisi",
    target = "lisi::compute_lisi",
    strict = FALSE,
    tolerance = 1e-12,
    target_intersection = TRUE,
    intersection_evidence = "sn_calculate_lisi() directly calls lisi::compute_lisi().",
    benchmark_value = "yes-current-installed"
  ),
  ucell = list(
    patch = "ucell",
    target = "UCell::AddModuleScore_UCell; UCell:::calculate_Uscore",
    strict = FALSE,
    tolerance = 1e-12,
    target_intersection = TRUE,
    intersection_evidence = "sn_score_programs(method='ucell') calls AddModuleScore_UCell().",
    benchmark_value = "yes-version-drift"
  ),
  scdblfinder = list(
    patch = "scdblfinder",
    target = "scDblFinder::scDblFinder plus validated internal default-processing targets",
    strict = FALSE,
    tolerance = 1e-8,
    target_intersection = TRUE,
    intersection_evidence = "sn_find_doublets() calls scDblFinder::scDblFinder() on its supported sparse default path.",
    benchmark_value = "yes-supported-envelope-only"
  ),
  seurat_normalize = list(
    patch = "seurat",
    target = "Seurat::NormalizeData.Seurat",
    strict = FALSE,
    tolerance = 1e-8,
    target_intersection = TRUE,
    intersection_evidence = "sn_normalize_data(method='seurat') calls Seurat::NormalizeData().",
    benchmark_value = "yes-operation-specific"
  ),
  merge = list(
    patch = "seurat_merge",
    target = "SeuratObject:::merge.Assay5",
    strict = FALSE,
    tolerance = 0,
    target_intersection = TRUE,
    intersection_evidence = "Shennong's Seurat scope observes merge(); Assay5 dispatch reaches merge.Assay5.",
    benchmark_value = "regression-existing-evidence"
  ),
  joinlayers = list(
    patch = "seurat_joinlayers",
    target = "SeuratObject:::JoinLayers.Assay5",
    strict = FALSE,
    tolerance = 0,
    target_intersection = TRUE,
    intersection_evidence = "Shennong's Seurat scope observes JoinLayers(); Assay5 dispatch reaches JoinLayers.Assay5.",
    benchmark_value = "regression-existing-evidence"
  ),
  clusterprofiler_go_cache = list(
    patch = "clusterprofiler",
    target = "clusterProfiler:::get_GO_data",
    strict = FALSE,
    tolerance = 1e-12,
    target_intersection = TRUE,
    intersection_evidence = "Repeated sn_enrich(database='GOBP') calls reach get_GO_data(); only the repeated identical request can be a cache hit.",
    benchmark_value = "yes-repeated-go-only"
  )
)

default_operations <- names(.operation_specs)

if (.flag("help")) {
  cat(
    "Usage: Rscript scripts/real-data/benchmark-autozyme-workflow-patches.R [options]\n",
    "\n",
    "Fresh-process arms:\n",
    "  direct_upstream       Direct upstream call with AutoZyme disabled\n",
    "  shennong_autozyme_off Shennong workflow with automatic acceleration disabled\n",
    "  shennong_autozyme_on  Shennong workflow with its real lazy AutoZyme scope\n",
    "\n",
    "Options:\n",
    "  --input PATH          Seurat fixture in .qs2 or .rds format.\n",
    "                        Default: SHENNONG_AUTOZYME_FIXTURE, validated Kotliarov,\n",
    "                        then known local 30k-capacity paths.\n",
    "  --output PATH         JSON report (default: INPUT.autozyme-workflow.json)\n",
    "  --operations CSV      Any of: ", paste(default_operations, collapse = ","), "\n",
    "  --repetitions N       Fresh-process repetitions per arm (default: 3)\n",
    "  --repetition-start N   First repetition/order index (default: 1)\n",
    "  --seed N              Deterministic seed (default: 717)\n",
    "  --max-cells N         Deterministic fixture cap; 0 keeps all (default: 30000)\n",
    "  --merge-samples N     Maximum sample/chunk objects for merge tests (default: 4)\n",
    "  --pca-features N      Features used only if benchmark PCA must be built (default: 2000)\n",
    "  --pca-dims N          PCA dimensions used by LISI preparation (default: 20)\n",
    "  --assay NAME          Seurat assay (default: RNA, then DefaultAssay fallback)\n",
    "  --layer NAME          Expression layer for UCell; counts are used where required\n",
    "                        (default: counts)\n",
    "  --reduction NAME      Existing reduction for LISI (default: pca)\n",
    "  --sample-by NAME      Sample metadata; auto-detected when omitted\n",
    "  --label-by NAME       LISI label metadata; auto-detected when omitted\n",
    "  --species NAME        human or mouse for GO cache operation (default: human)\n",
    "  --orgdb PACKAGE       OrgDb package (default: species-specific)\n",
    "  --go-genes N          Query genes for GO cache operation (default: 500)\n",
    "  --allow-mismatch      Write report and exit zero when equivalence fails\n",
    "  --help                Show this message\n",
    "\n",
    "Timing covers the analytical operation. GNU time peak RSS covers the complete\n",
    "fresh worker, including fixture load and untimed preparation. Patch activation\n",
    "is recorded but never reported as proof of an input-specific fast-path hit.\n",
    sep = ""
  )
  quit(status = 0L)
}

worker <- .flag("worker")
input <- .absolute_path(.option("input", .default_input()), must_work = TRUE)
operations <- .csv_option("operations", paste(default_operations, collapse = ","))
unknown_operations <- setdiff(operations, default_operations)
if (length(unknown_operations)) {
  stop("Unknown operation(s): ", paste(unknown_operations, collapse = ", "), ".", call. = FALSE)
}
repetitions <- .integer_option("repetitions", 3L, minimum = 1L)
repetition_start <- .integer_option("repetition-start", 1L, minimum = 1L)
seed <- .integer_option("seed", 717L, minimum = 0L)
max_cells <- .integer_option("max-cells", 30000L, minimum = 0L)
merge_samples <- .integer_option("merge-samples", 4L, minimum = 2L)
pca_features <- .integer_option("pca-features", 2000L, minimum = 20L)
pca_dims <- .integer_option("pca-dims", 20L, minimum = 2L)
go_genes <- .integer_option("go-genes", 500L, minimum = 10L)
assay_option <- .option("assay", "RNA")
layer_option <- .option("layer", "counts")
reduction_option <- .option("reduction", "pca")
sample_by_option <- .option("sample-by", "")
label_by_option <- .option("label-by", "")
species <- tolower(.option("species", "human"))
if (!species %in% c("human", "mouse")) stop("`--species` must be human or mouse.", call. = FALSE)
orgdb <- .option("orgdb", if (species == "human") "org.Hs.eg.db" else "org.Mm.eg.db")
output <- .absolute_path(.option(
  "output",
  paste0(sub("[.](qs2|rds)$", "", input, ignore.case = TRUE), ".autozyme-workflow.json")
))

.require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing benchmark package(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
}

.load_fixture <- function(path) {
  extension <- tolower(tools::file_ext(path))
  value <- switch(
    extension,
    qs2 = {
      .require_packages("qs2")
      qs2::qs_read(path)
    },
    rds = readRDS(path),
    stop("Unsupported fixture extension: .", extension, ". Use .qs2 or .rds.", call. = FALSE)
  )
  if (!inherits(value, "Seurat") && is.list(value) && inherits(value$object, "Seurat")) {
    value <- value$object
  }
  if (!inherits(value, "Seurat")) stop("`--input` must contain a Seurat object.", call. = FALSE)
  value
}

.resolve_metadata_column <- function(object, requested, candidates, require_multiple = FALSE) {
  metadata <- object[[]]
  if (nzchar(requested)) {
    if (!requested %in% colnames(metadata)) {
      stop("Requested metadata column `", requested, "` was not found.", call. = FALSE)
    }
    return(requested)
  }
  hits <- intersect(candidates, colnames(metadata))
  if (require_multiple) {
    hits <- hits[vapply(hits, function(column) {
      length(unique(as.character(metadata[[column]][!is.na(metadata[[column]])]))) > 1L
    }, logical(1))]
  }
  if (length(hits)) hits[[1L]] else ""
}

.subset_fixture <- function(object, maximum, seed) {
  if (maximum == 0L || ncol(object) <= maximum) return(object)
  set.seed(seed)
  cells <- sort(sample(colnames(object), maximum, replace = FALSE))
  object[, cells, drop = FALSE]
}

.ensure_lisi_input <- function(object, assay, reduction, dimensions, feature_cap, label_by) {
  if (!label_by %in% colnames(object[[]])) {
    object$.autozyme_benchmark_label <- rep(c("A", "B"), length.out = ncol(object))
    label_by <- ".autozyme_benchmark_label"
  }
  if (!reduction %in% names(object@reductions)) {
    old_options <- options(shennong.autozyme = FALSE)
    on.exit(options(old_options), add = TRUE)
    object <- autozyme::with_disabled({
      normalized <- Seurat::NormalizeData(object, assay = assay, verbose = FALSE)
      normalized <- Seurat::FindVariableFeatures(
        normalized,
        assay = assay,
        nfeatures = min(feature_cap, nrow(normalized)),
        verbose = FALSE
      )
      features <- SeuratObject::VariableFeatures(normalized[[assay]])
      normalized <- Seurat::ScaleData(normalized, assay = assay, features = features, verbose = FALSE)
      Seurat::RunPCA(
        normalized,
        assay = assay,
        features = features,
        npcs = min(dimensions, length(features) - 1L, ncol(normalized) - 1L),
        reduction.name = reduction,
        verbose = FALSE
      )
    })
  }
  list(object = object, label_by = label_by)
}

.prepare_sample_objects <- function(object, assay, sample_by, maximum_samples) {
  groups <- if (nzchar(sample_by)) as.character(object[[sample_by, drop = TRUE]]) else rep("sample", ncol(object))
  names(groups) <- colnames(object)
  valid_groups <- sort(unique(groups[!is.na(groups) & nzchar(groups)]))
  if (length(valid_groups) < 2L) {
    chunk <- cut(seq_len(ncol(object)), breaks = min(maximum_samples, max(2L, ncol(object))), labels = FALSE)
    groups <- paste0("chunk", chunk)
    names(groups) <- colnames(object)
    valid_groups <- sort(unique(groups))
  }
  valid_groups <- utils::head(valid_groups, maximum_samples)
  objects <- lapply(valid_groups, function(group) {
    cells <- names(groups)[groups == group]
    counts <- SeuratObject::LayerData(object, assay = assay, layer = "counts")[, cells, drop = FALSE]
    result <- SeuratObject::CreateSeuratObject(
      counts = counts,
      assay = assay,
      project = group,
      min.cells = 0,
      min.features = 0
    )
    result$.autozyme_benchmark_sample <- group
    result
  })
  names(objects) <- valid_groups
  objects
}

.prepare_context <- function(operation) {
  object <- .subset_fixture(.load_fixture(input), max_cells, seed)
  assay <- if (assay_option %in% names(object@assays)) assay_option else SeuratObject::DefaultAssay(object)
  layers <- SeuratObject::Layers(object[[assay]])
  layer <- if (layer_option %in% layers) layer_option else if ("counts" %in% layers) "counts" else layers[[1L]]
  if (!"counts" %in% layers && operation %in% c("scdblfinder", "seurat_normalize", "merge", "joinlayers")) {
    stop("Operation `", operation, "` requires a counts layer in assay `", assay, "`.", call. = FALSE)
  }
  object <- object[, colnames(object), drop = FALSE]
  sample_by <- .resolve_metadata_column(
    object,
    sample_by_option,
    c("real_sample", "sample", "sample_id", "donor", "orig.ident"),
    require_multiple = TRUE
  )
  label_by <- .resolve_metadata_column(
    object,
    label_by_option,
    c(sample_by, "batch", "real_response", "condition", "cell_type", "seurat_clusters", "orig.ident"),
    require_multiple = TRUE
  )
  if (operation == "lisi") {
    prepared <- .ensure_lisi_input(object, assay, reduction_option, pca_dims, pca_features, label_by)
    object <- prepared$object
    label_by <- prepared$label_by
  }
  counts <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  if (operation == "scdblfinder") {
    detected <- Matrix::colSums(counts > 0)
    totals <- Matrix::colSums(counts)
    keep <- colnames(counts)[totals > 0 & detected >= 0L]
    object <- object[, keep, drop = FALSE]
    counts <- counts[, keep, drop = FALSE]
  }
  signatures <- NULL
  if (operation == "ucell") {
    features <- rownames(SeuratObject::LayerData(object, assay = assay, layer = layer))
    if (length(features) < 10L) {
      stop("UCell benchmark requires at least ten expression features.", call. = FALSE)
    }
    width <- max(5L, min(25L, floor(length(features) / 4L)))
    signatures <- list(
      benchmark_program_1 = features[seq_len(width)],
      benchmark_program_2 = features[width + seq_len(width)]
    )
  }
  sample_objects <- if (operation %in% c("merge", "joinlayers")) {
    .prepare_sample_objects(object, assay, sample_by, merge_samples)
  } else {
    NULL
  }
  go_query <- NULL
  if (operation == "clusterprofiler_go_cache") {
    totals <- Matrix::rowSums(counts)
    go_query <- names(utils::head(sort(totals[is.finite(totals) & totals > 0], decreasing = TRUE), go_genes))
    if (length(go_query) < 10L) stop("GO cache benchmark requires at least ten nonzero features.", call. = FALSE)
  }
  operation_cells <- ncol(object)
  if (!is.null(sample_objects)) {
    operation_cells <- as.integer(sum(vapply(sample_objects, ncol, numeric(1))))
    samples <- length(sample_objects)
  } else if (nzchar(sample_by)) {
    sample_values <- as.character(object[[sample_by, drop = TRUE]])
    samples <- length(unique(sample_values[!is.na(sample_values) & nzchar(sample_values)]))
  } else {
    samples <- 1L
  }
  list(
    object = object,
    assay = assay,
    layer = layer,
    sample_by = sample_by,
    label_by = label_by,
    signatures = signatures,
    sample_objects = sample_objects,
    go_query = go_query,
    cells = operation_cells,
    features = nrow(object),
    samples = samples
  )
}

.merge_call <- function(sample_objects) {
  base::merge(
    x = sample_objects[[1L]],
    y = sample_objects[-1L],
    add.cell.ids = names(sample_objects),
    merge.data = FALSE
  )
}

.joinlayers_call <- function(sample_objects) {
  counts <- lapply(sample_objects, function(object) {
    SeuratObject::LayerData(object, assay = SeuratObject::DefaultAssay(object), layer = "counts")
  })
  assay <- SeuratObject::CreateAssay5Object(counts = counts)
  SeuratObject::JoinLayers(assay, layers = "counts", new = "counts")
}

.direct_ucell <- function(context) {
  matrix <- SeuratObject::LayerData(context$object, assay = context$assay, layer = context$layer)
  object <- SeuratObject::CreateSeuratObject(
    counts = methods::as(matrix, "CsparseMatrix"),
    min.cells = 0,
    min.features = 0
  )
  scored <- UCell::AddModuleScore_UCell(
    obj = object,
    features = context$signatures,
    name = "",
    ncores = 1L,
    BPPARAM = NULL,
    storeRanks = FALSE,
    slot = "counts"
  )
  metadata <- scored[[]]
  scores <- t(as.matrix(metadata[, names(context$signatures), drop = FALSE]))
  rownames(scores) <- names(context$signatures)
  colnames(scores) <- rownames(metadata)
  scores
}

.direct_scdblfinder <- function(context) {
  counts <- SeuratObject::LayerData(context$object, assay = context$assay, layer = "counts")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = context$object[[]]
  )
  scDblFinder::scDblFinder(
    sce = sce,
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  )
}

.direct_lisi <- function(context) {
  embedding <- Seurat::Embeddings(context$object, reduction = reduction_option)
  dimensions <- seq_len(min(pca_dims, ncol(embedding)))
  metadata <- context$object[[]][rownames(embedding), context$label_by, drop = FALSE]
  result <- lisi::compute_lisi(
    X = embedding[, dimensions, drop = FALSE],
    meta_data = metadata,
    label_colnames = context$label_by
  )
  data.frame(cell_id = rownames(result), result, row.names = NULL, check.names = FALSE)
}

.direct_go_cache <- function(context) {
  call <- function() clusterProfiler::enrichGO(
    gene = context$go_query,
    OrgDb = orgdb,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
  )
  invisible(call())
  result <- call()
  if (is.null(result)) {
    stop(
      "GO benchmark query produced no mapped enrichment result; use a symbol-based fixture or adjust `--go-genes`.",
      call. = FALSE
    )
  }
  result
}

.shennong_go_cache <- function(context) {
  call <- function() Shennong::sn_enrich(
    x = context$go_query,
    analysis = "ora",
    species = species,
    database = "GOBP",
    pvalue_cutoff = 1,
    p_adjust_method = "BH",
    qvalue_cutoff = 1,
    min_gs_size = 10,
    max_gs_size = 500,
    return_object = FALSE
  )
  invisible(call())
  call()
}

.run_direct <- function(operation, context) {
  switch(
    operation,
    lisi = .direct_lisi(context),
    ucell = .direct_ucell(context),
    scdblfinder = .direct_scdblfinder(context),
    seurat_normalize = Seurat::NormalizeData(
      context$object,
      assay = context$assay,
      verbose = FALSE
    ),
    merge = .merge_call(context$sample_objects),
    joinlayers = .joinlayers_call(context$sample_objects),
    clusterprofiler_go_cache = .direct_go_cache(context),
    stop("Unsupported direct operation: ", operation, call. = FALSE)
  )
}

.run_shennong <- function(operation, context) {
  switch(
    operation,
    lisi = Shennong::sn_calculate_lisi(
      context$object,
      reduction = reduction_option,
      label_by = context$label_by,
      dims = seq_len(min(pca_dims, ncol(Seurat::Embeddings(context$object, reduction_option)))),
      seed = seed
    ),
    ucell = Shennong::sn_score_programs(
      context$object,
      signatures = context$signatures,
      method = "ucell",
      assay = context$assay,
      layer = context$layer,
      name = "autozyme_benchmark_ucell",
      backend_control = list(ucell = list(
        name = "",
        ncores = 1L,
        storeRanks = FALSE,
        slot = "counts"
      )),
      return_object = FALSE
    ),
    scdblfinder = Shennong::sn_find_doublets(
      context$object,
      clusters = NULL,
      cluster_backend = "native",
      group_by = NULL,
      dbr_sd = NULL,
      ncores = 1L,
      assay = context$assay,
      layer = "counts",
      min_features = 0L
    ),
    seurat_normalize = Shennong::sn_normalize_data(
      context$object,
      method = "seurat",
      assay = context$assay,
      layer = "counts",
      verbose = FALSE
    ),
    merge = Shennong:::.sn_with_default_autozyme(
      .merge_call(context$sample_objects),
      patches = "seurat_merge",
      strict = FALSE,
      operation = "merge"
    ),
    joinlayers = Shennong:::.sn_with_default_autozyme(
      .joinlayers_call(context$sample_objects),
      patches = "seurat_joinlayers",
      strict = FALSE,
      operation = "joinlayers"
    ),
    clusterprofiler_go_cache = .shennong_go_cache(context),
    stop("Unsupported Shennong operation: ", operation, call. = FALSE)
  )
}

.assay_snapshot <- function(assay) {
  layers <- SeuratObject::Layers(assay)
  layer_payload <- lapply(layers, function(layer) {
    matrix <- SeuratObject::LayerData(assay, layer = layer)
    list(
      class = class(matrix),
      dim = dim(matrix),
      dimnames = dimnames(matrix),
      digest = digest::digest(matrix, algo = "sha256", serialize = TRUE)
    )
  })
  names(layer_payload) <- layers
  list(
    class = class(assay),
    dim = dim(assay),
    dimnames = dimnames(assay),
    layers = layer_payload,
    valid = isTRUE(methods::validObject(assay, test = TRUE))
  )
}

.canonical_table <- function(table, key_columns = character()) {
  table <- as.data.frame(table, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(table) && length(key_columns) && all(key_columns %in% names(table))) {
    table <- table[do.call(order, table[key_columns]), , drop = FALSE]
  }
  rownames(table) <- NULL
  table
}

.canonicalize <- function(operation, value, context, arm) {
  if (operation == "lisi") {
    return(list(kind = "table", value = .canonical_table(value, "cell_id")))
  }
  if (operation == "ucell") {
    table <- if (identical(arm, "direct_upstream")) {
      dplyr::bind_rows(lapply(seq_len(nrow(value)), function(index) {
        data.frame(
          entity = colnames(value),
          program = rownames(value)[[index]],
          score = as.numeric(value[index, ]),
          stringsAsFactors = FALSE
        )
      }))
    } else {
      value$tables$scores[, c("entity", "program", "score"), drop = FALSE]
    }
    return(list(kind = "table", value = .canonical_table(table, c("entity", "program"))))
  }
  if (operation == "scdblfinder") {
    table <- if (identical(arm, "direct_upstream")) {
      data.frame(
        cell = colnames(value),
        class = as.character(value$scDblFinder.class),
        score = as.numeric(value$scDblFinder.score),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        cell = colnames(value),
        class = as.character(value$scDblFinder.class),
        score = as.numeric(value$scDblFinder.score),
        stringsAsFactors = FALSE
      )
    }
    return(list(kind = "table", value = .canonical_table(table, "cell")))
  }
  if (operation == "seurat_normalize") {
    matrix <- SeuratObject::LayerData(value, assay = context$assay, layer = "data")
    return(list(kind = "matrix", value = methods::as(matrix, "CsparseMatrix")))
  }
  if (operation == "merge") {
    return(list(kind = "exact", value = list(
      dim = dim(value),
      cells = colnames(value),
      features = rownames(value),
      metadata = digest::digest(value[[]], algo = "sha256", serialize = TRUE),
      assay = .assay_snapshot(value[[context$assay]])
    )))
  }
  if (operation == "joinlayers") {
    return(list(kind = "exact", value = .assay_snapshot(value)))
  }
  if (operation == "clusterprofiler_go_cache") {
    table <- if (is.null(value)) data.frame() else as.data.frame(value)
    if ("geneID" %in% names(table)) {
      table$geneID <- vapply(strsplit(as.character(table$geneID), "/", fixed = TRUE), function(genes) {
        paste(sort(unique(genes)), collapse = "/")
      }, character(1))
    }
    return(list(kind = "table", value = .canonical_table(table, intersect(c("ID", "Description"), names(table)))))
  }
  stop("No canonicalizer for operation: ", operation, call. = FALSE)
}

.compare_tables <- function(reference, candidate, tolerance) {
  if (!identical(names(reference), names(candidate)) || nrow(reference) != nrow(candidate)) {
    return(list(pass = FALSE, max_abs_diff = Inf, reason = "table columns or row count differ"))
  }
  numeric_columns <- names(reference)[vapply(reference, is.numeric, logical(1))]
  other_columns <- setdiff(names(reference), numeric_columns)
  if (!identical(reference[other_columns], candidate[other_columns])) {
    return(list(pass = FALSE, max_abs_diff = Inf, reason = "non-numeric table values differ"))
  }
  if (!length(numeric_columns)) return(list(pass = TRUE, max_abs_diff = 0, reason = "exact non-numeric table"))
  a <- unlist(reference[numeric_columns], use.names = FALSE)
  b <- unlist(candidate[numeric_columns], use.names = FALSE)
  same_na <- identical(is.na(a), is.na(b))
  finite <- is.finite(a) & is.finite(b)
  max_diff <- if (any(finite)) max(abs(a[finite] - b[finite])) else 0
  pass <- same_na && is.finite(max_diff) && max_diff <= tolerance
  list(pass = pass, max_abs_diff = max_diff, reason = if (pass) "within tolerance" else "numeric drift exceeds tolerance")
}

.compare_canonical <- function(reference, candidate, tolerance) {
  if (!identical(reference$kind, candidate$kind)) {
    return(list(pass = FALSE, max_abs_diff = Inf, reason = "canonical kinds differ"))
  }
  if (reference$kind == "exact") {
    pass <- identical(reference$value, candidate$value)
    return(list(pass = pass, max_abs_diff = if (pass) 0 else Inf, reason = if (pass) "exact" else "exact snapshot differs"))
  }
  if (reference$kind == "table") {
    return(.compare_tables(reference$value, candidate$value, tolerance))
  }
  if (reference$kind == "matrix") {
    a <- methods::as(reference$value, "dgCMatrix")
    b <- methods::as(candidate$value, "dgCMatrix")
    structure_equal <- identical(methods::slot(a, "Dim"), methods::slot(b, "Dim")) &&
      identical(methods::slot(a, "Dimnames"), methods::slot(b, "Dimnames")) &&
      identical(a@i, b@i) && identical(a@p, b@p)
    if (!structure_equal) {
      return(list(pass = FALSE, max_abs_diff = Inf, reason = "sparse matrix structure differs"))
    }
    max_diff <- if (length(a@x)) max(abs(a@x - b@x)) else 0
    pass <- is.finite(max_diff) && max_diff <= tolerance
    return(list(pass = pass, max_abs_diff = max_diff, reason = if (pass) "within tolerance" else "matrix drift exceeds tolerance"))
  }
  list(pass = FALSE, max_abs_diff = Inf, reason = "unknown canonical kind")
}

.package_versions <- function(packages) {
  installed <- packages[vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  stats::setNames(
    lapply(installed, function(package) as.character(utils::packageVersion(package))),
    installed
  )
}

.active_patches <- function() {
  status <- tryCatch(autozyme::status(), error = function(error) character())
  names(status)[status == "active"]
}

.run_worker <- function() {
  operation <- .option("operation", "")
  arm <- .option("arm", "")
  result_path <- .absolute_path(.option("result", ""))
  if (!operation %in% names(.operation_specs)) stop("Worker received an invalid operation.", call. = FALSE)
  if (!arm %in% c("direct_upstream", "shennong_autozyme_off", "shennong_autozyme_on")) {
    stop("Worker received an invalid arm.", call. = FALSE)
  }
  spec <- .operation_specs[[operation]]
  required <- c("Shennong", "autozyme", "digest", "jsonlite", "Matrix", "Seurat", "SeuratObject")
  required <- c(required, switch(
    operation,
    lisi = "lisi",
    ucell = "UCell",
    scdblfinder = c("scDblFinder", "SingleCellExperiment", "BiocParallel"),
    clusterprofiler_go_cache = c("clusterProfiler", orgdb),
    character()
  ))
  .require_packages(unique(required))
  suppressPackageStartupMessages(library(Shennong))

  old_options <- options(shennong.autozyme = identical(arm, "shennong_autozyme_on"))
  on.exit(options(old_options), add = TRUE)
  try(Shennong::sn_disable_autozyme(spec$patch), silent = TRUE)
  on.exit(try(Shennong::sn_disable_autozyme(spec$patch), silent = TRUE), add = TRUE)

  check <- Shennong::sn_check_autozyme(
    patches = spec$patch,
    strict = isTRUE(spec$strict),
    allow_approximate = FALSE
  )
  active_before <- spec$patch %in% .active_patches()
  context <- .prepare_context(operation)
  set.seed(seed)
  invisible(gc())

  value <- NULL
  acceleration <- list()
  timing <- system.time({
    if (identical(arm, "direct_upstream")) {
      value <- autozyme::with_disabled(.run_direct(operation, context))
    } else if (identical(arm, "shennong_autozyme_off")) {
      value <- autozyme::with_disabled(.run_shennong(operation, context))
    } else {
      observed <- Shennong:::.sn_with_autozyme_provenance_context({
        result <- .run_shennong(operation, context)
        list(
          value = result,
          acceleration = Shennong:::.sn_autozyme_provenance()
        )
      }, patches = spec$patch)
      value <- observed$value
      acceleration <- observed$acceleration
    }
  })
  canonical <- .canonicalize(operation, value, context, arm)
  canonical_path <- file.path(dirname(result_path), "canonical.rds")
  saveRDS(canonical, canonical_path, compress = FALSE)
  fingerprint <- digest::digest(canonical, algo = "sha256", serialize = TRUE)
  active_observed <- spec$patch %in% (acceleration$active_patches %||% character())
  try(Shennong::sn_disable_autozyme(spec$patch), silent = TRUE)
  active_after <- spec$patch %in% .active_patches()

  description <- utils::packageDescription("autozyme")
  packages <- unique(c(
    "Shennong", "autozyme", "R", "Seurat", "SeuratObject", "Matrix",
    check$upstream,
    switch(
      operation,
      lisi = "lisi",
      ucell = "UCell",
      scdblfinder = c("scDblFinder", "SingleCellExperiment", "BiocParallel"),
      clusterprofiler_go_cache = c("clusterProfiler", orgdb, "enrichit"),
      character()
    )
  ))
  packages <- setdiff(packages, "R")
  run <- list(
    operation = operation,
    arm = arm,
    patch = spec$patch,
    target = spec$target,
    target_intersection = isTRUE(spec$target_intersection),
    intersection_evidence = spec$intersection_evidence,
    activation_is_fast_path_hit = FALSE,
    fast_path_hit = NA,
    fast_path_hit_evidence = paste(
      "Not instrumented. Eligibility, activation provenance, and a static target intersection",
      "do not prove that the patch's input guard accepted this call."
    ),
    eligible = isTRUE(check$eligible[[1L]]),
    eligibility_reason = as.character(check$reason[[1L]]),
    patch_registered = isTRUE(check$registered[[1L]]),
    patch_provider = as.character(check$patch_provider[[1L]] %||% NA_character_),
    patch_active_before = active_before,
    patch_active_observed = active_observed,
    patch_active_after = active_after,
    activation_evidence = if (active_observed) "Shennong provenance context recorded the patch as used" else "no activation observed",
    elapsed_seconds = unname(timing[["elapsed"]]),
    output_fingerprint = fingerprint,
    canonical_path = normalizePath(canonical_path, winslash = "/", mustWork = TRUE),
    cells = context$cells,
    features = context$features,
    samples = context$samples,
    package_versions = .package_versions(packages),
    autozyme = list(
      version = as.character(utils::packageVersion("autozyme")),
      remote_sha = description$RemoteSha %||% NA_character_,
      expected_sha = Shennong:::.sn_autozyme_expected_sha,
      source_match = check$autozyme_source_match[[1L]],
      patch_source_match = check$autozyme_patch_source_match[[1L]],
      installed_upstream = check$installed_version[[1L]],
      tested_upstream = check$tested_versions[[1L]]
    ),
    acceleration = acceleration
  )
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(run, result_path, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
}

if (worker) {
  .run_worker()
  quit(status = 0L)
}

.require_packages(c("jsonlite", "digest"))
time_bin <- "/usr/bin/time"
if (!file.exists(time_bin)) stop("GNU `/usr/bin/time` is required.", call. = FALSE)

.parse_peak_rss <- function(path) {
  lines <- readLines(path, warn = FALSE)
  matched <- grep("Maximum resident set size \\(kbytes\\):", lines, value = TRUE)
  if (length(matched) != 1L) return(NA_real_)
  as.numeric(trimws(sub("^.*:", "", matched[[1L]]))) / 1024
}

.run_child <- function(operation, arm, repetition, position) {
  run_dir <- tempfile(paste0("shennong-autozyme-workflow-", operation, "-"))
  dir.create(run_dir, recursive = TRUE)
  result_path <- file.path(run_dir, "result.json")
  time_path <- file.path(run_dir, "time.txt")
  log_path <- file.path(run_dir, "worker.log")
  worker_args <- c(
    "-v", "-o", shQuote(time_path),
    shQuote(file.path(R.home("bin"), "Rscript")),
    shQuote(script_file),
    "--worker",
    "--input", shQuote(input),
    "--operation", operation,
    "--arm", arm,
    "--result", shQuote(result_path),
    "--seed", as.character(seed),
    "--max-cells", as.character(max_cells),
    "--merge-samples", as.character(merge_samples),
    "--pca-features", as.character(pca_features),
    "--pca-dims", as.character(pca_dims),
    "--go-genes", as.character(go_genes),
    "--assay", assay_option,
    "--layer", layer_option,
    "--reduction", reduction_option,
    "--species", species,
    "--orgdb", orgdb
  )
  if (nzchar(sample_by_option)) {
    worker_args <- c(worker_args, "--sample-by", sample_by_option)
  }
  if (nzchar(label_by_option)) {
    worker_args <- c(worker_args, "--label-by", label_by_option)
  }
  status <- system2(
    time_bin,
    args = worker_args,
    stdout = log_path,
    stderr = log_path,
    wait = TRUE
  )
  if (!identical(as.integer(status), 0L) || !file.exists(result_path)) {
    log <- if (file.exists(log_path)) readLines(log_path, warn = FALSE) else character()
    stop(
      "Worker failed for ", operation, "/", arm, " repetition ", repetition,
      ":\n", paste(log, collapse = "\n"),
      call. = FALSE
    )
  }
  run <- jsonlite::read_json(result_path, simplifyVector = TRUE)
  canonical <- readRDS(run$canonical_path)
  run$canonical_path <- NULL
  run$peak_rss_mb <- .parse_peak_rss(time_path)
  run$repetition <- repetition
  run$order_position <- position
  unlink(run_dir, recursive = TRUE, force = TRUE)
  list(run = run, canonical = canonical)
}

.flatten_run <- function(run) {
  data.frame(
    operation = run$operation,
    arm = run$arm,
    patch = run$patch,
    repetition = as.integer(run$repetition),
    order_position = as.integer(run$order_position),
    eligible = isTRUE(run$eligible),
    patch_active_observed = isTRUE(run$patch_active_observed),
    patch_active_after = isTRUE(run$patch_active_after),
    target_intersection = isTRUE(run$target_intersection),
    fast_path_hit = NA,
    elapsed_seconds = as.numeric(run$elapsed_seconds),
    peak_rss_mb = as.numeric(run$peak_rss_mb),
    cells = as.integer(run$cells),
    features = as.integer(run$features),
    samples = as.integer(run$samples),
    output_fingerprint = run$output_fingerprint,
    equivalent_to_direct = isTRUE(run$equivalent_to_direct),
    max_abs_diff_to_direct = as.numeric(run$max_abs_diff_to_direct %||% NA_real_),
    autozyme_version = run$autozyme$version,
    autozyme_remote_sha = run$autozyme$remote_sha %||% NA_character_,
    autozyme_expected_sha = run$autozyme$expected_sha,
    package_versions = paste(
      paste(names(run$package_versions), unlist(run$package_versions), sep = "="),
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
}

arm_rotations <- list(
  c("direct_upstream", "shennong_autozyme_off", "shennong_autozyme_on"),
  c("shennong_autozyme_on", "shennong_autozyme_off", "direct_upstream"),
  c("shennong_autozyme_off", "direct_upstream", "shennong_autozyme_on")
)

all_runs <- list()
run_index <- 0L
all_equivalent <- TRUE
repetition_ids <- seq.int(repetition_start, length.out = repetitions)
for (repetition in repetition_ids) {
  arm_order <- arm_rotations[[(repetition - 1L) %% length(arm_rotations) + 1L]]
  operation_order <- if (repetition %% 2L) operations else rev(operations)
  for (operation in operation_order) {
    message(
      "Benchmarking ", operation, " repetition/order ", repetition,
      " (", match(repetition, repetition_ids), "/", repetitions, ")."
    )
    group <- list()
    for (position in seq_along(arm_order)) {
      arm <- arm_order[[position]]
      message("  ", arm, " (position ", position, ")")
      group[[arm]] <- .run_child(operation, arm, repetition, position)
    }
    reference <- group$direct_upstream$canonical
    tolerance <- .operation_specs[[operation]]$tolerance
    for (arm in names(group)) {
      comparison <- if (identical(arm, "direct_upstream")) {
        list(pass = TRUE, max_abs_diff = 0, reason = "reference arm")
      } else {
        .compare_canonical(reference, group[[arm]]$canonical, tolerance)
      }
      group[[arm]]$run$equivalent_to_direct <- isTRUE(comparison$pass)
      group[[arm]]$run$max_abs_diff_to_direct <- comparison$max_abs_diff
      group[[arm]]$run$equivalence_reason <- comparison$reason
      all_equivalent <- all_equivalent && isTRUE(comparison$pass)
      run_index <- run_index + 1L
      all_runs[[run_index]] <- group[[arm]]$run
    }
    rm(group, reference)
    invisible(gc())
  }
}

runs_table <- do.call(rbind, lapply(all_runs, .flatten_run))
summary_rows <- lapply(operations, function(operation) {
  selected <- runs_table[runs_table$operation == operation, , drop = FALSE]
  medians <- stats::aggregate(
    selected[, c("elapsed_seconds", "peak_rss_mb"), drop = FALSE],
    list(arm = selected$arm),
    stats::median,
    na.rm = TRUE
  )
  seconds <- stats::setNames(medians$elapsed_seconds, medians$arm)
  rss <- stats::setNames(medians$peak_rss_mb, medians$arm)
  on_runs <- selected[selected$arm == "shennong_autozyme_on", , drop = FALSE]
  data.frame(
    operation = operation,
    patch = .operation_specs[[operation]]$patch,
    target = .operation_specs[[operation]]$target,
    target_intersection = .operation_specs[[operation]]$target_intersection,
    intersection_evidence = .operation_specs[[operation]]$intersection_evidence,
    benchmark_value = .operation_specs[[operation]]$benchmark_value,
    repetitions = repetitions,
    cells = unique(selected$cells)[[1L]],
    features = unique(selected$features)[[1L]],
    samples = unique(selected$samples)[[1L]],
    direct_seconds = seconds[["direct_upstream"]],
    shennong_off_seconds = seconds[["shennong_autozyme_off"]],
    shennong_on_seconds = seconds[["shennong_autozyme_on"]],
    workflow_speedup = seconds[["shennong_autozyme_off"]] / seconds[["shennong_autozyme_on"]],
    direct_peak_rss_mb = rss[["direct_upstream"]],
    shennong_off_peak_rss_mb = rss[["shennong_autozyme_off"]],
    shennong_on_peak_rss_mb = rss[["shennong_autozyme_on"]],
    all_equivalent = all(selected$equivalent_to_direct),
    patch_eligible_all_on_runs = all(on_runs$eligible),
    patch_activation_observed_any = any(on_runs$patch_active_observed),
    fast_path_hit_proven = FALSE,
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)

fixture_sha256 <- digest::digest(file = input, algo = "sha256", serialize = FALSE)
report <- list(
  schema_version = "shennong.autozyme-workflow-benchmark/v1",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  input = list(path = input, sha256 = fixture_sha256),
  methodology = list(
    arms = c("direct_upstream", "shennong_autozyme_off", "shennong_autozyme_on"),
    fresh_process_per_arm = TRUE,
    repetitions_per_arm = repetitions,
    repetition_ids = repetition_ids,
    ordering = "three rotating arm orders; operation order reverses on even repetitions",
    timing = "system.time elapsed around the analytical operation",
    memory = "GNU time maximum resident set size for the complete worker",
    summary = "median",
    activation_warning = paste(
      "patch eligibility and activation are recorded, but neither proves that",
      "the patch's input-specific fast-path guard accepted the call"
    )
  ),
  configuration = list(
    seed = seed,
    repetition_start = repetition_start,
    max_cells = max_cells,
    merge_samples = merge_samples,
    assay = assay_option,
    layer = layer_option,
    reduction = reduction_option,
    pca_features = pca_features,
    pca_dims = pca_dims,
    sample_by = sample_by_option,
    label_by = label_by_option,
    species = species,
    orgdb = orgdb,
    go_genes = go_genes,
    operations = operations
  ),
  environment = list(
    os = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
    cpu = if (file.exists("/proc/cpuinfo")) {
      model <- sub(
        "^model name[[:space:]]*:[[:space:]]*",
        "",
        grep("^model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)
      )
      if (length(model)) model[[1L]] else NA_character_
    } else {
      NA_character_
    },
    logical_cores = parallel::detectCores(logical = TRUE),
    r = R.version.string
  ),
  summary = summary_table,
  runs = all_runs,
  all_equivalent = all_equivalent
)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  report,
  output,
  dataframe = "rows",
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  na = "null",
  digits = 10
)
utils::write.csv(summary_table, sub("[.]json$", ".summary.csv", output), row.names = FALSE, na = "")
utils::write.csv(runs_table, sub("[.]json$", ".runs.csv", output), row.names = FALSE, na = "")

print(summary_table, row.names = FALSE)
cat("Wrote ", normalizePath(output, winslash = "/", mustWork = TRUE), "\n", sep = "")
if (!all_equivalent && !.flag("allow-mismatch")) {
  stop("One or more Shennong arms did not match the direct-upstream arm; see the written report.", call. = FALSE)
}
