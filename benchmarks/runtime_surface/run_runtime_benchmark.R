script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
pkgload::load_all(repo_root, quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
stage_arg <- args[grepl("^--stages=", args)]
stages <- if (length(stage_arg)) {
  strsplit(sub("^--stages=", "", stage_arg), ",")[[1]]
} else {
  c("A", "B", "C", "D", "E")
}
with_cnv <- "--with-cnv" %in% args

data_root <- "/mnt/resources/pbmc"
nichenet_dir <- "/mnt/resources/nichenetr"
out_dir <- file.path(repo_root, "data-local", "runtime-benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
db_path <- file.path(out_dir, "shennong-usage.sqlite")
if (file.exists(db_path)) file.remove(db_path)

log_step <- function(name, expr) {
  t0 <- Sys.time()
  result <- tryCatch(
    list(ok = TRUE, value = force(expr)),
    error = function(e) list(ok = FALSE, error = conditionMessage(e))
  )
  seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- if (isTRUE(result$ok)) "OK" else "ERROR"
  suffix <- if (isTRUE(result$ok)) "" else paste0(" :", result$error)
  message(sprintf("[step] %-46s %-5s %9.1fs%s", name, status, seconds, suffix))
  flush.console()
  invisible(result)
}

options(future.globals.maxSize = 20 * 1024^3)
set.seed(717)

sn_enable_usage_tracking(db_path, mode = "benchmark")

dataset_names <- c("pbmc1k", "pbmc3k", "pbmc4k", "pbmc8k", "pbmc10k", "pbmc10k_5p")

load_dataset <- function(name) {
  h5 <- file.path(data_root, name, "alignments", "outs", "filtered_feature_bc_matrix.h5")
  counts <- tryCatch(
    sn_read(h5, format = "10x"),
    error = function(e) Seurat::Read10X_h5(h5)
  )
  object <- sn_initialize_seurat_object(counts, project = name)
  object$dataset <- name
  object
}

keep <- new.env()

if ("A" %in% stages) {
  for (name in dataset_names) {
    loaded <- log_step(paste0("load:", name), load_dataset(name))
    if (!isTRUE(loaded$ok)) next
    object <- loaded$value
    doubled <- log_step(paste0("doublets:", name), {
      result <- sn_find_doublets(object, ncores = 4)
      if (inherits(result, "Seurat")) result else object
    })
    object <- doubled$value
    clustered <- log_step(
      paste0("cluster-unintegrated:", name),
      sn_run_cluster(object, integration_method = "unintegrated", cluster_name = "unintegrated")
    )
    if (name %in% c("pbmc3k", "pbmc10k_5p") && isTRUE(clustered$ok)) {
      keep[[name]] <- clustered$value
    }
    rm(object); invisible(gc())
  }
}

if ("B" %in% stages) {
  for (name in c("pbmc3k", "pbmc10k_5p")) {
    object <- keep[[name]]
    if (is.null(object)) next
    for (nm in c("sctransform", "scran")) {
      log_step(
        paste0("cluster-", nm, ":", name),
        sn_run_cluster(object, integration_method = "unintegrated",
                       normalization_method = nm, cluster_name = nm)
      )
    }
  }
}

if ("C" %in% stages) {
  small <- keep$pbmc3k
  if (!is.null(small)) {
    for (alg in c("leiden", "slm", "louvain_multilevel")) {
      log_step(
        paste0("cluster-", alg, ":pbmc3k"),
        sn_run_cluster(small, integration_method = "unintegrated",
                       cluster_algorithm = alg, cluster_name = alg)
      )
    }
  }
}

if ("D" %in% stages) {
  integration_results <- list()
  for (m in c("unintegrated", "harmony", "coralysis", "seurat_rpca", "seurat_cca")) {
    res <- log_step(
      paste0("integration:", m),
      {
        objects <- lapply(c("pbmc1k", "pbmc3k", "pbmc4k"), load_dataset)
        merged <- merge(objects[[1]], y = list(objects[[2]], objects[[3]]),
                        add.cell.ids = c("pbmc1k", "pbmc3k", "pbmc4k"))
        merged$batch <- merged$dataset
        sn_run_cluster(merged, batch = "batch", integration_method = m,
                       cluster_name = paste0("c_", m))
      }
    )
    if (isTRUE(res$ok)) integration_results[[m]] <- res$value
    invisible(gc())
  }
  harmony_object <- integration_results$harmony
  if (!is.null(harmony_object)) {
    log_step("metrics:assess-harmony", {
      label_column <- intersect("c_harmony", colnames(harmony_object[[]]))
      label_column <- if (length(label_column)) label_column[[1]] else "seurat_clusters"
      sn_assess_integration(harmony_object, batch_by = "batch", label_by = label_column,
                            max_cells = 5000)
      invisible(NULL)
    })
  }
  rm(combined, integration_results, harmony_object)
  invisible(gc())
}

if ("E" %in% stages) {
  for (name in c("pbmc10k_5p", "pbmc10k", "pbmc8k", "pbmc3k")) {
    if (is.null(keep[[name]])) {
      loaded <- log_step(paste0("load:", name), load_dataset(name))
      if (!isTRUE(loaded$ok)) next
      object <- loaded$value
      doubled <- log_step(paste0("doublets:", name), {
        result <- sn_find_doublets(object, ncores = 4)
        if (inherits(result, "Seurat")) result else object
      })
      object <- doubled$value
      clustered <- log_step(
        paste0("cluster-unintegrated:", name),
        sn_run_cluster(object, integration_method = "unintegrated", cluster_name = "unintegrated")
      )
      if (isTRUE(clustered$ok)) keep[[name]] <- clustered$value
    }
  }

  big <- if (!is.null(keep$pbmc10k_5p)) {
    keep$pbmc10k_5p
  } else if (!is.null(keep$pbmc10k)) {
    keep$pbmc10k
  } else {
    keep$pbmc8k
  }
  small <- keep$pbmc3k

  if (!is.null(big)) {
    big$sample_id <- paste0("S", rep(seq_len(6), length.out = ncol(big)))
    big$condition <- factor(rep(c("A", "B"), length.out = ncol(big)))

    log_step("de:c0-vs-rest:pbmc10k_5p",
             sn_find_de(big, ident_1 = "0", group_by = "seurat_clusters",
                        method = "wilcoxon", store_name = "markers_c0"))

    cluster_ids <- head(as.character(sort(unique(as.numeric(as.character(big$seurat_clusters))))), 4L)
    signatures <- list()
    for (id in cluster_ids[-1]) {
      res <- log_step(paste0("de:c", id, "-vs-rest:pbmc10k_5p"),
                      sn_find_de(big, ident_1 = id, group_by = "seurat_clusters",
                                 method = "wilcoxon", store_name = paste0("markers_c", id)))
      if (!isTRUE(res$ok)) next
      stored <- sn_get_result(res$value, type = "de", name = paste0("markers_c", id))
      table <- stored$tables$primary
      genes <- head(table[order(-table$avg_log2FC), ]$gene, 100)
      signatures[[paste0("c", id)]] <- genes
    }

    log_step("de:pseudobulk-edger:pbmc10k_5p",
             sn_find_de(big, ident_1 = "0", group_by = "seurat_clusters",
                        sample_by = "sample_id", method = "edger",
                        store_name = "pseudobulk_edger_c0"))

    if (length(signatures)) {
      for (sm in c("ucell", "aucell", "gsva", "ssgsea", "mean")) {
        log_step(paste0("programs:", sm),
                 sn_score_programs(big, signatures, method = sm,
                                   group_by = "seurat_clusters", return_object = FALSE))
      }
    }

    log_step("enrich:ora-gobp",
             tryCatch(sn_enrich(big, analysis = "markers_c0", database = "GOBP"),
                      error = function(e) stop(conditionMessage(e))))
    log_step("enrich:gsea-hallmark",
             sn_enrich(big, analysis = "markers_c0", database = "Hallmark"))

    for (rm_method in c("dorothea", "progeny")) {
      log_step(paste0("regulatory:", rm_method),
               sn_run_regulatory_activity(big, method = rm_method,
                                          group_by = "seurat_clusters"))
    }

    log_step("metabolism:scmetabolism",
             sn_run_metabolism(big, method = "scmetabolism",
                               group_by = "seurat_clusters"))

    for (ab_method in c("propeller", "milo")) {
      log_step(paste0("abundance:", ab_method),
               sn_test_abundance(big, method = ab_method,
                                 sample_by = "sample_id",
                                 condition_by = "condition",
                                 cell_type_by = "seurat_clusters"))
    }
  }

  if (!is.null(small)) {
    for (cc_method in c("liana", "cellchat")) {
      log_step(paste0("communication:", cc_method),
               sn_run_cell_communication(small, method = cc_method,
                                         group_by = "seurat_clusters"))
    }
    log_step("communication:nichenet", {
      sn_run_cell_communication(small, method = "nichenet", group_by = "seurat_clusters",
                                ligand_target_matrix = file.path(nichenet_dir, "ligand_target_matrix_nsga2r_final.rds"),
                                lr_network = file.path(nichenet_dir, "lr_network_human_21122021.rds"))
    })
    log_step("communication:multinichenet", {
      sn_run_cell_communication(small, method = "multinichenet", group_by = "seurat_clusters",
                                ligand_target_matrix = file.path(nichenet_dir, "ligand_target_matrix_nsga2r_final.rds"),
                                lr_network = file.path(nichenet_dir, "lr_network_human_21122021.rds"))
    })
    log_step("trajectory:slingshot",
             sn_run_trajectory(small, method = "slingshot", test_dynamic = FALSE))
  }

  if (with_cnv && !is.null(big)) {
    cells <- sample(colnames(big), min(2000L, ncol(big)))
    subset_obj <- subset(big, cells = cells)
    log_step("cnv:copykat-2k-subset", sn_run_cnv(subset_obj, method = "copykat"))
    rm(subset_obj); invisible(gc())
  }
}

sn_disable_usage_tracking()

connection <- DBI::dbConnect(RSQLite::SQLite(), db_path)
on.exit(DBI::dbDisconnect(connection), add = TRUE)
summary_table <- DBI::dbGetQuery(connection, paste(
  "SELECT workflow, category, COUNT(*) AS calls,",
  "ROUND(SUM(elapsed_ms)/1000.0, 2) AS total_s,",
  "ROUND(AVG(elapsed_ms)/1000.0, 3) AS mean_s,",
  "ROUND(MAX(elapsed_ms)/1000.0, 2) AS max_s",
  "FROM workflow_runs WHERE depth = 0 AND status = 'ok'",
  "GROUP BY workflow ORDER BY total_s DESC"
))
errors_table <- DBI::dbGetQuery(connection, paste(
  "SELECT workflow, error_class, error_message_redacted FROM workflow_runs",
  "WHERE depth = 0 AND status = 'error'"
))
utils::write.csv(summary_table, file.path(out_dir, "top_functions.csv"), row.names = FALSE)
utils::write.csv(errors_table, file.path(out_dir, "errors.csv"), row.names = FALSE)
message("[summary] top functions by total elapsed time (depth-0 calls):")
print(head(summary_table, 20))
message("[summary] failed steps:")
print(errors_table)
