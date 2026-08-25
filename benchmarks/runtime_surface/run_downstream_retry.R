script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
pkgload::load_all(repo_root, quiet = TRUE)

data_root <- "/mnt/resources/pbmc"
nichenet_dir <- "/mnt/resources/nichenetr"
out_dir <- file.path(repo_root, "data-local", "runtime-benchmark")
db_path <- file.path(out_dir, "shennong-usage.sqlite")
stopifnot(file.exists(db_path))

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

prepare <- function(name, pseudo_samples = 6L) {
  loaded <- log_step(paste0("load:", name), load_dataset(name))
  if (!isTRUE(loaded$ok)) return(NULL)
  object <- loaded$value
  doubled <- log_step(paste0("doublets:", name), {
    result <- sn_find_doublets(object, ncores = 4)
    if (inherits(result, "Seurat")) result else object
  })
  object <- doubled$value
  object$sample_id <- paste0("S", rep(seq_len(pseudo_samples), length.out = ncol(object)))
  object$condition <- factor(rep(c("A", "B"), length.out = ncol(object)))
  clustered <- log_step(
    paste0("cluster-unintegrated:", name),
    sn_run_cluster(object, integration_method = "unintegrated", cluster_name = "unintegrated")
  )
  if (!isTRUE(clustered$ok)) return(NULL)
  invisible(gc())
  clustered$value
}

big <- prepare("pbmc10k")
small <- prepare("pbmc3k", pseudo_samples = 3L)

if (!is.null(big)) {
  big_name <- big$dataset[[1]]
  cluster_ids <- head(as.character(sort(unique(as.numeric(as.character(big$seurat_clusters))))), 4L)
  signatures <- list()
  for (id in cluster_ids[-1]) {
    store <- paste0("markers_c", id)
    res <- log_step(paste0("de:c", id, "-vs-c", cluster_ids[[1]], ":", big_name),
                    sn_find_de(big, ident_1 = id, ident_2 = cluster_ids[[1]],
                               group_by = "seurat_clusters", method = "wilcoxon",
                               store_name = store))
    if (!isTRUE(res$ok)) next
    stored <- sn_get_result(res$value, type = "de", name = store)
    table <- stored$tables$primary
    genes <- head(table[order(-table$avg_log2FC), ]$gene, 100)
    signatures[[paste0("c", id)]] <- genes
  }

  log_step(paste0("de:pseudobulk-edger:", big_name),
           sn_find_de(big, ident_1 = cluster_ids[[1]], ident_2 = cluster_ids[[2]],
                      group_by = "seurat_clusters", sample_by = "sample_id", method = "edger",
                      store_name = "pseudobulk_edger_c0"))

  if (length(signatures)) {
    for (sm in c("ucell", "aucell", "gsva", "ssgsea", "mean")) {
      log_step(paste0("programs:", sm),
               sn_score_programs(big, signatures, method = sm,
                                 group_by = "seurat_clusters", return_object = FALSE))
    }
  }

  if ("markers_c1" %in% names(signatures)) {
    log_step(paste0("enrich:ora-gobp:", big_name),
             sn_enrich(big, analysis = "markers_c1", database = "GOBP"))
    log_step(paste0("enrich:gsea-hallmark:", big_name),
             sn_enrich(big, analysis = "markers_c1", database = "Hallmark"))
  }

  log_step(paste0("regulatory:progeny:", big_name),
           sn_run_regulatory_activity(big, method = "progeny", group_by = "seurat_clusters"))
  log_step(paste0("metabolism:scmetabolism-aucell:", big_name),
           sn_run_metabolism(big, method = "scmetabolism", scoring_method = "aucell",
                             group_by = "seurat_clusters"))

  if (!is.null(small)) {
    small_name <- small$dataset[[1]]
    log_step(paste0("communication:liana:", small_name),
             sn_run_cell_communication(small, method = "liana", group_by = "seurat_clusters"))
    senders <- as.character(sort(unique(as.numeric(as.character(small$seurat_clusters)))))
    sender_cluster <- senders[[1]]
    receiver_cluster <- tail(senders, 1)
    log_step(paste0("communication:nichenet:", small_name), {
      sn_run_cell_communication(
        small,
        method = "nichenet",
        group_by = "seurat_clusters",
        sender = sender_cluster,
        receiver = receiver_cluster,
        ligand_target_matrix = file.path(nichenet_dir, "ligand_target_matrix_nsga2r_final.rds"),
        lr_network = file.path(nichenet_dir, "lr_network_human_21122021.rds")
      )
    })
    log_step(paste0("communication:multinichenet:", small_name), {
      sn_run_cell_communication(
        small,
        method = "multinichenet",
        group_by = "seurat_clusters",
        sample_by = "sample_id",
        condition_by = "condition",
        ligand_target_matrix = file.path(nichenet_dir, "ligand_target_matrix_nsga2r_final.rds"),
        lr_network = file.path(nichenet_dir, "lr_network_human_21122021.rds")
      )
    })
  }
}

sn_disable_usage_tracking()

connection <- DBI::dbConnect(RSQLite::SQLite(), db_path)
on.exit(DBI::dbDisconnect(connection), add = TRUE)
summary_table <- DBI::dbGetQuery(connection, paste(
  "SELECT workflow, COUNT(*) AS calls,",
  "ROUND(SUM(elapsed_ms)/1000.0, 2) AS total_s,",
  "ROUND(AVG(elapsed_ms)/1000.0, 3) AS mean_s,",
  "ROUND(MAX(elapsed_ms)/1000.0, 2) AS max_s",
  "FROM workflow_runs WHERE depth = 0 AND status = 'ok'",
  "GROUP BY workflow ORDER BY total_s DESC"
))
errors_table <- DBI::dbGetQuery(connection, paste(
  "SELECT workflow, error_message_redacted FROM workflow_runs",
  "WHERE depth = 0 AND status = 'error' GROUP BY workflow, error_message_redacted"
))
utils::write.csv(summary_table, file.path(out_dir, "top_functions.csv"), row.names = FALSE)
utils::write.csv(errors_table, file.path(out_dir, "errors.csv"), row.names = FALSE)
message("[summary] retry run finished")
print(summary_table[1:15, ])
