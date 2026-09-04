script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
pkgload::load_all(repo_root, quiet = TRUE)

data_root <- "/mnt/resources/pbmc"
nichenet_dir <- "/mnt/resources/nichenetr"
out_dir <- file.path(repo_root, "data-local", "runtime-benchmark")
db_path <- file.path(out_dir, "shennong-usage.sqlite")

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
  counts <- tryCatch(sn_read(h5, format = "10x"), error = function(e) Seurat::Read10X_h5(h5))
  object <- sn_initialize_seurat_object(counts, project = name)
  object$dataset <- name
  object
}

prepare_blocked <- function(name, n_samples) {
  loaded <- log_step(paste0("load:", name), load_dataset(name))
  if (!isTRUE(loaded$ok)) return(NULL)
  object <- loaded$value
  blocks <- cut(seq_len(ncol(object)), breaks = n_samples, labels = FALSE)
  object$sample_id <- factor(paste0("S", blocks), levels = paste0("S", seq_len(n_samples)))
  object$condition <- factor(ifelse(blocks <= floor(n_samples / 2), "A", "B"))
  clustered <- log_step(
    paste0("cluster-unintegrated:", name),
    sn_run_cluster(object, integration_method = "unintegrated", cluster_name = "unintegrated")
  )
  if (!isTRUE(clustered$ok)) return(NULL)
  invisible(gc())
  clustered$value
}

big <- prepare_blocked("pbmc10k", 6L)
small <- prepare_blocked("pbmc3k", 4L)

lr_network <- as.data.frame(readRDS(file.path(nichenet_dir, "lr_network_human_21122021.rds")))
ligand_target_matrix <- as.matrix(readRDS(file.path(nichenet_dir, "ligand_target_matrix_nsga2r_final.rds")))

if (!is.null(big)) {
  ids <- head(as.character(sort(unique(as.numeric(as.character(big$seurat_clusters))))), 2L)
  store <- "markers_c1"
  res <- log_step(paste0("de:wilcox:c1-vs-c0"),
                  sn_find_de(big, analysis = "markers",
                             ident_1 = ids[[2]], ident_2 = ids[[1]],
                             group_by = "seurat_clusters", method = "wilcox",
                             result_id = store))
  if (isTRUE(res$ok)) {
    big <- res$value
    log_step("enrich:ora-gobp", sn_enrich(big, analysis = store, database = "GOBP"))
    log_step("enrich:gsea-hallmark", sn_enrich(big, analysis = store, database = "Hallmark"))
  }
  log_step("de:pseudobulk-edger-fixed",
           sn_find_de(big, analysis = "pseudobulk",
                      ident_1 = ids[[1]], ident_2 = ids[[2]],
                      group_by = "seurat_clusters", sample_by = "sample_id",
                      method = "edgeR", result_id = "pseudobulk_edger_c0"))
}

if (!is.null(small)) {
  senders <- as.character(sort(unique(as.numeric(as.character(small$seurat_clusters)))))
  log_step("communication:nichenet-cond", {
    sn_run_cell_communication(
      small, method = "nichenet", group_by = "seurat_clusters",
      sender = senders[[1]], receiver = tail(senders, 1),
      condition_by = "condition", condition_oi = "B", condition_reference = "A",
      ligand_target_matrix = ligand_target_matrix,
      lr_network = lr_network
    )
  })
  log_step("communication:multinichenet-blocked", {
    sn_run_cell_communication(
      small, method = "multinichenet", group_by = "seurat_clusters",
      sample_by = "sample_id", condition_by = "condition",
      ligand_target_matrix = ligand_target_matrix,
      lr_network = lr_network
    )
  })
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
utils::write.csv(summary_table, file.path(out_dir, "top_functions.csv"), row.names = FALSE)
errors_table <- DBI::dbGetQuery(connection, paste(
  "SELECT workflow, error_message_redacted FROM workflow_runs",
  "WHERE depth = 0 AND status = 'error' GROUP BY workflow, error_message_redacted"
))
utils::write.csv(errors_table, file.path(out_dir, "errors.csv"), row.names = FALSE)
message("[summary] retry3 finished")
print(summary_table[1:18, ])
