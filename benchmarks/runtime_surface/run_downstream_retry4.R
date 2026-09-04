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

h5 <- file.path(data_root, "pbmc10k", "alignments", "outs", "filtered_feature_bc_matrix.h5")
big <- sn_initialize_seurat_object(Seurat::Read10X_h5(h5), project = "pbmc10k")
clustered <- log_step("prep:cluster-pbmc10k",
                      sn_run_cluster(big, integration_method = "unintegrated",
                                     cluster_name = "unintegrated"))
big <- clustered$value
ids <- head(as.character(sort(unique(as.numeric(as.character(big$seurat_clusters))))), 2L)
store <- "markers_c1"
res <- log_step("de:wilcox:c1-vs-c0",
                sn_find_de(big, analysis = "markers",
                           ident_1 = ids[[2]], ident_2 = ids[[1]],
                           group_by = "seurat_clusters", method = "wilcox",
                           result_id = store))
big <- res$value
stored <- sn_get_result(big, type = "de", result_id = store)
table <- stored$tables$primary
geneset <- head(table[order(-table$avg_log2FC), ]$gene, 100)

log_step("enrich:ora-gobp",
         sn_enrich(big, analysis = "ora", source_de_result_id = store, database = "GOBP"))
log_step("enrich:gsea-hallmark",
         sn_enrich(big, analysis = "gsea", source_de_result_id = store, database = "Hallmark"))

small <- readRDS("/tmp/opencode/clustered3k.rds")
log_step("communication:nichenet-geneset", {
  sn_run_cell_communication(
    small, method = "nichenet", group_by = "seurat_clusters",
    sender = "0", receiver = "8",
    geneset = geneset,
    ligand_target_matrix = as.matrix(readRDS(file.path(nichenet_dir, "ligand_target_matrix_nsga2r_final.rds"))),
    lr_network = as.data.frame(readRDS(file.path(nichenet_dir, "lr_network_human_21122021.rds")))
  )
})

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
message("[summary] retry4 finished")
print(summary_table[1:20, ])
