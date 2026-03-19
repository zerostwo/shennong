#!/usr/bin/env Rscript

if (!requireNamespace("Shennong", quietly = TRUE)) {
  stop("Package 'Shennong' is not installed.")
}

library(Shennong)

required_functions <- c(
  "sn_load_data",
  "sn_initialize_seurat_object",
  "sn_filter_cells",
  "sn_filter_genes",
  "sn_run_cluster",
  "sn_find_de",
  "sn_plot_dim",
  "sn_plot_dot",
  "sn_enrich",
  "sn_calculate_composition",
  "sn_calculate_lisi",
  "sn_get_codex_skill_path",
  "sn_install_codex_skill"
)

missing_functions <- required_functions[!vapply(required_functions, exists, logical(1), mode = "function")]
if (length(missing_functions) > 0) {
  stop("Missing expected functions: ", paste(missing_functions, collapse = ", "))
}

skill_dir <- sn_get_codex_skill_path()
if (!dir.exists(skill_dir)) {
  stop("Bundled skill directory was not found: ", skill_dir)
}

cat("Shennong skill smoke test passed.\n")
cat("Bundled skill path: ", skill_dir, "\n", sep = "")
