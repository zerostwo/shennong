#!/usr/bin/env Rscript
# Structural atlas gallery. Synthetic data demonstrate grammar, not paper results.
pkgload::load_all(quiet = TRUE)
out <- file.path(getwd(), "data-local", "atlas-layout-review")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(102)
m <- matrix(rnorm(48), 6, dimnames = list(paste0("State ", 1:6), paste0("Tissue ", 1:8)))
square <- stats::cor(t(m)); dimnames(square) <- rep(list(rownames(m)), 2)
other <- square * 0.65; diag(other) <- 1
sizes <- matrix(sample(0:100, 48, TRUE), 6, dimnames = dimnames(m))
ann <- data.frame(System = rep(c("Group I", "Group II"), each = 4), row.names = colnames(m))
row_ann <- data.frame(Lineage = rep(c("Lymphoid", "Myeloid", "Stromal"), each = 2), row.names = rownames(m))
rs <- stats::setNames(rep(c("Lymphoid", "Myeloid", "Stromal"), each = 2), rownames(m))
cs <- stats::setNames(rep(c("Normal", "Tumour"), each = 4), colnames(m))
plots <- list(
  upper_triangle = sn_plot_matrix(square, style = "upper_triangle", midpoint = 0, limits = c(-1, 1)),
  lower_triangle = sn_plot_matrix(square, style = "lower_triangle", midpoint = 0, limits = c(-1, 1)),
  split_triangle = sn_plot_matrix(square, style = "split_triangle", compare = other, midpoint = 0, limits = c(-1, 1), show_values = TRUE),
  circular = sn_plot_matrix(m, style = "circular", midpoint = 0),
  annotated = sn_plot_matrix(m, style = "annotated", midpoint = 0, column_annotation = ann, row_annotation = row_ann),
  dendrogram = sn_plot_matrix(m, style = "dendrogram", midpoint = 0),
  grouped_bubble = sn_plot_matrix(m, style = "grouped_bubble", size = sizes, row_split = rs, column_split = cs, midpoint = 0, value_label = "Mean expression", size_label = "Percent (%)"),
  heatmap_bar = sn_plot_matrix(m, style = "heatmap_bar", midpoint = 0, row_summary = stats::setNames(1:6, rownames(m)), summary_label = "Supplied count"),
  heatmap_density = sn_plot_matrix(m, style = "heatmap_density", midpoint = 0, column_summary = stats::setNames(c(0.1, 0.2, 0.3, 0.7, 1, 0.8, 0.4, 0.1), colnames(m)), summary_label = "Supplied density")
)
d <- data.frame(group = rep(c("Pre", "Post"), each = 36), donor = rep(paste0("D", 1:36), 2),
  IFN = c(rnorm(36, 1), rnorm(36, 2)), Cytotoxicity = c(rnorm(36, 2), rnorm(36, 1)), Proliferation = rnorm(72))
for (s in c("violin_box", "raincloud", "paired_box", "paired_slope", "split_violin", "stacked_violin")) {
  plots[[s]] <- sn_plot_distribution(d, if (s %in% c("split_violin", "stacked_violin")) c("IFN", "Cytotoxicity", "Proliferation") else "IFN",
    "group", style = s, pair_by = if (startsWith(s, "paired")) "donor" else NULL)
}
effects <- data.frame(label = c("NK state", "CD8 state", "B state", "Myeloid state", "Fibroblast"),
  estimate = c(-0.8, 0.5, 0.3, -0.3, 0.9), lower = c(-1.1, 0.2, -0.1, -0.6, 0.5), upper = c(-0.5, 0.8, 0.7, 0, 1.3))
for (s in c("forest", "diverging_bar", "lollipop")) plots[[s]] <- sn_plot_effect(effects, style = s, value_label = "Supplied coefficient")
plots$upset <- sn_plot_overlap(list(B = paste0("S", 1:12), T = paste0("S", 5:16), Myeloid = paste0("S", c(8:14, 17:20))))
e <- data.frame(source = c("A", "A", "B", "B", "C"), target = c("B", "C", "D", "E", "F"), weight = c(1, 3, 2, 1, 2), mutations = c(1, 2, 1, 2, 2))
n <- data.frame(node = LETTERS[1:6], lineage = c("Regulator", "Target", "Target", "Regulator", "Regulator", "Regulator"))
for (s in c("bipartite", "clone_tree", "chord", "rewiring")) plots[[s]] <- sn_plot_network(e, n, group_by = "lineage", style = s,
  compare_edges = if (s == "rewiring") e[c(1, 3, 5), ] else NULL,
  branch_length = if (s == "clone_tree") "mutations" else NULL, weight_label = "Supplied weight")
catalog <- sn_list_plot_styles(); catalog <- catalog[catalog$kind == "structure", ]
stopifnot(setequal(names(plots), unique(catalog$style)))
for (s in names(plots)) {
  p <- plots[[s]]
  heading <- paste(s, "| synthetic layout demonstration")
  if (inherits(p, "patchwork")) p <- p + patchwork::plot_annotation(title = heading,
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 8))) else p <- p + ggplot2::labs(title = heading)
  sn_export_figure_bundle(p, file.path(out, s), width = if (s == "rewiring") 235 else 185, height = 140, dpi = 150)
  cat("Rendered: ", s, "\n", sep = "")
}
# Original source data, separate from synthetic geometry demonstrations.
workbook <- "data-local/zhang-atlas-figures/41586_2025_9053_MOESM6_ESM.xlsx"
if (file.exists(workbook) && requireNamespace("readxl", quietly = TRUE)) {
  source <- readxl::read_excel(workbook, sheet = "Fig. 2c", skip = 1)
  roe <- vapply(source[-1], as.numeric, numeric(nrow(source))); rownames(roe) <- source[[1]]
  p <- sn_plot_matrix(roe, style = "dendrogram", limits = c(0, 3), palette = c("white", "#FCAE91", "#CB181D"),
    value_label = "Ro/e", title = "Nature Fig. 2c original values; independently recomputed clustering")
  sn_export_figure_bundle(p, file.path(out, "nature-original-values"), width = 300, height = 175)
  recovered <- as.matrix(utils::read.csv(file.path(out, "nature-original-values", "nature-original-values_data_matrix.csv"), row.names = 1, check.names = FALSE))
  stopifnot(isTRUE(all.equal(recovered, roe, tolerance = 1e-14)))
}
utils::write.csv(catalog, file.path(out, "layout-catalog.csv"), row.names = FALSE)
jsonlite::write_json(list(structural_styles = names(plots), n = length(plots),
  fixture = "synthetic inputs except nature-original-values", source_readback = "408 original Nature Ro/e values",
  semantics = "No p-values, modules, BCR trees, or network edges inferred"), file.path(out, "validation.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
cat("Gallery: ", out, "\n", sep = "")
