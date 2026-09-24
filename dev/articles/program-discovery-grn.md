# Program Discovery and Gene Regulatory Networks

Program discovery and GRN inference answer different questions. This
article uses author-normalized expression from real malignant cells in
GSE72056. The local NMF implementation is a core workflow; GENIE3 is an
optional R backend and never receives a fabricated edge table when it is
absent.

``` r

knitr::kable(data.frame(
  workflow = c("NMF program discovery", "GENIE3 GRN", "pySCENIC / GRNBoost2"),
  mode = c("core", "extended", "external import or runner"),
  status = c(
    if (core_ready) "executed below" else if (!run_vignette) "disabled: set SHENNONG_RUN_VIGNETTES=true" else "fixture missing",
    if (genie3_ready) "executed below" else if (!identical(real_profile, "all")) "disabled: requires SHENNONG_REAL_PROFILE=all" else if (!requireNamespace("GENIE3", quietly = TRUE)) "dependency missing: GENIE3" else "fixture missing",
    "not run: no audited external result supplied"
  ),
  check.names = FALSE
))
```

| workflow | mode | status |
|:---|:---|:---|
| NMF program discovery | core | disabled: set SHENNONG_RUN_VIGNETTES=true |
| GENIE3 GRN | extended | disabled: requires SHENNONG_REAL_PROFILE=all |
| pySCENIC / GRNBoost2 | external import or runner | not run: no audited external result supplied |

## Discover latent programs in real malignant cells

The fixture retains the study’s inferred-CNV malignancy class. We select
the 300 most variable genes among malignant cells and run four NMF
programs with three seeded restarts. Every weight and activity value
printed below is fitted during the article build.

``` r

library(Shennong)

melanoma <- qs2::qs_read(melanoma_path)
malignant <- subset(
  melanoma,
  cells = colnames(melanoma)[melanoma$malignant_call == "malignant"]
)
expression <- SeuratObject::LayerData(
  malignant, assay = "RNA", layer = "data"
)
gene_variance <- Matrix::rowMeans(expression ^ 2) -
  Matrix::rowMeans(expression) ^ 2
program_features <- names(head(sort(gene_variance, decreasing = TRUE), 300))

data.frame(
  tumors = length(unique(malignant$tumor)),
  malignant_cells = ncol(malignant),
  measured_genes = nrow(malignant),
  fitted_genes = length(program_features)
)
```

``` r

programs <- sn_discover_programs(
  malignant,
  method = "nmf",
  n_programs = 4,
  result_id = "melanoma_programs",
  assay = "RNA",
  layer = "data",
  features = program_features,
  backend_control = list(
    nrun = 3L,
    max_iter = 120L,
    top_genes = 25L,
    seed = 717L
  ),
  return_object = FALSE
)

head(programs$tables$gene_weights, 20)
programs$tables$fit_diagnostics
programs$tables$stability
sn_plot_discovered_programs(programs, type = "weights", n = 10)
sn_plot_discovered_programs(programs, type = "activity")
sn_plot_discovered_programs(programs, type = "stability")
```

## Extended: infer a GENIE3 network

When the `all` profile and GENIE3 are available, the backend uses the
same real malignant-cell matrix. The selected regulators must occur
among the 300 variable input genes. A backend failure is rendered as a
failed status and is not replaced by imported or simulated edges.

``` r

candidate_regulators <- c(
  "MITF", "SOX10", "STAT1", "IRF1", "JUN", "FOS", "CEBPB", "ETS1"
)
regulators <- intersect(candidate_regulators, program_features)

genie3_attempt <- tryCatch(
  sn_run_grn(
    malignant,
    method = "genie3",
    regulators = regulators,
    group_by = "tumor",
    assay = "RNA",
    layer = "data",
    result_id = "melanoma_genie3",
    backend_control = list(
      n_features = 300L,
      n_trees = 200L,
      n_workers = 2L,
      max_edges = 1000L,
      top_targets = 30L,
      seed = 717L
    ),
    return_object = FALSE
  ),
  error = identity
)

if (inherits(genie3_attempt, "error")) {
  knitr::kable(data.frame(
    backend = "GENIE3", status = "failed",
    detail = conditionMessage(genie3_attempt), check.names = FALSE
  ))
} else {
  grn <- genie3_attempt
  print(head(grn$tables$edges, 20))
  print(head(grn$tables$regulons, 20))
  print(head(grn$tables$specificity, 20))
  print(sn_plot_regulon(grn, type = "network", n = 30))
  print(sn_plot_regulon(grn, type = "specificity"))
}
```

pySCENIC and GRNBoost2 require an explicit audited runner or imported
result. The documentation build does not silently install their Python
environments or present NMF output as a GRN result.
