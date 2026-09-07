# Type-specific semantic validation for schema-v2 result tables.

.sn_analysis_result_builtin_types <- function() {
  c(
    "annotation", "de", "differential_abundance", "milo", "enrichment",
    "cell_communication", "regulatory_activity", "deconvolution", "cnv",
    "metabolism", "grn", "program_discovery", "program_scoring",
    "program_comparison", "state_priority", "scissor", "trajectory",
    "velocity", "fate", "spatial_features", "spatial_domains",
    "spatial_neighborhood", "spatial_communication", "spatial_integration",
    "integration_benchmark", "bulk_qc", "bulk_de", "bulk_pathway",
    "bulk_network", "bulk_survival", "bulk_clinical", "qc_assessment",
    "interpretation"
  )
}

.sn_analysis_result_type_specs <- function() {
  probability_columns <- c(
    "p_value", "adjusted_p_value", "pvalue", "p.adjust", "qvalue",
    "p_val", "p_val_adj", "PValue", "FDR", "SpatialFDR", "q_value"
  )
  probability_ranges <- stats::setNames(
    rep(list(c(0, 1)), length(probability_columns)),
    probability_columns
  )
  list(
    annotation = list(
      primary_columns = c("cell", "prediction"),
      key_columns = "cell",
      character_columns = c("cell", "prediction"),
      non_missing_columns = "prediction",
      numeric_columns = "prediction_score",
      finite_or_missing_columns = "prediction_score",
      logical_columns = "low_confidence"
    ),
    de = list(
      primary_columns = "gene",
      character_columns = "gene",
      non_missing_columns = "gene",
      numeric_columns = c(
        "log2_fold_change", "avg_log2FC", "avg_logFC", "logFC", "estimate",
        "statistic", probability_columns
      ),
      finite_or_missing_columns = c(
        "log2_fold_change", "avg_log2FC", "avg_logFC", "logFC", "estimate",
        "statistic", probability_columns
      ),
      range_columns = probability_ranges
    ),
    differential_abundance = list(
      primary_columns = c("feature", "estimate"),
      key_columns = "feature",
      character_columns = "feature",
      non_missing_columns = "feature",
      numeric_columns = c("estimate", "log2_fc", "inclusion_probability", probability_columns),
      finite_or_missing_columns = c("estimate", "log2_fc", "inclusion_probability", probability_columns),
      range_columns = c(probability_ranges, list(inclusion_probability = c(0, 1))),
      logical_columns = "credible"
    ),
    milo = list(
      primary_columns = character(),
      any_columns = c("Nhood", "feature"),
      numeric_columns = c("Nhood", "logFC", probability_columns),
      finite_or_missing_columns = c("Nhood", "logFC", probability_columns),
      range_columns = probability_ranges
    ),
    enrichment = list(
      primary_columns = character(),
      any_columns = c("ID", "Description", "feature"),
      character_columns = c("ID", "Description", "Cluster", "feature"),
      numeric_columns = c("NES", "Count", probability_columns),
      finite_or_missing_columns = c("NES", "Count", probability_columns),
      nonnegative_columns = "Count",
      range_columns = probability_ranges
    ),
    cell_communication = list(
      primary_columns = c("source", "target"),
      any_columns = c("ligand", "receptor", "pathway"),
      character_columns = c(
        "source", "target", "method", "ligand", "receptor", "pathway",
        "condition", "sample"
      ),
      non_missing_columns = c("source", "target", "method"),
      numeric_columns = c("score", "rank", "spatial_distance", probability_columns),
      finite_or_missing_columns = c("score", "rank", "spatial_distance", probability_columns),
      positive_columns = "rank",
      nonnegative_columns = "spatial_distance",
      range_columns = probability_ranges
    ),
    regulatory_activity = list(
      primary_columns = "source",
      any_columns = c("condition", "cell", "feature"),
      character_columns = c("source", "condition", "cell", "feature", "statistic"),
      non_missing_columns = "source",
      numeric_columns = c("score", "activity", probability_columns),
      finite_or_missing_columns = c("score", "activity", probability_columns),
      range_columns = probability_ranges
    ),
    deconvolution = list(
      primary_columns = character(),
      any_column_sets = list(c("sample", "cell_type", "fraction"), c("feature", "score")),
      character_columns = c("sample", "cell_type", "feature"),
      numeric_columns = c("fraction", "score"),
      finite_or_missing_columns = c("fraction", "score"),
      nonnegative_columns = "fraction"
    ),
    cnv = list(
      primary_columns = c("cell", "cnv_score", "malignant_score", "malignant_call"),
      key_columns = "cell",
      character_columns = c("cell", "sample", "malignant_call", "subclone", "method"),
      numeric_columns = c("cnv_score", "malignant_score"),
      finite_or_missing_columns = c("cnv_score", "malignant_score"),
      logical_columns = "is_reference"
    ),
    metabolism = list(
      primary_columns = c("cell", "pathway", "score"),
      key_columns = c("cell", "pathway"),
      character_columns = c("cell", "pathway", "method"),
      numeric_columns = "score",
      finite_columns = "score"
    ),
    grn = list(
      primary_columns = c("source", "target", "weight", "method"),
      key_columns = c("source", "target"),
      character_columns = c("source", "target", "method"),
      numeric_columns = "weight",
      finite_columns = "weight"
    ),
    program_discovery = list(
      primary_columns = c("cell", "program", "score"),
      key_columns = c("cell", "program"),
      character_columns = c("cell", "program", "group", "method"),
      numeric_columns = "score",
      finite_columns = "score"
    ),
    program_scoring = list(
      primary_columns = c("entity", "program", "score"),
      key_columns = c("entity", "program"),
      character_columns = c("entity", "program"),
      numeric_columns = "score",
      finite_columns = "score"
    ),
    program_comparison = list(
      primary_columns = c("group", "program", "comparison", "estimate"),
      key_columns = c("group", "program", "comparison"),
      character_columns = c("group", "program", "comparison"),
      numeric_columns = c("estimate", "mean_1", "mean_2", "n_1", "n_2", probability_columns),
      finite_or_missing_columns = c("estimate", "mean_1", "mean_2", "n_1", "n_2", probability_columns),
      nonnegative_columns = c("n_1", "n_2"),
      integer_columns = c("n_1", "n_2"),
      range_columns = probability_ranges
    ),
    trajectory = list(
      primary_columns = c("cell", "primary_pseudotime"),
      key_columns = "cell",
      character_columns = c("cell", "primary_lineage"),
      numeric_columns = "primary_pseudotime",
      numeric_column_patterns = c("^pseudotime_", "^weight_"),
      finite_or_missing_columns = "primary_pseudotime",
      finite_or_missing_column_patterns = c("^pseudotime_", "^weight_"),
      range_column_patterns = list("^weight_" = c(0, 1)),
      paired_missing_columns = list(c("primary_lineage", "primary_pseudotime"))
    ),
    state_priority = list(
      primary_columns = c("state", "priority_score"),
      key_columns = "state",
      character_columns = "state",
      numeric_columns = "priority_score",
      finite_columns = "priority_score"
    ),
    scissor = list(
      primary_columns = c("cell", "coefficient", "selection"),
      key_columns = "cell",
      character_columns = c("cell", "selection"),
      non_missing_columns = "selection",
      numeric_columns = "coefficient",
      finite_columns = "coefficient",
      allowed_values = list(
        selection = c("Scissor+", "Scissor-", "Unselected")
      )
    ),
    velocity = list(
      primary_columns = c("cell", "dimension_1", "dimension_2", "velocity_1", "velocity_2"),
      key_columns = "cell",
      character_columns = c("cell", "method"),
      numeric_columns = c(
        "dimension_1", "dimension_2", "velocity_1", "velocity_2",
        "pseudotime", "confidence", "velocity_length"
      ),
      finite_columns = c("dimension_1", "dimension_2"),
      finite_or_missing_columns = c(
        "velocity_1", "velocity_2", "pseudotime", "confidence", "velocity_length"
      ),
      paired_missing_columns = list(c("velocity_1", "velocity_2"))
    ),
    fate = list(
      primary_columns = c("cell", "state", "probability"),
      key_columns = c("cell", "state"),
      character_columns = c("cell", "state"),
      numeric_columns = "probability",
      finite_columns = "probability",
      range_columns = list(probability = c(0, 1)),
      custom_validator = "fate_probability_sums"
    ),
    spatial_features = list(
      primary_columns = c("feature", "statistic", "score"),
      key_columns = "feature",
      character_columns = c("feature", "statistic"),
      numeric_columns = c("score", "rank", probability_columns),
      finite_or_missing_columns = c("score", "rank", probability_columns),
      positive_columns = "rank",
      range_columns = probability_ranges
    ),
    spatial_domains = list(
      primary_columns = c("cell", "domain"),
      key_columns = "cell",
      character_columns = c("cell", "domain"),
      non_missing_columns = c("cell", "domain")
    ),
    spatial_neighborhood = list(
      primary_columns = c("source_group", "target_group"),
      key_columns = c("source_group", "target_group"),
      character_columns = c("source_group", "target_group"),
      numeric_columns = c("observed", "expected", "z_score", probability_columns),
      finite_or_missing_columns = c("observed", "expected", "z_score", probability_columns),
      nonnegative_columns = c("observed", "expected"),
      range_columns = probability_ranges
    ),
    spatial_communication = list(
      primary_columns = c("source", "target", "spatial_distance", "within_distance"),
      any_columns = c("ligand", "receptor", "pathway"),
      character_columns = c("source", "target", "ligand", "receptor", "pathway"),
      numeric_columns = c("score", "spatial_distance", probability_columns),
      finite_or_missing_columns = c("score", "spatial_distance", probability_columns),
      nonnegative_columns = "spatial_distance",
      logical_columns = "within_distance",
      range_columns = probability_ranges
    ),
    spatial_integration = list(
      primary_columns = "cell",
      key_columns = "cell",
      character_columns = "cell",
      custom_validator = "spatial_integration_dimensions"
    ),
    integration_benchmark = list(
      primary_columns = c("run_id", "benchmark_id"),
      key_columns = "run_id",
      character_columns = c("run_id", "benchmark_id", "preprocess_id", "reason"),
      logical_columns = c("graph_only", "comparable", "supervised_label_leakage"),
      numeric_columns = c("Total", "rank_within_preprocess", "rank_overall"),
      finite_or_missing_columns = c("Total", "rank_within_preprocess", "rank_overall"),
      positive_columns = c("rank_within_preprocess", "rank_overall")
    ),
    bulk_qc = list(
      primary_columns = c("sample", "library_size", "detected_features", "outlier"),
      key_columns = "sample",
      character_columns = "sample",
      numeric_columns = c(
        "library_size", "detected_features", "median_expression", "mean_correlation",
        "library_z", "detected_z", "correlation_z"
      ),
      finite_or_missing_columns = c(
        "library_size", "detected_features", "median_expression", "mean_correlation",
        "library_z", "detected_z", "correlation_z"
      ),
      nonnegative_columns = c("library_size", "detected_features"),
      integer_columns = "detected_features",
      range_columns = list(mean_correlation = c(-1, 1)),
      logical_columns = "outlier"
    ),
    bulk_de = list(
      primary_columns = c("gene", "log2_fold_change", "p_value", "adjusted_p_value", "method"),
      key_columns = "gene",
      character_columns = c("gene", "method"),
      numeric_columns = c("log2_fold_change", "statistic", "base_mean", probability_columns),
      finite_or_missing_columns = c("log2_fold_change", "statistic", "base_mean", probability_columns),
      range_columns = probability_ranges
    ),
    bulk_pathway = list(
      primary_columns = c("sample", "pathway", "score"),
      key_columns = c("sample", "pathway"),
      character_columns = c("sample", "pathway"),
      numeric_columns = "score",
      finite_columns = "score"
    ),
    bulk_network = list(
      primary_columns = c("gene", "module"),
      key_columns = "gene",
      character_columns = c("gene", "module"),
      non_missing_columns = c("gene", "module")
    ),
    bulk_survival = list(
      primary_columns = c(
        "feature", "hazard_ratio", "conf_low", "conf_high", "p_value"
      ),
      key_columns = "feature",
      character_columns = "feature",
      numeric_columns = c("hazard_ratio", "conf_low", "conf_high", "p_value"),
      finite_columns = c("hazard_ratio", "conf_low", "conf_high", "p_value"),
      positive_columns = c("hazard_ratio", "conf_low", "conf_high"),
      range_columns = list(p_value = c(0, 1)),
      custom_validator = "bulk_survival_interval"
    ),
    bulk_clinical = list(
      primary_columns = c("feature", "clinical_variable", "test", "p_value", "n"),
      key_columns = c("feature", "clinical_variable"),
      character_columns = c("feature", "clinical_variable", "test"),
      numeric_columns = c("estimate", "statistic", "n", probability_columns),
      finite_or_missing_columns = c("estimate", "statistic", "n", probability_columns),
      nonnegative_columns = "n",
      integer_columns = "n",
      range_columns = probability_ranges,
      allowed_values = list(test = c("linear_model", "welch_t", "anova", "adjusted_anova"))
    ),
    qc_assessment = list(
      primary_columns = c("sample", "n_cells", "qc_score", "qc_label"),
      key_columns = "sample",
      character_columns = c("sample", "qc_label"),
      numeric_columns = c(
        "n_cells", "qc_score", "comparison_score", "retention_fraction",
        "low_quality_removed_fraction", "doublet_removed_fraction"
      ),
      finite_or_missing_columns = c(
        "n_cells", "qc_score", "comparison_score", "retention_fraction",
        "low_quality_removed_fraction", "doublet_removed_fraction"
      ),
      positive_columns = "n_cells",
      integer_columns = "n_cells",
      range_columns = list(
        qc_score = c(0, 100), comparison_score = c(0, 100),
        retention_fraction = c(0, 1), low_quality_removed_fraction = c(0, 1),
        doublet_removed_fraction = c(0, 1)
      )
    ),
    interpretation = list(primary_columns = character())
  )
}


.sn_result_primary_columns_label <- function(columns) {
  paste0("`", columns, "`", collapse = ", ")
}

.sn_result_primary_missing <- function(value, blank = FALSE) {
  missing <- is.na(value)
  if (isTRUE(blank) && (is.character(value) || is.factor(value))) {
    missing <- missing | !nzchar(trimws(as.character(value)))
  }
  missing
}

.sn_result_columns_matching <- function(primary, columns = character(),
                                        patterns = character()) {
  pattern_matches <- unique(unlist(lapply(patterns, function(pattern) {
    grep(pattern, names(primary), value = TRUE)
  }), use.names = FALSE))
  unique(c(intersect(columns, names(primary)), pattern_matches))
}

.sn_result_primary_semantic_errors <- function(primary, analysis_type, spec) {
  errors <- character()
  prefix <- paste0(
    "`tables$primary` for analysis type '", analysis_type, "'"
  )

  any_columns <- spec$any_columns %||% character()
  if (length(any_columns) > 0L && !any(any_columns %in% names(primary))) {
    errors <- c(
      errors,
      paste0(
        prefix, " must contain at least one of: ",
        .sn_result_primary_columns_label(any_columns), "."
      )
    )
  }

  any_column_sets <- spec$any_column_sets %||% list()
  if (length(any_column_sets) > 0L &&
      !any(vapply(any_column_sets, function(columns) {
        all(columns %in% names(primary))
      }, logical(1)))) {
    alternatives <- vapply(
      any_column_sets,
      .sn_result_primary_columns_label,
      character(1)
    )
    errors <- c(
      errors,
      paste0(
        prefix, " must contain one complete column set: ",
        paste(alternatives, collapse = " or "), "."
      )
    )
  }

  key_columns <- spec$key_columns %||% character()
  if (length(key_columns) > 0L && all(key_columns %in% names(primary))) {
    key_values <- primary[key_columns]
    valid_key_types <- vapply(
      key_values,
      function(value) is.atomic(value) && length(value) == nrow(primary),
      logical(1)
    )
    key_label <- .sn_result_primary_columns_label(key_columns)
    if (!all(valid_key_types)) {
      errors <- c(
        errors,
        paste0(prefix, " key column(s) ", key_label, " must be atomic vectors.")
      )
    } else {
      missing_keys <- vapply(
        key_values,
        function(value) any(.sn_result_primary_missing(value, blank = TRUE)),
        logical(1)
      )
      if (any(missing_keys)) {
        errors <- c(
          errors,
          paste0(
            prefix, " key column(s) ",
            .sn_result_primary_columns_label(key_columns[missing_keys]),
            " must not contain missing or blank values."
          )
        )
      }
      if (!any(missing_keys) && anyDuplicated(key_values)) {
        errors <- c(
          errors,
          paste0(prefix, " must contain unique rows for key column(s) ", key_label, ".")
        )
      }
    }
  }

  non_missing_columns <- intersect(
    spec$non_missing_columns %||% character(),
    names(primary)
  )
  for (column in non_missing_columns) {
    if (any(.sn_result_primary_missing(primary[[column]], blank = TRUE))) {
      errors <- c(
        errors,
        paste0(
          prefix, " column `", column,
          "` must not contain missing or blank values."
        )
      )
    }
  }

  character_columns <- intersect(
    spec$character_columns %||% character(),
    names(primary)
  )
  for (column in character_columns) {
    if (!is.character(primary[[column]])) {
      errors <- c(
        errors,
        paste0(prefix, " column `", column, "` must be character.")
      )
    }
  }

  numeric_columns <- .sn_result_columns_matching(
    primary,
    spec$numeric_columns %||% character(),
    spec$numeric_column_patterns %||% character()
  )
  numeric_status <- stats::setNames(
    vapply(numeric_columns, function(column) is.numeric(primary[[column]]), logical(1)),
    numeric_columns
  )
  for (column in names(numeric_status)[!numeric_status]) {
    errors <- c(
      errors,
      paste0(prefix, " column `", column, "` must be numeric.")
    )
  }

  logical_columns <- intersect(spec$logical_columns %||% character(), names(primary))
  for (column in logical_columns) {
    if (!is.logical(primary[[column]])) {
      errors <- c(
        errors,
        paste0(prefix, " column `", column, "` must be logical.")
      )
    }
  }

  integer_columns <- intersect(spec$integer_columns %||% character(), names(primary))
  for (column in integer_columns) {
    value <- primary[[column]]
    if (is.numeric(value) && any(
      !is.na(value) & (!is.finite(value) | value != trunc(value))
    )) {
      errors <- c(
        errors,
        paste0(prefix, " column `", column, "` must contain integer-like values.")
      )
    }
  }

  finite_columns <- intersect(spec$finite_columns %||% character(), names(primary))
  for (column in finite_columns) {
    value <- primary[[column]]
    if (is.numeric(value) && any(!is.finite(value))) {
      errors <- c(
        errors,
        paste0(prefix, " column `", column, "` must contain only finite values.")
      )
    }
  }

  finite_or_missing_columns <- .sn_result_columns_matching(
    primary,
    spec$finite_or_missing_columns %||% character(),
    spec$finite_or_missing_column_patterns %||% character()
  )
  for (column in finite_or_missing_columns) {
    value <- primary[[column]]
    if (is.numeric(value) && any(is.nan(value) | is.infinite(value))) {
      errors <- c(
        errors,
        paste0(
          prefix, " column `", column,
          "` may contain missing values but not NaN or infinite values."
        )
      )
    }
  }

  positive_columns <- intersect(spec$positive_columns %||% character(), names(primary))
  for (column in positive_columns) {
    value <- primary[[column]]
    if (is.numeric(value) && any(value <= 0, na.rm = TRUE)) {
      errors <- c(
        errors,
        paste0(prefix, " column `", column, "` must contain only positive values.")
      )
    }
  }

  nonnegative_columns <- intersect(
    spec$nonnegative_columns %||% character(),
    names(primary)
  )
  for (column in nonnegative_columns) {
    value <- primary[[column]]
    if (is.numeric(value) && any(value < 0, na.rm = TRUE)) {
      errors <- c(
        errors,
        paste0(prefix, " column `", column, "` must contain only non-negative values.")
      )
    }
  }

  range_columns <- spec$range_columns %||% list()
  for (column in intersect(names(range_columns), names(primary))) {
    value <- primary[[column]]
    bounds <- range_columns[[column]]
    if (is.numeric(value) && any(value < bounds[[1]] | value > bounds[[2]], na.rm = TRUE)) {
      errors <- c(
        errors,
        paste0(
          prefix, " column `", column, "` must lie in [",
          bounds[[1]], ", ", bounds[[2]], "]."
        )
      )
    }
  }

  range_column_patterns <- spec$range_column_patterns %||% list()
  for (pattern in names(range_column_patterns)) {
    bounds <- range_column_patterns[[pattern]]
    for (column in grep(pattern, names(primary), value = TRUE)) {
      value <- primary[[column]]
      if (is.numeric(value) &&
          any(value < bounds[[1]] | value > bounds[[2]], na.rm = TRUE)) {
        errors <- c(
          errors,
          paste0(
            prefix, " column `", column, "` must lie in [",
            bounds[[1]], ", ", bounds[[2]], "]."
          )
        )
      }
    }
  }

  paired_missing_columns <- spec$paired_missing_columns %||% list()
  for (columns in paired_missing_columns) {
    if (length(columns) == 2L && all(columns %in% names(primary))) {
      first_missing <- .sn_result_primary_missing(primary[[columns[[1L]]]], blank = TRUE)
      second_missing <- .sn_result_primary_missing(primary[[columns[[2L]]]], blank = TRUE)
      if (any(xor(first_missing, second_missing))) {
        errors <- c(
          errors,
          paste0(
            prefix, " columns `", columns[[1L]], "` and `", columns[[2L]],
            "` must be missing for the same rows."
          )
        )
      }
    }
  }

  allowed_values <- spec$allowed_values %||% list()
  for (column in intersect(names(allowed_values), names(primary))) {
    value <- primary[[column]]
    allowed <- allowed_values[[column]]
    invalid <- !is.na(value) & !value %in% allowed
    if (any(invalid)) {
      errors <- c(
        errors,
        paste0(
          prefix, " column `", column, "` must contain only: ",
          paste(allowed, collapse = ", "), "."
        )
      )
    }
  }

  custom_validator <- spec$custom_validator %||% NULL
  if (identical(custom_validator, "fate_probability_sums") &&
      all(c("cell", "probability") %in% names(primary)) &&
      is.numeric(primary$probability) && nrow(primary) > 0L &&
      !anyNA(primary$cell) && !anyNA(primary$probability) &&
      all(is.finite(primary$probability))) {
    sums <- tapply(primary$probability, primary$cell, sum)
    if (any(abs(sums - 1) > 1e-6)) {
      errors <- c(errors, paste0(prefix, " probabilities must sum to 1 for every cell."))
    }
  }
  if (identical(custom_validator, "bulk_survival_interval") &&
      all(c("hazard_ratio", "conf_low", "conf_high") %in% names(primary)) &&
      all(vapply(primary[c("hazard_ratio", "conf_low", "conf_high")], is.numeric, logical(1)))) {
    invalid <- with(
      primary,
      !is.na(hazard_ratio) & !is.na(conf_low) & !is.na(conf_high) &
        (conf_low > hazard_ratio | hazard_ratio > conf_high)
    )
    if (any(invalid)) {
      errors <- c(
        errors,
        paste0(prefix, " confidence intervals must satisfy conf_low <= hazard_ratio <= conf_high.")
      )
    }
  }
  if (identical(custom_validator, "spatial_integration_dimensions")) {
    dimensions <- setdiff(names(primary), "cell")
    numeric_dimensions <- dimensions[vapply(primary[dimensions], is.numeric, logical(1))]
    if (length(dimensions) < 2L || length(numeric_dimensions) != length(dimensions)) {
      errors <- c(
        errors,
        paste0(prefix, " must contain at least two numeric embedding dimensions after `cell`.")
      )
    } else if (any(!is.finite(as.matrix(primary[numeric_dimensions])))) {
      errors <- c(errors, paste0(prefix, " embedding dimensions must contain only finite values."))
    }
  }

  errors
}
