# Bulk-model design and contrast construction helpers.

.sn_bulk_has_random_effect <- function(design) grepl("\\|", paste(deparse(design), collapse = ""), fixed = FALSE)

.sn_bulk_fixed_effect_formula <- function(design) {
  design_terms <- stats::terms(design)
  labels <- attr(design_terms, "term.labels")
  labels <- labels[!grepl("|", labels, fixed = TRUE)]
  fixed <- stats::reformulate(
    labels,
    intercept = attr(design_terms, "intercept")
  )
  environment(fixed) <- environment(design)
  fixed
}

.sn_bulk_design <- function(metadata, design) {
  if (!inherits(design, "formula")) stop("`design` must be a formula.", call. = FALSE)
  variables <- all.vars(design)
  missing <- setdiff(variables, colnames(metadata))
  if (length(missing) > 0L) stop("Design variable(s) missing from metadata: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  if (.sn_bulk_has_random_effect(design)) return(list(formula = design, matrix = NULL, rank = NA_integer_))
  design_matrix <- stats::model.matrix(design, data = metadata)
  rank <- qr(design_matrix)$rank
  if (rank < ncol(design_matrix)) stop("The bulk design matrix is not full rank.", call. = FALSE)
  list(formula = design, matrix = design_matrix, rank = rank)
}

.sn_bulk_validate_contrast <- function(metadata, contrast, design = NULL) {
  if (!is.character(contrast) || length(contrast) != 3L ||
      anyNA(contrast) || any(!nzchar(contrast)) || identical(contrast[[2]], contrast[[3]])) {
    stop("`contrast` must be c(variable, numerator, denominator).", call. = FALSE)
  }
  variable <- contrast[[1]]
  if (!variable %in% colnames(metadata)) stop("Contrast variable '", variable, "' was not found.", call. = FALSE)
  if (!is.null(design)) {
    if (!inherits(design, "formula")) stop("`design` must be a formula.", call. = FALSE)
    if (!variable %in% all.vars(design)) {
      stop(
        "Contrast variable '", variable,
        "' is not included in `design`; refusing an undefined comparison.",
        call. = FALSE
      )
    }
  }
  if (!is.factor(metadata[[variable]]) &&
      !is.character(metadata[[variable]]) &&
      !is.logical(metadata[[variable]])) {
    stop("Contrast variable '", variable, "' must be categorical (factor, character, or logical).", call. = FALSE)
  }
  values <- as.character(metadata[[variable]])
  missing <- setdiff(contrast[2:3], unique(values))
  if (length(missing) > 0L) stop("Contrast level(s) missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  contrast
}

.sn_bulk_contrast_estimand <- function(metadata, design, contrast) {
  metadata <- as.data.frame(metadata)
  for (column in colnames(metadata)) {
    if (is.character(metadata[[column]]) || is.logical(metadata[[column]])) {
      metadata[[column]] <- factor(metadata[[column]])
    }
  }
  reference <- metadata[1, , drop = FALSE]
  for (column in colnames(metadata)) {
    if (is.numeric(metadata[[column]])) reference[[column]] <- stats::median(metadata[[column]], na.rm = TRUE)
    if (is.factor(metadata[[column]])) reference[[column]] <- factor(levels(metadata[[column]])[[1]], levels = levels(metadata[[column]]))
  }
  numerator <- denominator <- reference
  variable <- contrast[[1]]
  if (is.factor(metadata[[variable]])) {
    numerator[[variable]] <- factor(contrast[[2]], levels = levels(metadata[[variable]]))
    denominator[[variable]] <- factor(contrast[[3]], levels = levels(metadata[[variable]]))
  } else {
    numerator[[variable]] <- contrast[[2]]
    denominator[[variable]] <- contrast[[3]]
  }
  design_frame <- stats::model.frame(design, data = metadata, na.action = stats::na.fail)
  design_terms <- stats::terms(design_frame)
  design_matrix <- stats::model.matrix(design_terms, data = design_frame)
  factor_levels <- stats::.getXlevels(design_terms, design_frame)
  contrasts <- attr(design_matrix, "contrasts")
  matrix_for <- function(profile) {
    frame <- stats::model.frame(
      design_terms,
      data = profile,
      xlev = factor_levels,
      na.action = stats::na.pass
    )
    stats::model.matrix(design_terms, data = frame, contrasts.arg = contrasts)
  }
  numerator_matrix <- matrix_for(numerator)
  denominator_matrix <- matrix_for(denominator)
  if (!identical(colnames(numerator_matrix), colnames(design_matrix)) ||
      !identical(colnames(denominator_matrix), colnames(design_matrix))) {
    stop("Could not align the requested contrast with the fitted bulk design matrix.", call. = FALSE)
  }
  vector <- as.numeric(numerator_matrix[1L, ] - denominator_matrix[1L, ])
  names(vector) <- colnames(design_matrix)
  if (any(!is.finite(vector)) || !any(abs(vector) > sqrt(.Machine$double.eps))) {
    stop("The requested bulk contrast produced an invalid or zero estimand.", call. = FALSE)
  }
  profile_values <- function(profile) {
    values <- profile[, all.vars(design), drop = FALSE]
    lapply(values, function(value) {
      value <- value[[1L]]
      if (is.factor(value)) as.character(value) else unname(value)
    })
  }
  list(
    vector = vector,
    numerator_profile = profile_values(numerator),
    denominator_profile = profile_values(denominator),
    covariate_policy = "numeric covariates at observed medians; categorical covariates at first observed factor levels"
  )
}

.sn_bulk_contrast_vector <- function(metadata, design, contrast) {
  .sn_bulk_contrast_estimand(metadata, design, contrast)$vector
}
