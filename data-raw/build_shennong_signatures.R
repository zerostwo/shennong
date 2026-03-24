# Build the editable Shennong signature registry and package data catalog.
# Run manually during development.
#
# Default behavior:
# - if `data-raw/shennong_signature_registry.json` does not exist, bootstrap it
#   from the installed SignatuR package while preserving the full tree.
# - always rebuild `data/shennong_signature_catalog.rda` from the current
#   registry.
#
# Set `SHENNONG_REFRESH_SIGNATURE_REGISTRY=true` when you explicitly want to
# replace the editable registry with a fresh import from SignatuR.

refresh_from_signatur <- identical(
  tolower(Sys.getenv("SHENNONG_REFRESH_SIGNATURE_REGISTRY", unset = "false")),
  "true"
)
registry_path <- file.path("data-raw", "shennong_signature_registry.json")
data_path <- file.path("data", "shennong_signature_catalog.rda")

extract_signatur_genes <- function(node) {
  raw_signature <- tryCatch(
    SignatuR::GetSignature(node),
    error = function(e) NULL
  )
  genes <- unique(unname(as.character(unlist(raw_signature, use.names = FALSE))))
  genes[!is.na(genes) & nzchar(genes)]
}

signatur_to_tree <- function(node, source = "SignatuR") {
  children <- node$children
  child_nodes <- if (length(children) > 0) {
    stats::setNames(
      lapply(children, signatur_to_tree, source = source),
      vapply(children, function(child) child$name, character(1))
    )
  } else {
    list()
  }

  genes <- extract_signatur_genes(node)
  kind <- if (length(child_nodes) > 0) "group" else if (length(genes) > 0) "signature" else "group"

  out <- list(
    name = node$name,
    kind = kind,
    source = source,
    children = child_nodes
  )

  if (identical(kind, "signature") && length(genes) > 0) {
    out$genes <- genes
  }

  out
}

bootstrap_registry_from_signatur <- function() {
  if (!requireNamespace("SignatuR", quietly = TRUE)) {
    stop(
      "Install SignatuR before bootstrapping the signature registry from the upstream tree.",
      call. = FALSE
    )
  }

  signature_env <- new.env(parent = emptyenv())
  utils::data("SignatuR", package = "SignatuR", envir = signature_env)
  signatur_data <- get("SignatuR", envir = signature_env, inherits = FALSE)

  human_tree <- signatur_to_tree(signatur_data$Hs)
  human_tree$name <- "human"
  human_tree$source_name <- signatur_data$Hs$name

  mouse_tree <- signatur_to_tree(signatur_data$Mm)
  mouse_tree$name <- "mouse"
  mouse_tree$source_name <- signatur_data$Mm$name

  list(
    metadata = list(
      registry_name = "shennong_signature_registry",
      source_package = "SignatuR",
      source_version = as.character(utils::packageVersion("SignatuR")),
      built_on = as.character(Sys.Date())
    ),
    tree = list(
      human = human_tree,
      mouse = mouse_tree
    )
  )
}

read_registry <- function(path) {
  jsonlite::read_json(path, simplifyVector = FALSE)
}

write_registry <- function(registry, path) {
  jsonlite::write_json(
    registry,
    path = path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
}

if (isTRUE(refresh_from_signatur) || !file.exists(registry_path)) {
  registry <- bootstrap_registry_from_signatur()
  dir.create(dirname(registry_path), recursive = TRUE, showWarnings = FALSE)
  write_registry(registry, registry_path)
}

registry <- read_registry(registry_path)
shennong_signature_catalog <- list(
  metadata = registry$metadata,
  tree = registry$tree
)

dir.create(dirname(data_path), recursive = TRUE, showWarnings = FALSE)
save(
  shennong_signature_catalog,
  file = data_path,
  version = 2,
  compress = "xz"
)
