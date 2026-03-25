# Build the packaged Shennong signature snapshot from SignatuR.
# Run manually during development.

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

build_catalog_from_signatur <- function() {
  if (!requireNamespace("SignatuR", quietly = TRUE)) {
    stop(
      "Install SignatuR before building the packaged signature catalog from the upstream tree.",
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
      snapshot_name = "shennong_signature_catalog",
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

shennong_signature_catalog <- build_catalog_from_signatur()

dir.create(dirname(data_path), recursive = TRUE, showWarnings = FALSE)
save(
  shennong_signature_catalog,
  file = data_path,
  version = 2,
  compress = "xz"
)
