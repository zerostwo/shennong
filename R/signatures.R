.sn_signature_catalog <- function(source = c("package", "registry"), registry_path = NULL) {
  source <- match.arg(source)

  if (source == "package") {
    if (exists("shennong_signature_catalog", inherits = TRUE)) {
      return(get("shennong_signature_catalog", inherits = TRUE))
    }

    data_env <- new.env(parent = emptyenv())
    utils::data("shennong_signature_catalog", package = "Shennong", envir = data_env)
    return(get("shennong_signature_catalog", envir = data_env, inherits = FALSE))
  }

  .sn_read_signature_registry(path = registry_path)
}

.sn_signature_tree <- function(source = c("package", "registry"), registry_path = NULL) {
  .sn_signature_catalog(source = source, registry_path = registry_path)$tree
}

.sn_signature_registry_path <- function(path = NULL, create = FALSE) {
  if (!is_null(path)) {
    return(path.expand(path))
  }

  candidates <- unique(c(
    file.path(getwd(), "data-raw", "shennong_signature_registry.json"),
    file.path(.sn_namespace_path(), "data-raw", "shennong_signature_registry.json")
  ))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(normalizePath(existing[[1]], winslash = "/", mustWork = TRUE))
  }

  if (isTRUE(create)) {
    dir_path <- file.path(getwd(), "data-raw")
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    return(file.path(dir_path, "shennong_signature_registry.json"))
  }

  stop(
    "Could not locate `data-raw/shennong_signature_registry.json`. ",
    "Supply `registry_path` explicitly when using signature-maintenance helpers.",
    call. = FALSE
  )
}

.sn_read_signature_registry <- function(path = NULL) {
  registry_path <- .sn_signature_registry_path(path = path)
  jsonlite::read_json(registry_path, simplifyVector = FALSE)
}

.sn_write_signature_registry <- function(registry, path = NULL) {
  registry_path <- .sn_signature_registry_path(path = path, create = TRUE)
  jsonlite::write_json(
    registry,
    path = registry_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  invisible(normalizePath(registry_path, winslash = "/", mustWork = FALSE))
}

.sn_signature_normalize_key <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", x))
}

.sn_signature_as_character <- function(x) {
  if (is_null(x)) {
    return(character(0))
  }

  unique(unname(as.character(unlist(x, use.names = FALSE))))
}

.sn_signature_species_node_name <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  if (identical(species, "human")) "Hs" else "Mm"
}

.sn_signature_signatur_to_genes <- function(node) {
  raw_signature <- node$Get("Signature")
  if (length(raw_signature) == 0 || is.na(raw_signature) || !nzchar(raw_signature)) {
    return(character(0))
  }

  genes <- strsplit(raw_signature, ",", fixed = TRUE)[[1]]
  genes <- trimws(genes)
  genes[nzchar(genes)]
}

.sn_signature_tree_node_to_signatur <- function(parent, node) {
  children <- node$children %||% list()
  reference <- node$source %||% node$source_name %||% NA_character_
  child <- parent$AddChild(name = node$name, Reference = reference)

  if (length(children) == 0) {
    genes <- .sn_signature_as_character(node$genes)
    if (length(genes) > 0) {
      child$Signature <- paste(genes, collapse = ",")
    }
    return(invisible(child))
  }

  invisible(lapply(children, .sn_signature_tree_node_to_signatur, parent = child))
}

.sn_signature_tree_to_signatur_db <- function(tree) {
  check_installed_github("SignatuR", "carmonalab/SignatuR", reason = "to edit signatures with the upstream SignatuR API.")
  check_installed("data.tree")

  root <- data.tree::Node$new("SignatuR")
  species_map <- c(human = "Hs", mouse = "Mm")

  for (species in names(species_map)) {
    species_node <- tree[[species]]
    signatur_root <- root$AddChild(
      name = species_map[[species]],
      Reference = species_node$source_name %||% species_node$name %||% species
    )
    invisible(lapply(
      species_node$children %||% list(),
      .sn_signature_tree_node_to_signatur,
      parent = signatur_root
    ))
  }

  root
}

.sn_signature_signatur_node_to_tree <- function(node) {
  children <- node$children %||% list()
  out <- list(
    name = node$name,
    kind = if (length(children) > 0) "group" else "signature",
    children = if (length(children) > 0) {
      stats::setNames(
        lapply(children, .sn_signature_signatur_node_to_tree),
        vapply(children, function(child) child$name, character(1))
      )
    } else {
      list()
    }
  )

  reference <- node$Get("Reference")
  if (!is.null(reference) && length(reference) > 0 && !all(is.na(reference))) {
    out$source <- as.character(reference[[1]])
  }

  if (length(children) == 0) {
    out$genes <- .sn_signature_signatur_to_genes(node)
  }

  out
}

.sn_signature_signatur_db_to_tree <- function(db) {
  species_map <- c(human = "Hs", mouse = "Mm")
  tree <- list()

  for (species in names(species_map)) {
    signatur_name <- species_map[[species]]
    species_node <- db[[signatur_name]]
    if (is.null(species_node)) {
      stop(glue("The SignatuR object is missing the expected species root '{signatur_name}'."), call. = FALSE)
    }

    tree[[species]] <- list(
      name = species,
      kind = "group",
      source_name = species_node$name,
      children = stats::setNames(
        lapply(species_node$children %||% list(), .sn_signature_signatur_node_to_tree),
        vapply(species_node$children %||% list(), function(child) child$name, character(1))
      )
    )
  }

  tree
}

.sn_signature_signatur_get_node <- function(db, species, path_parts = character(0)) {
  current <- db[[.sn_signature_species_node_name(species)]]
  if (is.null(current)) {
    return(NULL)
  }

  if (length(path_parts) == 0) {
    return(current)
  }

  for (part in path_parts) {
    current <- current[[part]]
    if (is.null(current)) {
      return(NULL)
    }
  }

  current
}

.sn_signature_signatur_ensure_path <- function(db, species, path_parts, reference = "custom") {
  current_parts <- character(0)

  for (part in path_parts) {
    current_parts <- c(current_parts, part)
    existing <- .sn_signature_signatur_get_node(db, species, current_parts)
    if (!is.null(existing)) {
      next
    }

    parent_node <- .sn_signature_signatur_get_node(db, species, current_parts[-length(current_parts)])
    db <- SignatuR::AddNode(
      db = db,
      parent_node = parent_node,
      name = part,
      reference = reference %||% NA_character_
    )
  }

  db
}

.sn_signature_node_kind <- function(node) {
  node$kind %||% if (length(node$children %||% list()) > 0) "group" else "signature"
}

.sn_flatten_signature_node <- function(node,
                                       species,
                                       parent_path = NULL,
                                       include_groups = FALSE) {
  current_path <- if (is_null(parent_path) || !nzchar(parent_path)) {
    node$name
  } else {
    paste(parent_path, node$name, sep = "/")
  }
  kind <- .sn_signature_node_kind(node)
  genes <- .sn_signature_as_character(node$genes)

  rows <- list()
  if (isTRUE(include_groups) || identical(kind, "signature")) {
    rows[[length(rows) + 1]] <- tibble::tibble(
      species = species,
      path = current_path,
      name = node$name,
      kind = kind,
      n_genes = length(genes),
      genes = list(genes)
    )
  }

  children <- node$children %||% list()
  if (length(children) == 0) {
    return(dplyr::bind_rows(rows))
  }

  child_rows <- lapply(children, .sn_flatten_signature_node,
    species = species,
    parent_path = current_path,
    include_groups = include_groups
  )
  dplyr::bind_rows(c(rows, child_rows))
}

.sn_flatten_signature_tree <- function(tree, include_groups = FALSE) {
  species_names <- names(tree)
  rows <- lapply(species_names, function(species) {
    species_node <- tree[[species]]
    children <- species_node$children %||% list()
    dplyr::bind_rows(lapply(children, .sn_flatten_signature_node,
      species = species,
      parent_path = NULL,
      include_groups = include_groups
    ))
  })
  dplyr::bind_rows(rows)
}

.sn_signature_leaf_table <- function(species = c("human", "mouse"),
                                     source = c("package", "registry"),
                                     registry_path = NULL) {
  selected_species <- match.arg(species)
  tree <- .sn_signature_tree(source = source, registry_path = registry_path)
  .sn_flatten_signature_tree(tree, include_groups = FALSE) |>
    dplyr::filter(.data$species == selected_species)
}

.sn_signature_legacy_alias_map <- function() {
  c(
    g1s = "Programs/cellCycle.G1S",
    g2m = "Programs/cellCycle.G2M"
  )
}

.sn_signature_resolve_queries <- function(query, leaf_table) {
  if (length(query) == 0) {
    return(integer(0))
  }

  exact_path <- leaf_table$path
  exact_name <- leaf_table$name
  normalized_path <- .sn_signature_normalize_key(exact_path)
  normalized_name <- .sn_signature_normalize_key(exact_name)
  legacy_aliases <- .sn_signature_legacy_alias_map()

  matches <- integer(0)
  for (current_query in unique(query)) {
    current_matches <- which(exact_path == current_query)
    if (length(current_matches) == 0) {
      current_matches <- which(exact_name == current_query)
    }
    if (length(current_matches) == 0) {
      normalized_query <- .sn_signature_normalize_key(current_query)
      current_matches <- which(normalized_path == normalized_query)
      if (length(current_matches) == 0) {
        current_matches <- which(normalized_name == normalized_query)
      }
      if (length(current_matches) == 0 && normalized_query %in% names(legacy_aliases)) {
        current_matches <- which(exact_path == legacy_aliases[[normalized_query]])
      }
    }

    if (length(current_matches) == 0) {
      warning(
        glue(
          "Unknown signature query '{current_query}'. ",
          "Use `sn_list_signatures()` to inspect the available paths."
        ),
        call. = FALSE
      )
      next
    }

    if (length(current_matches) > 1) {
      stop(
        glue(
          "Signature query '{current_query}' is ambiguous. ",
          "Use one of the full paths instead: ",
          "{paste(leaf_table$path[current_matches], collapse = ', ')}."
        ),
        call. = FALSE
      )
    }

    matches <- c(matches, current_matches)
  }

  unique(matches)
}

.sn_signature_split_path <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

#' List bundled Shennong signatures
#'
#' @param species Optional species filter. Use \code{NULL} to return all
#'   species.
#' @param include_groups If \code{TRUE}, include non-leaf group nodes from the
#'   signature tree.
#' @param source One of \code{"package"} for the built package catalog or
#'   \code{"registry"} for the editable source registry under
#'   \code{data-raw/}.
#' @param registry_path Optional path to a signature registry JSON file. Used
#'   only when \code{source = "registry"}.
#'
#' @return A tibble with the available signature paths, node kinds, and gene
#'   counts.
#'
#' @examples
#' sn_list_signatures(species = "human")
#' sn_list_signatures(species = "human", include_groups = TRUE)
#' @export
sn_list_signatures <- function(species = NULL,
                               include_groups = FALSE,
                               source = c("package", "registry"),
                               registry_path = NULL) {
  source <- match.arg(source)
  signature_table <- .sn_flatten_signature_tree(
    tree = .sn_signature_tree(source = source, registry_path = registry_path),
    include_groups = include_groups
  ) |>
    dplyr::select(-"genes")

  if (is_null(species)) {
    return(signature_table)
  }

  selected_species <- rlang::arg_match(species, c("human", "mouse"))
  dplyr::filter(signature_table, .data$species == selected_species)
}

#' Retrieve bundled Shennong signature genes by category or tree path
#'
#' Shennong ships a package-owned signature catalog under \code{data/}, built
#' from the upstream \pkg{SignatuR} tree plus any package-maintained additions.
#' Queries can use either leaf names such as \code{"mito"} or full tree paths
#' such as \code{"Compartments/Mito"}.
#'
#' @param species One of \code{"human"} or \code{"mouse"}.
#' @param category A character vector of signature names or full tree paths.
#'
#' @return A unique character vector of signature gene symbols.
#'
#' @examples
#' sn_get_signatures(
#'   species = "human",
#'   category = c("mito", "Compartments/Ribo")
#' )
#' @export
sn_get_signatures <- function(species = "human",
                              category = NULL) {
  species <- rlang::arg_match(species, c("human", "mouse"))
  if (is_null(category) || length(category) == 0) {
    stop("`category` must contain at least one signature category or path.", call. = FALSE)
  }

  leaf_table <- .sn_signature_leaf_table(species = species, source = "package")
  matched_rows <- .sn_signature_resolve_queries(category, leaf_table)
  unique(unlist(leaf_table$genes[matched_rows], use.names = FALSE))
}

#' Add a signature to the editable source registry
#'
#' @param species One of \code{"human"} or \code{"mouse"}.
#' @param path Slash-delimited signature path relative to the species root, for
#'   example \code{"Programs/MyProgram/MySignature"}.
#' @param genes Character vector of gene symbols stored at the leaf node.
#' @param registry_path Optional path to the editable registry JSON file.
#' @param overwrite If \code{TRUE}, replace an existing signature at the same
#'   path.
#' @param source Optional source label recorded on the signature node.
#'
#' @return Invisibly returns the normalized registry path.
#'
#' @examples
#' \dontrun{
#' sn_add_signature(
#'   species = "human",
#'   path = "Programs/custom/MySignature",
#'   genes = c("GENE1", "GENE2")
#' )
#' }
#' @export
sn_add_signature <- function(species = "human",
                             path,
                             genes,
                             registry_path = NULL,
                             overwrite = FALSE,
                             source = "custom") {
  species <- rlang::arg_match(species, c("human", "mouse"))
  path_parts <- .sn_signature_split_path(path)
  if (length(path_parts) == 0) {
    stop("`path` must contain at least one non-empty path segment.", call. = FALSE)
  }
  if (length(genes) == 0) {
    stop("`genes` must contain at least one gene symbol.", call. = FALSE)
  }

  registry <- .sn_read_signature_registry(path = registry_path)
  signatur_db <- .sn_signature_tree_to_signatur_db(registry$tree)
  parent_parts <- path_parts[-length(path_parts)]
  signatur_db <- .sn_signature_signatur_ensure_path(
    db = signatur_db,
    species = species,
    path_parts = parent_parts,
    reference = source
  )
  parent_node <- .sn_signature_signatur_get_node(signatur_db, species, parent_parts)
  signatur_db <- SignatuR::AddSignature(
    db = signatur_db,
    node = parent_node,
    name = path_parts[[length(path_parts)]],
    signature = .sn_signature_as_character(genes),
    reference = source %||% NA_character_,
    overwrite = overwrite
  )
  registry$tree <- .sn_signature_signatur_db_to_tree(signatur_db)
  .sn_write_signature_registry(registry, path = registry_path)
}

#' Update a signature in the editable source registry
#'
#' @param species One of \code{"human"} or \code{"mouse"}.
#' @param path Slash-delimited signature path relative to the species root.
#' @param genes Optional replacement gene vector. If \code{NULL}, keep the
#'   existing genes.
#' @param registry_path Optional path to the editable registry JSON file.
#' @param rename_to Optional new terminal node name.
#' @param source Optional source label recorded on the signature node.
#'
#' @return Invisibly returns the normalized registry path.
#'
#' @examples
#' \dontrun{
#' sn_update_signature(
#'   species = "human",
#'   path = "Programs/custom/MySignature",
#'   genes = c("GENE1", "GENE2", "GENE3")
#' )
#' }
#' @export
sn_update_signature <- function(species = "human",
                                path,
                                genes = NULL,
                                registry_path = NULL,
                                rename_to = NULL,
                                source = NULL) {
  species <- rlang::arg_match(species, c("human", "mouse"))
  path_parts <- .sn_signature_split_path(path)
  if (length(path_parts) == 0) {
    stop("`path` must contain at least one non-empty path segment.", call. = FALSE)
  }

  registry <- .sn_read_signature_registry(path = registry_path)
  signatur_db <- .sn_signature_tree_to_signatur_db(registry$tree)
  existing_node <- .sn_signature_signatur_get_node(signatur_db, species, path_parts)
  if (is.null(existing_node)) {
    stop(glue("Signature path '{paste(path_parts, collapse = '/')}' was not found."), call. = FALSE)
  }
  if (length(existing_node$children %||% list()) > 0) {
    stop(
      glue("Path '{paste(path_parts, collapse = '/')}' refers to a group node, not a leaf signature."),
      call. = FALSE
    )
  }

  updated_genes <- if (is_null(genes)) {
    unique(unname(as.character(unlist(SignatuR::GetSignature(existing_node), use.names = FALSE))))
  } else {
    .sn_signature_as_character(genes)
  }
  if (length(updated_genes) == 0) {
    stop("`genes` must contain at least one gene symbol.", call. = FALSE)
  }

  parent_parts <- path_parts[-length(path_parts)]
  old_name <- path_parts[[length(path_parts)]]
  new_name <- if (!is_null(rename_to) && nzchar(rename_to)) rename_to else old_name
  parent_node <- .sn_signature_signatur_get_node(signatur_db, species, parent_parts)
  current_reference <- existing_node$Get("Reference")
  updated_reference <- source %||% if (length(current_reference) > 0 && !is.na(current_reference)) current_reference else NA_character_

  signatur_db <- SignatuR::AddSignature(
    db = signatur_db,
    node = parent_node,
    name = new_name,
    signature = updated_genes,
    reference = updated_reference,
    overwrite = TRUE
  )
  if (!identical(new_name, old_name)) {
    old_node <- .sn_signature_signatur_get_node(signatur_db, species, path_parts)
    SignatuR::RemoveSignature(old_node)
  }

  registry$tree <- .sn_signature_signatur_db_to_tree(signatur_db)
  .sn_write_signature_registry(registry, path = registry_path)
}

#' Delete a signature from the editable source registry
#'
#' @param species One of \code{"human"} or \code{"mouse"}.
#' @param path Slash-delimited signature path relative to the species root.
#' @param registry_path Optional path to the editable registry JSON file.
#'
#' @return Invisibly returns the normalized registry path.
#'
#' @examples
#' \dontrun{
#' sn_delete_signature(
#'   species = "human",
#'   path = "Programs/custom/MySignature"
#' )
#' }
#' @export
sn_delete_signature <- function(species = "human",
                                path,
                                registry_path = NULL) {
  species <- rlang::arg_match(species, c("human", "mouse"))
  path_parts <- .sn_signature_split_path(path)
  if (length(path_parts) == 0) {
    stop("`path` must contain at least one non-empty path segment.", call. = FALSE)
  }

  registry <- .sn_read_signature_registry(path = registry_path)
  signatur_db <- .sn_signature_tree_to_signatur_db(registry$tree)
  target_node <- .sn_signature_signatur_get_node(signatur_db, species, path_parts)
  if (is.null(target_node)) {
    stop(glue("Signature path '{paste(path_parts, collapse = '/')}' was not found."), call. = FALSE)
  }
  SignatuR::RemoveSignature(target_node)
  registry$tree <- .sn_signature_signatur_db_to_tree(signatur_db)
  .sn_write_signature_registry(registry, path = registry_path)
}
