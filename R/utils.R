#' Retrieve blocklist genes from SignatuR by preset gene set names
#'
#' This function uses the \code{SignatuR} package to fetch blocklist (or related)
#' genes (e.g. ribosomal, mitochondrial, heatshock, etc.). It then checks
#' those symbols via \code{HGNChelper} to remove invalid or outdated entries.
#'
#' @param genesets A character vector of predefined gene sets. Must be keys in
#'   \code{SignatuR::SignatuR$Hs$...}.
#'
#' @return A unique character vector of valid gene symbols.
#'
#' @examples
#' \dontrun{
#' block_genes <- sn_get_signatures(c("mito", "ribo", "tcr", "immunoglobulins", "pseudogenes", "noncoding"))
#' }
#' @export
sn_get_signatures <- function(
    species = "human",
    category = NULL) {
  # TODO: Add simplify arg
  check_installed(pkg = "HGNChelper")
  check_installed_github(pkg = "SignatuR", repo = "carmonalab/SignatuR")
  arg_match(arg = species, values = c("human", "mouse"))

  if (species == "human") {
    gene_list <- list(
      g1s = SignatuR::SignatuR$Hs$Programs$cellCycle.G1S,
      g2m = SignatuR::SignatuR$Hs$Programs$cellCycle.G2M,
      pseudogenes = SignatuR::SignatuR$Hs$Blocklists$Pseudogenes,
      noncoding = SignatuR::SignatuR$Hs$Blocklists$`Non-coding`,
      tcr = SignatuR::SignatuR$Hs$Compartments$TCR,
      immunoglobulins = SignatuR::SignatuR$Hs$Compartments$Immunoglobulins,
      ribo = SignatuR::SignatuR$Hs$Compartments$Ribo,
      mito = SignatuR::SignatuR$Hs$Compartments$Mito,
      heatshock = SignatuR::SignatuR$Hs$Programs$HeatShock
    )
  } else if (species == "mouse") {
    gene_list <- list(
      g1s = SignatuR::SignatuR$Mm$Programs$cellCycle.G1S,
      g2m = SignatuR::SignatuR$Mm$Programs$cellCycle.G2M,
      pseudogenes = SignatuR::SignatuR$Mm$Blocklists$Pseudogenes,
      noncoding = SignatuR::SignatuR$Mm$Blocklists$`Non-coding`,
      tcr = SignatuR::SignatuR$Mm$Compartments$TCR,
      immunoglobulins = SignatuR::SignatuR$Mm$Compartments$Immunoglobulins,
      ribo = SignatuR::SignatuR$Mm$Compartments$Ribo,
      mito = SignatuR::SignatuR$Mm$Compartments$Mito,
      heatshock = SignatuR::SignatuR$Mm$Programs$HeatShock
    )
  }

  genes <- unlist(lapply(category, function(x) {
    if (!x %in% names(gene_list)) {
      warning("[g_block_genes] Unknown gene set: ", x)
      return(NULL)
    }
    SignatuR::GetSignature(gene_list[[x]]) |> unlist()
  }))

  checked <- HGNChelper::checkGeneSymbols(genes, species = species)
  valid_genes <- checked$Suggested.Symbol[!is_na(checked$Suggested.Symbol)]

  if (length(valid_genes) < length(genes)) {
    invalid <- genes[is_na(checked$Suggested.Symbol)]
    message(
      "[g_block_genes] Removed ", length(invalid),
      " invalid symbols (e.g. ", paste(head(invalid, 3), collapse = ", "), ")"
    )
  }

  return(unique(valid_genes))
}

#' Check if files exist
#'
#' This function takes a vector of file paths as input and checks if each file exists.
#' If a file does not exist, the function returns a message indicating which file(s) do(es) not exist.
#'
#' @param x A vector of file paths.
#' @param stop A logical value indicating whether to stop the function if a file does not exist.
#' Default is TRUE.
#'
#' @return If stop is TRUE, the function stops and returns a message indicating which file(s) do(es) not exist.
#' If stop is FALSE, the function returns a vector of file paths that do not exist.
#'
#' @examples
#' sn_check_file(c("file1.txt", "file2.txt", "file3.txt"))
#' sn_check_file(c("file1.txt", "file2.txt", "file3.txt"), stop = FALSE)
#'
#' @export
#' @keywords file, existence
#' @seealso \code{\link{file.exists}}, \code{\link{stop}}
sn_check_file <- function(x, stop = TRUE) {
  not_exist <- x[!file.exists(x)]
  if (length(x = not_exist) > 0) {
    if (stop) {
      stop(
        paste0(
          "Please check the file path, the file listed below does not exist:\n"
        ),
        paste0(not_exist, collapse = "\n")
      )
    } else {
      return(not_exist)
    }
  }
}

#' @export
sn_set_path <- function(path) {
  path <- glue(path)
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  return(path)
}

check_installed_github <- function(pkg, repo) {
  check_installed(pkg = "pak")
  check_installed(pkg, action = function(pkg, ...) pak::pak(github_repo))
}
