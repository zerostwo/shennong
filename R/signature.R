#' Retrieve blocklist genes from SignatuR by preset gene set names
#'
#' This function uses the \code{SignatuR} package to fetch blocklist (or related)
#' genes (e.g. ribosomal, mitochondrial, heatshock, etc.). It then checks
#' those symbols via \code{HGNChelper} to remove invalid or outdated entries.
#'
#' @param genesets A character vector of predefined gene sets. Must be keys in
#'   \code{SignatuR$Hs$...}.
#'
#' @return A unique character vector of valid gene symbols.
#'
#' @examples
#' \dontrun{
#' block_genes <- sn_get_signatures(c("mito", "ribo", "tcr", "immunoglobulins", "pseudogenes", "noncoding"))
#' }
#' @export
sn_get_signatures <- function(species = "human",
                              category = NULL) {
  # TODO: Add simplify arg
  check_installed(pkg = "HGNChelper")
  check_installed_github(pkg = "SignatuR", repo = "carmonalab/SignatuR")
  arg_match(arg = species, values = c("human", "mouse"))
  SignatuR <- get("SignatuR", envir = asNamespace("SignatuR"))

  if (species == "human") {
    gene_list <- list(
      g1s = SignatuR$Hs$Programs$cellCycle.G1S,
      g2m = SignatuR$Hs$Programs$cellCycle.G2M,
      pseudogenes = SignatuR$Hs$Blocklists$Pseudogenes,
      noncoding = SignatuR$Hs$Blocklists$`Non-coding`,
      tcr = SignatuR$Hs$Compartments$TCR,
      immunoglobulins = SignatuR$Hs$Compartments$Immunoglobulins,
      ribo = SignatuR$Hs$Compartments$Ribo,
      mito = SignatuR$Hs$Compartments$Mito,
      heatshock = SignatuR$Hs$Programs$HeatShock
    )
  } else if (species == "mouse") {
    gene_list <- list(
      g1s = SignatuR$Mm$Programs$cellCycle.G1S,
      g2m = SignatuR$Mm$Programs$cellCycle.G2M,
      pseudogenes = SignatuR$Mm$Blocklists$Pseudogenes,
      noncoding = SignatuR$Mm$Blocklists$`Non-coding`,
      tcr = SignatuR$Mm$Compartments$TCR,
      immunoglobulins = SignatuR$Mm$Compartments$Immunoglobulins,
      ribo = SignatuR$Mm$Compartments$Ribo,
      mito = SignatuR$Mm$Compartments$Mito,
      heatshock = SignatuR$Mm$Programs$HeatShock
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
