# Internal fallback signatures used when the optional SignatuR package is not
# available. These cover the categories required by Shennong core workflows.
.sn_fallback_signatures <- function(species = c("human", "mouse")) {
  species <- match.arg(species)

  human <- list(
    g1s = Seurat::cc.genes.updated.2019$s.genes,
    g2m = Seurat::cc.genes.updated.2019$g2m.genes,
    mito = c(
      "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
      "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5",
      "MT-ND6", "MT-RNR1", "MT-RNR2"
    ),
    ribo = c(
      "RPLP0", "RPLP1", "RPLP2", "RPL3", "RPL4", "RPL5", "RPL6", "RPL7",
      "RPL7A", "RPL8", "RPL9", "RPL10", "RPL10A", "RPL11", "RPL12",
      "RPL13", "RPL13A", "RPL14", "RPL15", "RPS3", "RPS4X", "RPS6",
      "RPS7", "RPS8", "RPS9", "RPS10", "RPS11", "RPS12", "RPS13",
      "RPS14", "RPS15", "RPS16", "RPS18", "RPS19"
    ),
    tcr = c("TRAC", "TRBC1", "TRBC2", "TRDC", "CD3D", "CD3E", "CD3G"),
    immunoglobulins = c("IGHM", "IGHD", "IGHG1", "IGKC", "IGLC2", "JCHAIN"),
    pseudogenes = character(),
    noncoding = character()
  )
  human$heatshock <- c("HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "DNAJB1", "HSPB1")

  mouse <- list(
    g1s = stringr::str_to_title(human$g1s),
    g2m = stringr::str_to_title(human$g2m),
    mito = stringr::str_replace(human$mito, "^MT-", "mt-"),
    ribo = stringr::str_to_title(human$ribo),
    tcr = stringr::str_to_title(human$tcr),
    immunoglobulins = stringr::str_to_title(human$immunoglobulins),
    pseudogenes = character(),
    noncoding = character(),
    heatshock = stringr::str_to_title(human$heatshock)
  )

  if (species == "human") human else mouse
}

.sn_supported_fallback_signature_categories <- function() {
  c("g1s", "g2m", "mito", "ribo", "heatshock", "tcr", "immunoglobulins", "pseudogenes", "noncoding")
}

#' Retrieve blocklist genes from SignatuR by preset gene set names
#'
#' This function uses the \code{SignatuR} package to fetch blocklist (or related)
#' genes (e.g. ribosomal, mitochondrial, heatshock, etc.). It then checks
#' those symbols via \code{HGNChelper} to remove invalid or outdated entries.
#'
#' @param species One of \code{"human"} or \code{"mouse"}.
#' @param category A character vector of predefined signature categories.
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
  arg_match(arg = species, values = c("human", "mouse"))
  category <- unique(category)

  use_signatur <- rlang::is_installed("SignatuR")
  if (use_signatur) {
    signature_env <- new.env(parent = emptyenv())
    utils::data("SignatuR", package = "SignatuR", envir = signature_env)
    use_signatur <- exists("SignatuR", envir = signature_env, inherits = FALSE)
  }

  if (use_signatur) {
    SignatuR <- get("SignatuR", envir = signature_env, inherits = FALSE)

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
    } else {
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
  } else {
    fallback_categories <- .sn_supported_fallback_signature_categories()
    unsupported <- setdiff(category, fallback_categories)
    if (length(unsupported) > 0) {
      stop(
        glue(
          "Package 'SignatuR' is required for signature categories: {paste(unsupported, collapse = ', ')}.\n",
          "Install it with:\n",
          "  remotes::install_github('carmonalab/SignatuR')"
        ),
        call. = FALSE
      )
    }
    gene_list <- .sn_fallback_signatures(species = species)
  }

  genes <- unlist(lapply(category, function(x) {
    if (!x %in% names(gene_list)) {
      warning("[g_block_genes] Unknown gene set: ", x)
      return(NULL)
    }
    if (use_signatur) {
      SignatuR::GetSignature(gene_list[[x]]) |> unlist()
    } else {
      unlist(gene_list[[x]], use.names = FALSE)
    }
  }))

  checked <- suppressWarnings(HGNChelper::checkGeneSymbols(genes, species = species))
  valid_genes <- checked$Suggested.Symbol[!is_na(checked$Suggested.Symbol)]

  if (length(valid_genes) < length(genes)) {
    invalid <- genes[is_na(checked$Suggested.Symbol)]
    message(
      "[g_block_genes] Removed ", length(invalid),
      " invalid symbols (e.g. ", paste(utils::head(invalid, 3), collapse = ", "), ")"
    )
  }

  return(unique(valid_genes))
}
