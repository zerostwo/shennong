# Score Cell Cycle Phases

This function scores the cell cycle phase of single-cell RNA-seq data
using S and G2M phase marker genes.

## Usage

``` r
sn_score_cell_cycle(object, species = NULL)
```

## Arguments

- object:

  A Seurat object containing single-cell RNA-seq data.

- species:

  (Optional) A character string indicating the species (e.g., "human" or
  "mouse"). If NULL, the function will attempt to retrieve species
  information from `Seurat::Misc(object)`.

## Value

A Seurat object with cell cycle scores added, including: - `S.Score`:
Score for the S phase - `G2M.Score`: Score for the G2M phase - `Phase`:
Assigned cell cycle phase - `CC.Difference`: Difference between S and
G2M scores

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(rpois(400 * 20, lambda = 3), nrow = 400, ncol = 20)
  rownames(counts) <- c(
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2",
    "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS",
    "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN",
    "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3",
    "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45",
    "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
    "BRIP1", "E2F8",
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5",
    "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
    "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3",
    "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
    "KIF20B", "HJURP", "CDCA3", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5",
    "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3",
    "CBX5", "CENPA",
    paste0("GENE", 1:316)
  )
  colnames(counts) <- paste0("cell", 1:20)
  obj <- sn_initialize_seurat_object(counts, species = "human")
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- sn_score_cell_cycle(obj, species = "human")
  head(obj[[]][, c("S.Score", "G2M.Score", "Phase")])
}
} # }
```
