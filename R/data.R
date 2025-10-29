#' Signature genes
"SignatuR"

#' Load PBMC datasets from Zenodo
#'
#' This function downloads (if not already cached) and loads processed
#' single-cell RNA-seq datasets of peripheral blood mononuclear cells (PBMCs)
#' from a Zenodo record. Each dataset (1k, 3k, 4k, or 8k) includes both
#' \emph{filtered} and \emph{raw} gene–barcode matrices in standard Cell Ranger
#' HDF5 format.
#'
#' Instead of storing serialized R objects, the function caches the original
#' Zenodo \code{.h5} files locally for reproducibility and version consistency.
#' When requested, it dynamically constructs and returns either:
#' \itemize{
#'   \item a Seurat object (for filtered data), or
#'   \item a sparse count matrix (for raw data).
#' }
#' Users can also choose to only download/cache the data without loading it
#' into memory.
#'
#' @param dataset Character scalar. Which PBMC dataset to load.
#'   One of \code{"pbmc1k"}, \code{"pbmc3k"}, \code{"pbmc4k"}, \code{"pbmc8k"}.
#'   Default: \code{"pbmc3k"}.
#'
#' @param matrix_type Character scalar. Which matrix type to load.
#'   One of:
#'   \itemize{
#'     \item \code{"filtered"}: the \code{filtered_feature_bc_matrix.h5}
#'           (Cell Ranger filtered barcodes).
#'     \item \code{"raw"}: the \code{raw_feature_bc_matrix.h5}
#'           (unfiltered barcodes; useful for ambient RNA correction tools
#'           such as SoupX).
#'   }
#'   Default: \code{"filtered"}.
#'
#' @param save_dir Character scalar. Local cache directory where the downloaded
#'   Zenodo files will be stored. The directory will be created if it does not
#'   exist. Default: \code{"~/.shennong/data"}.
#'
#' @param return_object Logical. If \code{TRUE} (default), the function returns
#'   an in-memory object. If \code{FALSE}, the function only ensures that the
#'   file is downloaded locally and then returns (invisibly) the local file
#'   path.
#'
#' @param species Character scalar. Species label passed to
#'   \code{sn_initialize_seurat_object()} when constructing Seurat objects
#'   (i.e. when \code{matrix_type == "filtered"}). Ignored if
#'   \code{matrix_type == "raw"}. Default: \code{"human"}.
#'
#' @details
#' All datasets were re-aligned using Cell Ranger v9.0.1 with a custom reference
#' based on GENCODE v48 (GRCh38.p14).
#'
#' \strong{Caching strategy}
#'
#' The function caches the original files from Zenodo on first use:
#'
#' \preformatted{
#' {dataset}_{matrix_type}_feature_bc_matrix.h5
#' }
#'
#' For example, after loading \code{pbmc3k} you may see:
#'
#' \preformatted{
#' ~/.shennong/data/
#' ├─ pbmc1k_filtered_feature_bc_matrix.h5
#' ├─ pbmc1k_raw_feature_bc_matrix.h5
#' ├─ pbmc3k_filtered_feature_bc_matrix.h5
#' ├─ pbmc3k_raw_feature_bc_matrix.h5
#' ├─ pbmc4k_filtered_feature_bc_matrix.h5
#' ├─ pbmc4k_raw_feature_bc_matrix.h5
#' ├─ pbmc8k_filtered_feature_bc_matrix.h5
#' └─ pbmc8k_raw_feature_bc_matrix.h5
#' }
#'
#' When \code{return_object = TRUE}, the cached \code{.h5} file is read on the
#' fly:
#' \itemize{
#'   \item If \code{matrix_type == "filtered"}:
#'         a Seurat object is constructed via
#'         \code{sn_initialize_seurat_object()}.
#'   \item If \code{matrix_type == "raw"}:
#'         the function returns a sparse count matrix
#'         (typically a \code{dgCMatrix}), suitable for ambient RNA
#'         estimation / SoupX workflows.
#' }
#'
#' When \code{return_object = FALSE}, no object is constructed;
#' only the file is ensured to exist locally.
#'
#' @return
#' One of:
#'
#' \itemize{
#'   \item If \code{matrix_type == "filtered"} and
#'         \code{return_object = TRUE}:
#'         a Seurat object.
#'
#'   \item If \code{matrix_type == "raw"} and
#'         \code{return_object = TRUE}:
#'         a sparse count matrix (e.g. \code{dgCMatrix}).
#'
#'   \item If \code{return_object = FALSE}:
#'         the local file path to the cached \code{.h5} file,
#'         returned invisibly.
#' }
#'
#' @examples
#' \dontrun{
#' # 1. Load filtered PBMC3k as a Seurat object:
#' pbmc <- sn_load_pbmc()
#'
#' # 2. Load raw PBMC3k counts as a sparse matrix (for SoupX etc.):
#' pbmc_raw <- sn_load_pbmc(matrix_type = "raw")
#'
#' # 3. Only download/cache PBMC8k, don't construct anything in-memory:
#' sn_load_pbmc(dataset = "pbmc8k", return_object = FALSE)
#'
#' # 4. Use a custom cache directory:
#' pbmc4k <- sn_load_pbmc(
#'   dataset  = "pbmc4k",
#'   save_dir = "~/datasets/pbmc_cache"
#' )
#' }
#'
#' @references
#' Duan S. (2025).
#' Processed PBMC datasets re-aligned with Cell Ranger v9.0.1
#' (GENCODE v48, GRCh38.p14). Zenodo.
#' DOI: 10.5281/zenodo.14884845
#'
#' @seealso
#' \code{\link{sn_initialize_seurat_object}},
#' \code{\link{sn_read_matrix_h5}},
#' \code{\link{sn_write}},
#' \code{\link{sn_read}}
#'
#' @export
sn_load_pbmc <- function(dataset = "pbmc3k",
                         matrix_type = c("filtered", "raw"),
                         save_dir = "~/.shennong/data",
                         return_object = TRUE,
                         species = "human") {
  matrix_type <- match.arg(matrix_type)
  zenodo_record <- "14884845"

  save_dir <- sn_set_path(save_dir)

  local_h5 <- glue("{save_dir}/{dataset}_{matrix_type}_feature_bc_matrix.h5")

  if (file.exists(local_h5)) {
    counts <- sn_read(local_h5)
  } else {
    remote_url <- glue::glue(
      "https://zenodo.org/records/{zenodo_record}/files/{dataset}_{matrix_type}_feature_bc_matrix.h5"
    )
    cli::cli_inform("✨ Fetching {dataset} ({matrix_type}) from Zenodo cache...")
    curl::curl_download(url = remote_url, destfile = local_h5)
  }

  if (matrix_type == "filtered") {
    x <- sn_initialize_seurat_object(
      x = counts,
      species = species
    )
  } else {
    x <- counts
  }

  return(x)
}


#' Pan-Immune CellTypist Metadata
#'
#' This dataset contains hierarchical immune cell type metadata from the Pan-Immune CellTypist atlas.
#' It includes high- and low-hierarchy cell type annotations, human and mouse marker genes, and curated marker gene classification.
#'
#' **Reference**:
#' C. Domínguez Conde et al., Science, 2022. DOI: [10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197)
#'
#' **Original Data Source**:
#' The dataset was obtained from the CellTypist repository at:
#' [GitHub: Teichlab/celltypist_wiki](https://github.com/Teichlab/celltypist_wiki/blob/main/atlases/Pan_Immune_CellTypist/v2/tables/Basic_celltype_information.xlsx)
#'
#' **Description**:
#' A dataset providing hierarchical immune cell types with their ontology identifiers, curated human and mouse marker genes, and a classification flag for high-hierarchy marker genes.
#'
#' @format A data frame with 350 rows and 5 columns:
#' \describe{
#'   \item{high_hierarchy_cell_type}{Character. Broad cell type category (e.g., "B cells").}
#'   \item{low_hierarchy_cell_type}{Character. More specific cell type classification (e.g., "Follicular B cells").}
#'   \item{human}{Character. Human marker gene associated with the cell type (e.g., "CD79A").}
#'   \item{mouse}{Character. Corresponding mouse homolog of the human marker gene (e.g., "Cd79a").}
#'   \item{high_hierarchy_marker_gene}{Logical. Indicates whether the marker gene is a high-hierarchy marker gene (TRUE/FALSE).}
#' }
#'
#' @source CellTypist Atlas: Pan-Immune dataset, Teichlab
#' [GitHub Repository](https://github.com/Teichlab/celltypist_wiki)
#'
#' @usage data(pan_immune_celltypist)
#'
"marker_genes"

#' Human-Mouse Homologous Genes Table
#'
#' This dataset contains homologous gene mappings between human and mouse,
#' including gene symbols, HGNC IDs, and MGI IDs.
#'
#' The data was sourced from the Mouse Genome Informatics (MGI) Homology report:
#' **HOM_MouseHumanSequence.rpt**.
#'
#' **Reference Databases**:
#' - **MGI (Mouse Genome Informatics)**: <http://www.informatics.jax.org/>
#' - **Original Data Source**: [HOM_MouseHumanSequence.rpt](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt)
#' - **HGNC (HUGO Gene Nomenclature Committee)**: <https://www.genenames.org/>
#' - **Ensembl Biomart**: <https://www.ensembl.org/biomart>
#'
#' @format A data frame with 5 columns:
#' \describe{
#'   \item{DB.Class.Key}{Integer. A unique key for each gene homology group in the MGI database.}
#'   \item{human}{Character. Approved human gene symbol from HGNC.}
#'   \item{HGNC.ID}{Character. The corresponding HGNC ID for the human gene.}
#'   \item{mouse}{Character. Approved mouse gene symbol from MGI.}
#'   \item{Mouse.MGI.ID}{Character. The corresponding MGI ID for the mouse gene.}
#' }
#'
#' @source MGI Homology Report: [HOM_MouseHumanSequence.rpt](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt)
#'
#' @usage data(hom_genes)
#'
"hom_genes"
