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
