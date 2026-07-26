#' Pan-Immune CellTypist Metadata
#'
#' Hierarchical immune-cell annotations and curated marker genes from the
#' Pan-Immune CellTypist atlas.
#'
#' @format A data frame with 350 rows and 5 columns:
#' \describe{
#'   \item{\code{high_hierarchy_cell_type}}{Broad cell type category.}
#'   \item{\code{low_hierarchy_cell_type}}{More specific cell type label.}
#'   \item{\code{human}}{Human marker gene.}
#'   \item{\code{mouse}}{Mouse homolog marker gene.}
#'   \item{\code{high_hierarchy_marker_gene}}{Whether the marker is considered
#'   high-hierarchy in the source atlas.}
#' }
#'
#' @source CellTypist Pan-Immune atlas table:
#'   \url{https://github.com/Teichlab/celltypist_wiki/blob/main/atlases/Pan_Immune_CellTypist/v2/tables/Basic_celltype_information.xlsx}
"marker_genes"

#' Human-Mouse Homologous Genes Table
#'
#' Human and mouse homolog mappings derived from the MGI homology report.
#'
#' @format A data frame with 5 columns:
#' \describe{
#'   \item{\code{DB.Class.Key}}{MGI homology-group identifier.}
#'   \item{\code{human}}{Approved human gene symbol.}
#'   \item{\code{HGNC.ID}}{HGNC identifier for the human gene.}
#'   \item{\code{mouse}}{Approved mouse gene symbol.}
#'   \item{\code{Mouse.MGI.ID}}{MGI identifier for the mouse gene.}
#' }
#'
#' @source MGI homology report:
#'   \url{https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt}
"hom_genes"

#' Bundled Shennong Signature Catalog
#'
#' A package-owned, tree-structured signature catalog built from the upstream
#' \pkg{SignatuR} hierarchy and stored as the current Shennong-controlled
#' snapshot under \code{data/}.
#'
#' @format A named list with two top-level entries:
#' \describe{
#'   \item{\code{metadata}}{Build metadata such as source package, source
#'   version, and build date.}
#'   \item{\code{tree}}{A nested tree for \code{human} and \code{mouse}, where
#'   each node records its \code{name}, \code{kind}, optional \code{genes}, and
#'   optional \code{children}.}
#' }
#'
#' @source Imported from \pkg{SignatuR} during development and built into
#'   package data with \code{data-raw/build_shennong_signatures.R}.
"shennong_signature_catalog"

#' Bundled Shennong Gene Annotations
#'
#' A package-owned snapshot of gene-level annotations parsed from the GENCODE
#' human GRCh38 v48 and mouse GRCm39 vM37 primary-assembly GTF files.
#'
#' @format A data frame with one row per gene and the following columns:
#' \describe{
#'   \item{\code{species}}{Species label used by Shennong.}
#'   \item{\code{release}}{GENCODE release identifier.}
#'   \item{\code{seqname}}{Sequence or chromosome name from the GTF.}
#'   \item{\code{source}}{Annotation source from the GTF.}
#'   \item{\code{start}}{1-based gene start coordinate.}
#'   \item{\code{end}}{1-based gene end coordinate.}
#'   \item{\code{strand}}{Gene strand.}
#'   \item{\code{gene_id}}{Original GENCODE gene identifier, including version
#'   suffix when present.}
#'   \item{\code{gene_id_base}}{Version-stripped GENCODE gene identifier.}
#'   \item{\code{gene_name}}{Gene symbol from the GTF.}
#'   \item{\code{gene_type}}{GENCODE gene biotype.}
#'   \item{\code{gene_status}}{Optional GTF \code{gene_status} attribute.}
#'   \item{\code{gene_source}}{Optional GTF \code{gene_source} attribute.}
#'   \item{\code{level}}{Optional GENCODE annotation level.}
#'   \item{\code{gene_class}}{Shennong coarse classification used by gene
#'   filtering helpers.}
#' }
#'
#' @source Built from local GENCODE sources via
#'   \code{data-raw/build_shennong_gene_annotations.R}.
"shennong_gene_annotations"
