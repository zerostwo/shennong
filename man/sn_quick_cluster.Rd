% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integrate.R
\name{sn_quick_cluster}
\alias{sn_quick_cluster}
\title{Quickly run a clustering pipeline (standard or SCTransform-based)}
\usage{
sn_quick_cluster(
  object,
  pipeline = "standard",
  nfeatures = 3000,
  vars_to_regress = NULL,
  resolution = 0.8,
  dims = 1:30,
  return_cluster = FALSE,
  verbose = TRUE,
  species = NULL
)
}
\arguments{
\item{object}{A \code{Seurat} object.}

\item{pipeline}{One of \code{"standard"} or \code{"sctransform"} (alias:
\code{"sct"}).}

\item{nfeatures}{Number of variable features to find.}

\item{vars_to_regress}{Covariates for regression in \code{ScaleData} or
\code{SCTransform}.}

\item{resolution}{Resolution parameter for \code{FindClusters}.}

\item{dims}{Dimensions to use for PCA/UMAP/FindNeighbors.}

\item{return_cluster}{If \code{TRUE}, returns just the cluster vector instead
of the entire Seurat object.}

\item{verbose}{Whether to print/log progress messages.}
}
\value{
A \code{Seurat} object with clustering and UMAP performed, or a
  cluster vector if \code{return_cluster=TRUE}.
}
\description{
This function provides a simplified clustering workflow for a Seurat object,
including data normalization, variable feature selection, scaling, PCA,
neighbor finding, clustering, and optional UMAP embedding.
}
\examples{
\dontrun{
seurat_obj <- sn_quick_cluster(seurat_obj, pipeline = "standard", resolution = 0.8)
}
}
