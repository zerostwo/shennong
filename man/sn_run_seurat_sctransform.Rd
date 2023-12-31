% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scrna-seq.R
\name{sn_run_seurat_sctransform}
\alias{sn_run_seurat_sctransform}
\title{SCTransform Workflow for Seurat}
\usage{
sn_run_seurat_sctransform(
  object,
  nfeatures = 3000,
  vars_to_regress = NULL,
  vst_flavor = "v2",
  seed = 717,
  verbose = TRUE
)
}
\arguments{
\item{object}{An object of class Seurat, representing single-cell RNA-seq
data.}

\item{nfeatures}{Use this many features as variable features after ranking by
residual variance; default is 3000. Only applied if residual.features is
not set.}

\item{vars_to_regress}{Variables to regress out in a second non-regularized
linear regression. For example, percent.mito. Default is NULL.}

\item{vst_flavor}{When set to 'v2' sets method = glmGamPoi_offset,
n_cells=2000, and exclude_poisson = TRUE which causes the model to learn
theta and intercept only besides excluding poisson genes from learning and
regularization.}

\item{seed}{Set a random seed. By default, sets the seed to 1448145. Setting
NULL will not set a seed.}

\item{verbose}{Whether to print messages and progress bars.}
}
\value{
Returns a Seurat object with a new assay (named SCT by default) with
  counts being (corrected) counts, data being log1p(counts), scale.data being
  pearson residuals; sctransform::vst intermediate results are saved in misc
  slot of the new assay.
}
\description{
Applies SCTransform workflow for single-cell RNA-seq data analysis. This
includes normalization and variance stabilization. Requires `glmGamPoi`
package for improved speed.
}
