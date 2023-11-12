#' Run Seurat Workflow
#'
#' Executes a Seurat workflow for single-cell RNA-seq data analysis. Supports
#' both standard and SCTransform workflows, followed by PCA reduction and
#' post-processing.
#'
#' @inheritParams sn_run_seurat_sctransform
#' @inheritParams sn_run_seurat_postprocessing
#' @param workflow The workflow to use. Can be either "standard" or
#'   "sctransform".
#'
#' @return An updated Seurat object after performing the specified workflow and
#'   post-processing.
#' @importFrom rlang arg_match
#' @importFrom Seurat RunPCA
#' @export
sn_run_seurat <- function(object,
                          workflow = "sctransform",
                          nfeatures = 3000,
                          vars_to_regress = NULL,
                          vst_flavor = "v2",
                          reduction = "pca",
                          dims = 1:30,
                          resolution = 0.8,
                          seed = 717,
                          verbose = TRUE,
                          run_tsne = FALSE) {
  arg_match(workflow, values = c("standard", "sctransform"))

  if (workflow == "standard") {
    object <- sn_run_seurat_standard(
      object = object,
      nfeatures = nfeatures,
      vars_to_regress = vars_to_regress,
      verbose = verbose
    )
  } else {
    object <- sn_run_seurat_sctransform(
      object = object,
      nfeatures = nfeatures,
      vars_to_regress = vars_to_regress,
      vst_flavor = vst_flavor,
      seed = seed,
      verbose = verbose
    )
  }

  # Perform linear dimensional reduction
  object <- RunPCA(
    object = object,
    seed.use = seed,
    verbose = verbose
  )

  # Perform Seurat post processing
  object <- sn_run_seurat_postprocessing(
    object = object,
    reduction = reduction,
    dims = dims,
    resolution = resolution,
    seed = seed,
    verbose = verbose,
    run_tsne = run_tsne
  )

  return(object)
}

#' Standard Seurat Workflow
#'
#' Performs the standard workflow for Seurat single-cell RNA-seq data analysis.
#' This includes normalization, identification of variable features, and
#' scaling.
#'
#' @param object An object of class Seurat, representing single-cell RNA-seq
#'   data.
#' @param nfeatures Use this many features as variable features after ranking by
#'   residual variance; default is 3000. Only applied if residual.features is
#'   not set.
#' @param vars_to_regress Variables to regress out in a second non-regularized
#'   linear regression. For example, percent.mito. Default is NULL.
#' @param verbose Whether to print messages and progress bars.
#'
#' @return Seurat object after applying the standard workflow.
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData
#' @export
sn_run_seurat_standard <- function(object,
                                   nfeatures = 2000,
                                   vars_to_regress = NULL,
                                   verbose = TRUE) {
  object <- NormalizeData(object = object, verbose = verbose)
  object <- FindVariableFeatures(
    object = object,
    nfeatures = nfeatures,
    verbose = verbose
  )
  object <- ScaleData(
    object = object,
    vars.to.regress = vars_to_regress,
    verbose = verbose
  )
  return(object)
}


#' SCTransform Workflow for Seurat
#'
#' Applies SCTransform workflow for single-cell RNA-seq data analysis. This
#' includes normalization and variance stabilization. Requires `glmGamPoi`
#' package for improved speed.
#'
#' @inheritParams sn_run_seurat_standard
#' @param vst_flavor When set to 'v2' sets method = glmGamPoi_offset,
#'   n_cells=2000, and exclude_poisson = TRUE which causes the model to learn
#'   theta and intercept only besides excluding poisson genes from learning and
#'   regularization.
#' @param seed Set a random seed. By default, sets the seed to 1448145. Setting
#'   NULL will not set a seed.
#'
#' @return Returns a Seurat object with a new assay (named SCT by default) with
#'   counts being (corrected) counts, data being log1p(counts), scale.data being
#'   pearson residuals; sctransform::vst intermediate results are saved in misc
#'   slot of the new assay.
#' @importFrom Seurat SCTransform
#' @export
sn_run_seurat_sctransform <- function(object,
                                      nfeatures = 3000,
                                      vars_to_regress = NULL,
                                      vst_flavor = "v2",
                                      seed = 717,
                                      verbose = TRUE) {
  # Check if glmGamPoi package is installed
  sn_check_package("glmGamPoi", error = FALSE)

  object <- SCTransform(
    object = object,
    variable.features.n = nfeatures,
    vars.to.regress = vars_to_regress,
    vst.flavor = vst_flavor,
    seed.use = seed,
    verbose = verbose
  )
  return(object)
}

#' Post-Processing for Seurat Workflow
#'
#' Executes post-processing steps for Seurat workflow including neighbor
#' finding, clustering, and UMAP/t-SNE reduction.
#'
#' @inheritParams sn_run_seurat_sctransform
#' @param reduction Reduction to use as input for building the (S)NN.
#' @param dims Dimensions of reduction to use as input.
#' @param resolution Value of the resolution parameter, use a value above
#'   (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param run_tsne Whether to run t-SNE reduction.
#'
#' @return Seurat object after post-processing.
#' @importFrom Seurat FindNeighbors FindClusters RunUMAP
#' @export
sn_run_seurat_postprocessing <- function(
    object,
    reduction = "pca",
    dims = 1:30,
    resolution = 0.8,
    seed = 717,
    verbose = TRUE,
    run_tsne = FALSE) {
  object <- FindNeighbors(
    object = object,
    reduction = reduction,
    dims = dims,
    verbose = verbose
  )
  object <- FindClusters(
    object = object,
    resolution = resolution,
    random.seed = seed,
    verbose = verbose
  )
  object <- RunUMAP(
    object = object,
    dims = dims,
    reduction = reduction,
    seed.use = seed,
    verbose = verbose
  )
  if (run_tsne) {
    object <- RunTSNE(
      object = object,
      dims = dims,
      reduction = reduction,
      seed.use = seed
    )
  }
  return(object)
}
