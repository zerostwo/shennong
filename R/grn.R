#' Run PySCENIC on a Seurat object
#'
#' This function runs PySCENIC on a Seurat object and generates a new Seurat
#' object with PySCENIC output.
#'
#' @param x A Seurat object
#' @param assay The name of the assay to use for PySCENIC analysis (default is
#'   "RNA")
#' @param slot The slot to use for PySCENIC analysis (default is "data")
#' @param outdir The output directory for PySCENIC analysis (default is
#'   "./pyscenic")
#' @param pyscenic The path to the PySCENIC executable (default is NULL, which
#'   uses the default path)
#' @param species The species to use for PySCENIC analysis (default is "human")
#' @param features A vector of gene names to use for PySCENIC analysis (default
#'   is NULL, which uses all genes)
#' @param tfs_fname The path to the file containing the list of transcription
#'   factors (default is NULL, which uses the default file)
#' @param database_fname The path to the file containing the database of motifs
#'   (default is NULL, which uses the default file)
#' @param annotations_fname The path to the file containing the annotations of
#'   motifs (default is NULL, which uses the default file)
#' @param mode The mode to use for PySCENIC analysis (default is
#'   "custom_multiprocessing")
#' @param ncores The number of cores to use for PySCENIC analysis (default is
#'   10)
#' @param seed The seed to use for PySCENIC analysis (default is 717)
#' @return A Seurat object with PySCENIC output
#' @importFrom stringr str_split_fixed
#' @importFrom tictoc tic toc
#' @importFrom logger log_info
#' @importFrom dplyr mutate select rename left_join count
#' @export
sn_run_pyscenic <-
  function(x,
           assay = "RNA",
           layer = "data",
           outdir = "../pyscenic",
           pyscenic = NULL,
           species = "human",
           features = NULL,
           tfs_fname = NULL,
           database_fname = NULL,
           annotations_fname = NULL,
           mode = "custom_multiprocessing",
           ncores = 10,
           seed = 717) {
    check_installed("data.table")
    if (is.null(pyscenic)) {
      pyscenic <- default_pyscenic_path
      sn_check_file(pyscenic)
    }
    if (species == "human") {
      tfs_fname <- "/mnt/resources/cistarget/tf_lists/allTFs_hg38.txt"
      database_fname <-
        c(
          "/mnt/resources/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
          "/mnt/resources/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
        )
      annotations_fname <-
        "/mnt/resources/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    } else if (species == "mouse") {
      tfs_fname <- "/mnt/resources/cistarget/tf_lists/allTFs_mm.txt"
      database_fname <-
        c(
          "/mnt/resources/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
          "/mnt/resources/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
        )
      annotations_fname <-
        "/mnt/resources/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
    } else {
      stop("Please check if the species is human or mouse!")
    }
    sn_check_file(x = c(tfs_fname, database_fname, annotations_fname))
    outdir <- sn_set_path(outdir)
    expression_mtx_fname <- file.path(outdir, "expression_mtx.csv")
    module_fname <-
      file.path(outdir, "expression_mtx.adjacencies.tsv")
    signatures_fname <- file.path(outdir, "regulons.gmt")
    auc_mtx_fname <- file.path(outdir, "auc_mtx.csv")
    tfs_target_fname <- file.path(outdir, "tfs_targer.tsv")
    auc_g_mtx_fname <- file.path(outdir, "auc_g_mtx.csv")
    # pre-processing
    if (!is.null(features)) {
      # x <- x[features,]
      tfs <- read.table(tfs_fname)[, 1]
      filtered_tfs <- intersect(tfs, features)
      tfs_fname <- file.path(outdir, "tfs.txt")
      write.table(
        x = filtered_tfs,
        file = tfs_fname,
        row.names = FALSE,
        quote = FALSE,
        col.names = FALSE
      )
      counts <-
        SeuratObject::LayerData(
          object = x,
          assay = assay,
          layer = layer
        )[features, ] |>
        Matrix::as.matrix() |>
        as.data.frame()
    } else {
      counts <- SeuratObject::LayerData(
        object = x,
        assay = assay,
        layer = layer
      ) |>
        Matrix::as.matrix() |>
        as.data.frame()
    }
    if (!file.exists(expression_mtx_fname)) {
      log_info("Write counts to a csv file...")
      tic("Write csv")
      data.table::fwrite(
        x = counts,
        file = expression_mtx_fname,
        quote = FALSE,
        row.names = TRUE
      )
      toc()
    }
    # grn
    if (!file.exists(... = module_fname)) {
      log_info("Derive co-expression modules from expression matrix...")
      tic("Run grn")
      system2(
        command = pyscenic,
        args = c(
          "grn",
          "--transpose",
          "--num_workers",
          ncores,
          "--seed",
          seed,
          "--output",
          module_fname,
          expression_mtx_fname,
          tfs_fname
        )
      )
      toc()
    }
    # ctx
    if (!file.exists(signatures_fname)) {
      log_info(
        "Find enriched motifs for a gene signature and optionally prune targets from this signature based on cis-regulatory cues..."
      )
      tic("Run ctx")
      system2(
        command = pyscenic,
        args = c(
          "ctx",
          "--transpose",
          "--output",
          signatures_fname,
          "--mode",
          mode,
          "--annotations_fname",
          annotations_fname,
          "--num_workers",
          ncores,
          "--expression_mtx_fname",
          expression_mtx_fname,
          module_fname,
          database_fname
        )
      )
      toc()
    }
    # aucell
    if (!file.exists(auc_mtx_fname)) {
      log_info("Quantify activity of gene signatures across single cells...")
      tic("Run aucell")
      system2(
        command = pyscenic,
        args = c(
          "aucell",
          "--transpose",
          "--output",
          auc_mtx_fname,
          "--num_workers",
          ncores,
          expression_mtx_fname,
          signatures_fname
        )
      )
      toc()
    }

    # tfs_target
    if (!file.exists(tfs_target_fname)) {
      signatures <-
        clusterProfiler::read.gmt(signatures_fname)
      signatures_count <- signatures |>
        count(term)
      signatures <-
        left_join(
          x = signatures, y = signatures_count,
          by = "term"
        )
      signatures <- signatures |>
        mutate(term = str_split_fixed(
          string = term,
          pattern = "\\(",
          n = 2
        )[, 1]) |>
        mutate(symbol = paste0(term, " (", n, "g)")) |>
        select(symbol, term, gene) |>
        rename(target_gene = gene, tf = term)
      write.csv(signatures,
        tfs_target_fname,
        row.names = FALSE
      )
    } else {
      signatures <- read.csv(tfs_target_fname)
    }
    # # auc_g_mtx
    # if (!file.exists(auc_g_mtx_fname)) {
    #   auc_mtx <-
    #     read.csv(auc_mtx_fname,
    #              row.names = 1,
    #              check.names = FALSE)
    #   rownames(auc_mtx) <- unique(signatures$tf)
    #   fwrite(
    #     x = auc_mtx,
    #     file = auc_g_mtx_fname,
    #     quote = FALSE,
    #     row.names = TRUE
    #   )
    # } else {
    #   auc_mtx <- read.csv(auc_g_mtx_fname,
    #                       row.names = 1,
    #                       check.names = FALSE)
    # }
    # # colnames(auc_mtx) <-
    # #   gsub(pattern = "\\.",
    # #        replacement = "-",
    # #        x = colnames(auc_mtx))
    # x[["scenic"]] <- CreateAssayObject(counts = auc_mtx)
    #
    # return(x)
  }
