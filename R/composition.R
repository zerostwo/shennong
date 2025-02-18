#' @export
sn_calculate_composition <- function(x, group_by, variable, min_cells = 20, additional_cols = NULL) {
  if (inherits(x, what = "Seurat")) {
    metadata <- x@meta.data
  } else {
    metadata <- x
  }

  if (min_cells > 0) {
    keep_group_bys <- metadata |>
      count(.data[[group_by]]) |>
      dplyr::filter(n >= min_cells) |>
      pull(.data[[group_by]])
    metadata <- metadata |>
      dplyr::filter(.data[[group_by]] %in% keep_group_bys)
  }
  adata <- table(metadata[[group_by]], metadata[[variable]]) |>
    prop.table(margin = 1) |>
    as.data.frame() |>
    mutate(Freq = Freq * 100) |>
    na.omit()
  colnames(x = adata) <- c(group_by, variable, "proportion")

  # Add additional columns to the output if specified
  if (!is.null(additional_cols)) {

  }

  return(adata)
}


#' calculate Startrac.dist (tissue distribution preference)
#' @import data.table
#' @importFrom plyr aaply
#' @importFrom stats chisq.test
#' @param dat.tb data.frame. Each line for a cell, and these columns as required: `majorCluster`, `loc`
#' @param byPatient logical. whether calculate the index for each patient. (default: FALSE)
#' @param colname.cluster character. which column specify the cluster (default: "majorCluster")
#' @param colname.patient character. which column specify the patient  (default: "patient")
#' @param colname.tissue character. which column specify the tissue  (default: "loc")
#' @param method character. method to use, one of "chisq", "fisher", and "freq"  (default: "chisq")
#' @param min.rowSum integer. rows with rowSum <= this value will be filtered out (default: 0)
#' @details calculate Startrac.dist (tissue distribution preference).
#' @return an array full of R_{o/e} (method="chisq") or list with components of OR, p.value etc. from fisher.test (method="fisher")
#' @export
sn_calculate_odds_ratio <- function(dat.tb,
                                    byPatient = F,
                                    colname.cluster = "majorCluster",
                                    colname.patient = "patient",
                                    colname.tissue = "loc",
                                    method = "chisq",
                                    min.rowSum = 0) {
  if (method == "freq") {
    ncount.sampleID <- dat.tb[, .(N.tot = .N), by = c(colname.patient, colname.tissue)]
    ncount.sampleID_mcls <- dat.tb[, .(N.mcls = .N), by = c(colname.patient, colname.tissue, colname.cluster)]
    ncount.sampleID_mcls <- reshape2::melt(
        reshape2::dcast(
            ncount.sampleID_mcls,
            as.formula(paste(colname.patient, "+", colname.tissue, "~", colname.cluster)),
            value.var = "N.mcls",
            fill = 0
        ),
        id.vars = c(colname.patient, colname.tissue),
        variable.name = colname.cluster,
        value.name = "N.mcls"
    )
    ncount.sampleID_mcls <- merge(ncount.sampleID_mcls, ncount.sampleID, by = c(colname.patient, colname.tissue))
    ncount.sampleID_mcls <- data.table::as.data.table(ncount.sampleID_mcls)
    ncount.sampleID_mcls[, freq.mcls := N.mcls / N.tot]
    res <- ncount.sampleID_mcls[,
      {
        loc.vec <- unique(.SD$loc)
        o.tb <- data.table::as.data.table(plyr::ldply(loc.vec, function(xx) {
          freq.x <- .SD[loc == xx, ][["freq.mcls"]]
          freq.y <- .SD[loc != xx, ][["freq.mcls"]]
          if (length(freq.x) >= 3) {
            oo.t <- t.test(freq.x, freq.y)
            oo.w <- wilcox.test(freq.x, freq.y)
            data.table(
              loc = xx,
              logFC = log2(mean(freq.x) / mean(freq.y)),
              p.value.t = oo.t$p.value,
              p.value.w = oo.w$p.value
            )
          } else {
            NULL
          }
        }))
      },
      by = c(colname.cluster)
    ][order(colname.cluster, colname.tissue), ]
    res[, FDR.t := p.adjust(p.value.t, "BH")]
    res[, FDR.w := p.adjust(p.value.w, "BH")]

    res[, char.sig := ""]
    res[FDR.t < 0.01 & p.value.t < 0.05, char.sig := "\U2020"]
    res[FDR.t < 0.05, char.sig := "\U2731"]
    res[FDR.t < 0.01, char.sig := "\U2731\U2731"]
    # res[FDR.t < 0.05,char.sig:="*"]
    # res[FDR.t < 0.01,char.sig:="**"]
    dist.charSig.tb <- dcast(res, majorCluster ~ loc, value.var = "char.sig")
    dist.logFC.tb <- dcast(res, majorCluster ~ loc, value.var = "logFC")
    dist.FDR.tb <- dcast(res, majorCluster ~ loc, value.var = "FDR.t")
    ret <- list(
      "res" = res,
      "dist.logFC.tb" = dist.logFC.tb,
      "dist.FDR.tb" = dist.FDR.tb,
      "dist.charSig.tb" = dist.charSig.tb
    )
  } else {
    if (byPatient == F) {
      N.o <- table(dat.tb[[colname.cluster]], dat.tb[[colname.tissue]])
      if (method == "chisq") {
        res.chisq <- chisq.test(N.o)
        R.oe <- (res.chisq$observed) / (res.chisq$expected)
        ret <- R.oe
      } else if (method == "fisher") {
        count.dist <- N.o[rowSums(N.o) > min.rowSum, , drop = F]
        count.dist.melt.ext.tb <- .table.fisher(count.dist)
        dist.p.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "p.value")
        dist.padj.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "p.adj")
        dist.OR.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "OR")
        ret <- list(
          "count.dist" = count.dist.melt.ext.tb,
          "p.tb" = dist.p.tb,
          "padj.tb" = dist.padj.tb,
          "OR.tb" = dist.OR.tb
        )
      }
    } else {
      N.o.byPatient <- table(
        dat.tb[[colname.patient]],
        dat.tb[[colname.cluster]], dat.tb[[colname.tissue]]
      )
      ret <- aaply(N.o.byPatient, 1, function(x) {
        if (method == "chisq") {
          res.chisq <- chisq.test(x)
          return((res.chisq$observed) / (res.chisq$expected))
        } else {
          res.fisher <- .table.fisher(x)
          return(dcast(res.fisher, rid ~ cid, value.var = "OR"))
        }
      })
    }
  }
  return(ret)
}

.table.fisher <- function(count.dist) {
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(unclass(count.dist))
  data.table::setDT(count.dist.tb, keep.rownames = T)
  count.dist.melt.tb <- reshape2::melt(count.dist.tb, id.vars = "rn")
  colnames(count.dist.melt.tb) <- c("rid", "cid", "count")
  count.dist.melt.ext.tb <- data.table::as.data.table(plyr::ldply(seq_len(nrow(count.dist.melt.tb)), function(i) {
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col] - this.c
    this.m <- matrix(
      c(
        this.c,
        sum.row[this.row] - this.c,
        other.col.c,
        sum(sum.col) - sum.row[this.row] - other.col.c
      ),
      ncol = 2
    )
    res.test <- fisher.test(this.m)
    data.frame(
      rid = this.row,
      cid = this.col,
      p.value = res.test$p.value,
      OR = res.test$estimate
    )
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb, count.dist.melt.ext.tb,
    by = c("rid", "cid")
  )
  count.dist.melt.ext.tb[, p.adj := p.adjust(p.value, "BH")]
  return(count.dist.melt.ext.tb)
}
