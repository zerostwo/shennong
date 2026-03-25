set.seed(20260325)

cache_dir <- path.expand("~/.shennong/data")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

zenodo_record <- "14884845"
datasets <- c("pbmc1k", "pbmc3k")
filtered_sizes <- c(pbmc1k = 80, pbmc3k = 120)
extra_raw_sizes <- c(pbmc1k = 80, pbmc3k = 120)

fetch_h5 <- function(dataset, matrix_type, cache_dir, zenodo_record) {
  local_h5 <- file.path(cache_dir, paste0(dataset, "_", matrix_type, "_feature_bc_matrix.h5"))
  if (!file.exists(local_h5)) {
    remote_url <- glue::glue(
      "https://zenodo.org/records/{zenodo_record}/files/{dataset}_{matrix_type}_feature_bc_matrix.h5"
    )
    curl::curl_download(url = remote_url, destfile = local_h5)
  }
  local_h5
}

rename_with_prefix <- function(mat, prefix) {
  colnames(mat) <- paste(prefix, colnames(mat), sep = "_")
  mat
}

filtered_list <- list()
raw_list <- list()
meta_list <- list()

for (dataset in datasets) {
  filtered_path <- fetch_h5(dataset, "filtered", cache_dir, zenodo_record)
  raw_path <- fetch_h5(dataset, "raw", cache_dir, zenodo_record)

  filtered <- Seurat::Read10X_h5(filtered_path)
  raw <- Seurat::Read10X_h5(raw_path)

  keep_filtered <- sample(colnames(filtered), size = filtered_sizes[[dataset]])
  filtered_subset <- filtered[, keep_filtered, drop = FALSE]

  raw_extra_pool <- setdiff(colnames(raw), keep_filtered)
  keep_raw <- c(
    keep_filtered,
    sample(raw_extra_pool, size = extra_raw_sizes[[dataset]])
  )
  raw_subset <- raw[, keep_raw, drop = FALSE]

  filtered_list[[dataset]] <- rename_with_prefix(filtered_subset, dataset)
  raw_list[[dataset]] <- rename_with_prefix(raw_subset, dataset)
  meta_list[[dataset]] <- data.frame(
    cell = colnames(filtered_list[[dataset]]),
    sample = dataset,
    source_dataset = dataset,
    source_barcode = keep_filtered,
    stringsAsFactors = FALSE
  )
}

pbmc_small_counts <- do.call(cbind, filtered_list)
pbmc_small_raw <- do.call(cbind, raw_list)
pbmc_small_meta <- dplyr::bind_rows(meta_list)

pbmc_small <- SeuratObject::CreateSeuratObject(
  counts = pbmc_small_counts,
  meta.data = pbmc_small_meta,
  project = "pbmc_small"
)

save(
  pbmc_small,
  file = "data/pbmc_small.rda",
  compress = "bzip2"
)
save(
  pbmc_small_raw,
  file = "data/pbmc_small_raw.rda",
  compress = "bzip2"
)
