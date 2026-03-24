# Build bundled Shennong gene annotations from local GENCODE GTF files.
# Run manually during development.

human_gtf <- Sys.getenv(
  "SHENNONG_GTF_HUMAN",
  "/mnt/reference_genomes/gencode/human/48/GRCh38.v48.primary_assembly.annotation.gtf"
)
mouse_gtf <- Sys.getenv(
  "SHENNONG_GTF_MOUSE",
  "/mnt/reference_genomes/gencode/mouse/M37/GRCm39.vM37.primary_assembly.annotation.gtf"
)
data_path <- file.path("data", "shennong_gene_annotations.rda")

extract_gtf_attribute <- function(x, key) {
  pattern <- sprintf('(?:^|; )%s "([^"]+)"', key)
  matches <- regexec(pattern, x, perl = TRUE)
  captures <- regmatches(x, matches)
  vapply(
    captures,
    function(hit) {
      if (length(hit) >= 2) hit[[2]] else NA_character_
    },
    character(1)
  )
}

classify_gene_type <- function(gene_type) {
  coding_types <- c(
    "protein_coding",
    "IG_C_gene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_LV_gene",
    "IG_V_gene",
    "TR_C_gene",
    "TR_D_gene",
    "TR_J_gene",
    "TR_V_gene"
  )

  out <- ifelse(
    is.na(gene_type) | !nzchar(gene_type),
    NA_character_,
    ifelse(gene_type %in% coding_types, "coding", "noncoding")
  )
  unname(out)
}

read_gene_rows <- function(path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install data.table before building gene annotations.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("GTF file does not exist: %s", path), call. = FALSE)
  }

  cmd <- sprintf("grep -v '^#' %s | awk '$3 == \"gene\"'", shQuote(path))
  data.table::fread(
    cmd = cmd,
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    data.table = FALSE,
    showProgress = TRUE
  )
}

parse_gtf_genes <- function(path, species, release) {
  raw <- read_gene_rows(path)
  colnames(raw) <- c(
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes"
  )

  gene_id <- extract_gtf_attribute(raw$attributes, "gene_id")
  gene_name <- extract_gtf_attribute(raw$attributes, "gene_name")
  gene_type <- extract_gtf_attribute(raw$attributes, "gene_type")
  gene_biotype <- extract_gtf_attribute(raw$attributes, "gene_biotype")
  gene_status <- extract_gtf_attribute(raw$attributes, "gene_status")
  gene_source <- extract_gtf_attribute(raw$attributes, "gene_source")
  level <- extract_gtf_attribute(raw$attributes, "level")
  resolved_gene_type <- ifelse(is.na(gene_type) | !nzchar(gene_type), gene_biotype, gene_type)

  out <- data.frame(
    species = species,
    release = release,
    seqname = raw$seqname,
    source = raw$source,
    start = as.integer(raw$start),
    end = as.integer(raw$end),
    strand = raw$strand,
    gene_id = gene_id,
    gene_id_base = sub("\\..*$", "", gene_id),
    gene_name = gene_name,
    gene_type = resolved_gene_type,
    gene_status = gene_status,
    gene_source = gene_source,
    level = level,
    gene_class = classify_gene_type(resolved_gene_type),
    stringsAsFactors = FALSE
  )

  out[!duplicated(out[c("species", "gene_id", "gene_name")]), , drop = FALSE]
}

shennong_gene_annotations <- rbind(
  parse_gtf_genes(human_gtf, species = "human", release = "GENCODE v48"),
  parse_gtf_genes(mouse_gtf, species = "mouse", release = "GENCODE vM37")
)

shennong_gene_annotations <- shennong_gene_annotations[
  order(
    shennong_gene_annotations$species,
    shennong_gene_annotations$seqname,
    shennong_gene_annotations$start,
    shennong_gene_annotations$gene_name,
    shennong_gene_annotations$gene_id
  ),
  ,
  drop = FALSE
]

dir.create(dirname(data_path), recursive = TRUE, showWarnings = FALSE)
save(
  shennong_gene_annotations,
  file = data_path,
  version = 2,
  compress = "xz"
)
