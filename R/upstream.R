#' Run fastp for a sample (single-end or paired-end)
#'
#' @param sample_id Sample ID used as output prefix. If NULL, inferred from input filename.
#' @param reads Character vector of one (single-end) or two (paired-end) FASTQ file paths.
#' @param outdir Output directory. Default is "./results".
#' @param ncores Number of threads to use. Default is 4.
#' @param execute Whether to run the command (TRUE) or return it as a string (FALSE).
#' @param fastp_path Path to the fastp executable.
#' @param overwrite Whether to overwrite existing output files. Default is FALSE.
#' @param extra_args Additional arguments to pass directly to fastp (e.g. "--trim_front1 5").
#'
#' @return If `execute = TRUE`, returns (invisibly) the output file paths; otherwise, the command string.
#' @export
run_fastp <- function(
  sample_id = NULL,
  reads,
  outdir = "./results",
  ncores = 4,
  execute = FALSE,
  fastp_path = "/opt/fastp/0.24.0/bin/fastp",
  overwrite = FALSE,
  extra_args = ""
) {
  stopifnot(length(reads) %in% c(1, 2))

  is_paired <- length(reads) == 2
  is_gzipped <- grepl("\\.gz$", reads[1])
  suffix <- if (is_gzipped) ".gz" else ""

  # Auto sample_id if not provided
  if (is_null(sample_id) || sample_id == "") {
    sample_id <- sub("(_R[12])?(_[12])?\\.(fq|fastq)(\\.gz)?$", "", basename(reads[1]))
  }

  # Create output directory
  outdir_fastp <- file.path(outdir, "fastp")
  if (!dir.exists(outdir_fastp) && execute) {
    dir.create(outdir_fastp, recursive = TRUE)
  }

  # Construct output file paths
  if (is_paired) {
    output1 <- file.path(outdir_fastp, paste0(sample_id, "_1.fastq", suffix))
    output2 <- file.path(outdir_fastp, paste0(sample_id, "_2.fastq", suffix))
    json_file <- file.path(outdir_fastp, paste0(sample_id, ".fastp.json"))
    html_file <- file.path(outdir_fastp, paste0(sample_id, ".fastp.html"))

    cmd <- glue(
      "{fastp_path} ",
      "--in1 {reads[1]} ",
      "--in2 {reads[2]} ",
      "--out1 {output1} ",
      "--out2 {output2} ",
      "--json {json_file} ",
      "--html {html_file} ",
      "--thread {ncores} ",
      "{extra_args}"
    )

    out_files <- c(output1, output2, json_file, html_file)
  } else {
    output <- file.path(outdir_fastp, paste0(sample_id, ".fastq", suffix))
    json_file <- file.path(outdir_fastp, paste0(sample_id, ".fastp.json"))
    html_file <- file.path(outdir_fastp, paste0(sample_id, ".fastp.html"))

    cmd <- glue(
      "{fastp_path} ",
      "--in1 {reads[1]} ",
      "--out1 {output} ",
      "--json {json_file} ",
      "--html {html_file} ",
      "--thread {ncores} ",
      "{extra_args}"
    )

    out_files <- c(output, json_file, html_file)
  }

  # Execute or return command
  if (execute) {
    existing <- out_files[file.exists(out_files)]
    if (length(existing) > 0 && !overwrite) {
      message("Output files already exist and overwrite = FALSE: ", paste(existing, collapse = ", "))
      return(invisible(NULL))
    }

    message("Running fastp for sample: ", sample_id)
    system(cmd)
    message("Finished processing sample: ", sample_id)
    return(invisible(out_files))
  } else {
    return(as.character(cmd))
  }
}

#' Filter rRNA reads using HISAT2
#'
#' @param sample_id Sample ID used for output naming. If NULL, inferred from input reads.
#' @param reads A character vector of 1 (single-end) or 2 (paired-end) fastq files.
#' @param rRNA_index Path to rRNA index base name (without .ht2 extension).
#' @param outdir Output directory for filtered fastq and summary.
#' @param ncores Number of threads to use.
#' @param hisat2_path Path to the hisat2 executable.
#' @param samtools_path Path to the samtools executable.
#' @param execute Whether to execute the command or return it as a string.
#' @param overwrite Whether to overwrite existing output.
#'
#' @return Filtered fastq files and summary file (or command string).
#' @export
filter_rRNA <- function(
  sample_id = NULL,
  reads,
  rRNA_index = "/mnt/reference_genomes/gencode/mouse/M36/index/hisat2/rRNA/genome",
  outdir = "./results/alignment/rRNA_dup",
  ncores = 4,
  hisat2_path = "/opt/hisat2/2.2.1/bin/hisat2",
  samtools_path = "/opt/samtools/1.21/bin/samtools",
  execute = FALSE,
  overwrite = FALSE
) {
  stopifnot(length(reads) %in% c(1, 2))
  is_paired <- length(reads) == 2
  is_gz <- grepl("\\.gz$", reads[1])
  suffix <- if (is_gz) ".gz" else ""

  # 自动推断 sample_id
  if (is_null(sample_id) || sample_id == "") {
    sample_id <- sub("(_R[12])?(_[12])?\\.(fq|fastq)(\\.gz)?$", "", basename(reads[1]))
  }

  if (!dir.exists(outdir) && execute) {
    dir.create(outdir, recursive = TRUE)
  }

  summary_file <- file.path(outdir, glue("{sample_id}_rRNA_summary.txt"))

  if (is_paired) {
    out1 <- file.path(outdir, glue("{sample_id}_1.fastq{suffix}"))
    out2 <- file.path(outdir, glue("{sample_id}_2.fastq{suffix}"))
    out_bam <- file.path(outdir, glue("{sample_id}_rRNA_sort.bam"))

    tmp_prefix <- file.path(outdir, glue("{sample_id}_fastq"))

    cmd <- glue(
      "{hisat2_path} ",
      "--summary-file {summary_file} ",
      "--no-spliced-alignment --no-softclip --norc --no-unal ",
      "--threads {ncores} --dta ",
      "--un-conc-gz {tmp_prefix}.gz ",
      "-x {rRNA_index} ",
      "-1 {reads[1]} -2 {reads[2]} | ",
      "{samtools_path} view -@ {ncores} -Shub - | ",
      "{samtools_path} sort -@ {ncores} -o {out_bam} - && ",
      "mv {tmp_prefix}.1.gz {out1} && mv {tmp_prefix}.2.gz {out2}"
    )

    out_files <- c(out1, out2, summary_file, out_bam)
  } else {
    out <- file.path(outdir, glue("{sample_id}.fastq{suffix}"))
    out_bam <- file.path(outdir, glue("{sample_id}_rRNA_sort.bam"))

    cmd <- glue(
      "{hisat2_path} ",
      "--summary-file {summary_file} ",
      "--no-spliced-alignment --no-softclip --norc --no-unal ",
      "--threads {ncores} --dta ",
      "--un-gz {out} ",
      "-x {rRNA_index} ",
      "-U {reads[1]} | ",
      "{samtools_path} view -@ {ncores} -Shub - | ",
      "{samtools_path} sort -@ {ncores} -o {out_bam} -"
    )

    out_files <- c(out, summary_file, out_bam)
  }

  if (execute) {
    existing <- out_files[file.exists(out_files)]
    if (length(existing) > 0 && !overwrite) {
      message("Output files exist and overwrite = FALSE: ", paste(existing, collapse = ", "))
      return(invisible(NULL))
    }

    message("Filtering rRNA reads for sample: ", sample_id)
    system(cmd)
    message("Finished filtering rRNA for sample: ", sample_id)
    return(invisible(out_files))
  } else {
    return(as.character(cmd))
  }
}

#' Align reads using HISAT2 and output BAM file
#'
#' @param sample_id Sample ID used as prefix for output files. If NULL, inferred from reads.
#' @param reads A character vector of 1 (single-end) or 2 (paired-end) fastq files.
#' @param index_base Base path to HISAT2 index (without .ht2 extension).
#' @param outdir Output directory for alignment results.
#' @param ncores Number of threads to use.
#' @param execute Whether to execute or return command as string.
#' @param overwrite Whether to overwrite existing output files.
#' @param hisat2_path Path to hisat2 executable.
#' @param samtools_path Path to samtools executable.
#' @param extra_args Additional arguments passed to hisat2.
#'
#' @return BAM file and summary file path (or command string).
#' @export
hisat2_align <- function(
  sample_id = NULL,
  reads,
  index_base = "/mnt/reference_genomes/gencode/mouse/M36/index/hisat2/GRCm39.MeRIP-seq/genome",
  outdir = "./results/alignment/hisat2",
  ncores = 8,
  execute = FALSE,
  overwrite = FALSE,
  hisat2_path = "/opt/hisat2/2.2.1/bin/hisat2",
  samtools_path = "/opt/samtools/1.21/bin/samtools",
  extra_args = ""
) {
  stopifnot(length(reads) %in% c(1, 2))
  is_paired <- length(reads) == 2
  is_gzipped <- grepl("\\.gz$", reads[1])

  if (is_null(sample_id) || sample_id == "") {
    sample_id <- sub("(_R[12])?(_[12])?\\.(fq|fastq)(\\.gz)?$", "", basename(reads[1]))
  }

  if (!dir.exists(outdir) && execute) {
    dir.create(outdir, recursive = TRUE)
  }

  bam_file <- file.path(outdir, glue("{sample_id}.hisat2.bam"))
  summary_file <- file.path(outdir, glue("{sample_id}.hisat2.summary.txt"))

  if (is_paired) {
    cmd <- glue(
      "{hisat2_path} ",
      "--summary-file {summary_file} ",
      "-p {ncores} --dta {extra_args} ",
      "-x {index_base} ",
      "-1 {reads[1]} -2 {reads[2]} | ",
      "{samtools_path} view -@ {ncores} -hbS - > {bam_file}"
    )
  } else {
    cmd <- glue(
      "{hisat2_path} ",
      "--summary-file {summary_file} ",
      "-p {ncores} --dta {extra_args} ",
      "-x {index_base} ",
      "-U {reads[1]} | ",
      "{samtools_path} view -@ {ncores} -hbS - > {bam_file}"
    )
  }

  out_files <- c(bam_file, summary_file)

  if (execute) {
    existing <- out_files[file.exists(out_files)]
    if (length(existing) > 0 && !overwrite) {
      message("Output files already exist and overwrite = FALSE: ", paste(existing, collapse = ", "))
      return(invisible(NULL))
    }

    message("Running HISAT2 alignment for sample: ", sample_id)
    system(cmd)
    message("Finished alignment for sample: ", sample_id)
    return(invisible(out_files))
  } else {
    return(as.character(cmd))
  }
}
