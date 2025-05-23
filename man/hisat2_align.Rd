% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/upstream.R
\name{hisat2_align}
\alias{hisat2_align}
\title{Align reads using HISAT2 and output BAM file}
\usage{
hisat2_align(
  sample_id = NULL,
  reads,
  index_base =
    "/mnt/reference_genomes/gencode/mouse/M36/index/hisat2/GRCm39.MeRIP-seq/genome",
  outdir = "./results/alignment/hisat2",
  ncores = 8,
  execute = FALSE,
  overwrite = FALSE,
  hisat2_path = "/opt/hisat2/2.2.1/bin/hisat2",
  samtools_path = "/opt/samtools/1.21/bin/samtools",
  extra_args = ""
)
}
\arguments{
\item{sample_id}{Sample ID used as prefix for output files. If NULL, inferred from reads.}

\item{reads}{A character vector of 1 (single-end) or 2 (paired-end) fastq files.}

\item{index_base}{Base path to HISAT2 index (without .ht2 extension).}

\item{outdir}{Output directory for alignment results.}

\item{ncores}{Number of threads to use.}

\item{execute}{Whether to execute or return command as string.}

\item{overwrite}{Whether to overwrite existing output files.}

\item{hisat2_path}{Path to hisat2 executable.}

\item{samtools_path}{Path to samtools executable.}

\item{extra_args}{Additional arguments passed to hisat2.}
}
\value{
BAM file and summary file path (or command string).
}
\description{
Align reads using HISAT2 and output BAM file
}
