# Bundled Shennong Gene Annotations

A package-owned snapshot of gene-level annotations parsed from the
GENCODE human GRCh38 v48 and mouse GRCm39 vM37 primary-assembly GTF
files.

## Usage

``` r
shennong_gene_annotations
```

## Format

A data frame with one row per gene and the following columns:

- `species`:

  Species label used by Shennong (`"human"` or `"mouse"`).

- `release`:

  GENCODE release identifier used to build the row.

- `seqname`:

  Sequence/chromosome name from the GTF.

- `source`:

  Second GTF column describing the annotation source.

- `start`:

  1-based gene start coordinate.

- `end`:

  1-based gene end coordinate.

- `strand`:

  Gene strand.

- `gene_id`:

  Original GENCODE gene identifier, including version suffix when
  present.

- `gene_id_base`:

  GENCODE gene identifier with the version suffix removed.

- `gene_name`:

  Gene symbol/name from the GTF.

- `gene_type`:

  GENCODE `gene_type` / `gene_biotype` annotation.

- `gene_status`:

  Optional GTF `gene_status` attribute when available.

- `gene_source`:

  Optional GTF `gene_source` attribute when available.

- `level`:

  Optional GENCODE annotation level.

- `gene_class`:

  Shennong-level coarse classification used by filtering helpers.
  Currently `"coding"` marks protein-coding, IG, and TCR genes; all
  other annotated types are treated as `"noncoding"`.

## Source

Built during development from:

- `/mnt/reference_genomes/gencode/human/48/GRCh38.v48.primary_assembly.annotation.gtf`

- `/mnt/reference_genomes/gencode/mouse/M37/GRCm39.vM37.primary_assembly.annotation.gtf`

via `data-raw/build_shennong_gene_annotations.R`.
