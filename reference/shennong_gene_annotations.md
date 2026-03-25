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

  Species label used by Shennong.

- `release`:

  GENCODE release identifier.

- `seqname`:

  Sequence or chromosome name from the GTF.

- `source`:

  Annotation source from the GTF.

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

  Version-stripped GENCODE gene identifier.

- `gene_name`:

  Gene symbol from the GTF.

- `gene_type`:

  GENCODE gene biotype.

- `gene_status`:

  Optional GTF `gene_status` attribute.

- `gene_source`:

  Optional GTF `gene_source` attribute.

- `level`:

  Optional GENCODE annotation level.

- `gene_class`:

  Shennong coarse classification used by gene filtering helpers.

## Source

Built from local GENCODE sources via
`data-raw/build_shennong_gene_annotations.R`.
