# Human-Mouse Homologous Genes Table

This dataset contains homologous gene mappings between human and mouse,
including gene symbols, HGNC IDs, and MGI IDs.

## Usage

``` r
data(hom_genes)
```

## Format

A data frame with 5 columns:

- DB.Class.Key:

  Integer. A unique key for each gene homology group in the MGI
  database.

- human:

  Character. Approved human gene symbol from HGNC.

- HGNC.ID:

  Character. The corresponding HGNC ID for the human gene.

- mouse:

  Character. Approved mouse gene symbol from MGI.

- Mouse.MGI.ID:

  Character. The corresponding MGI ID for the mouse gene.

## Source

MGI Homology Report:
[HOM_MouseHumanSequence.rpt](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt)

## Details

The data was sourced from the Mouse Genome Informatics (MGI) Homology
report: **HOM_MouseHumanSequence.rpt**.

**Reference Databases**:

- **MGI (Mouse Genome Informatics)**: <http://www.informatics.jax.org/>

- **Original Data Source**:
  [HOM_MouseHumanSequence.rpt](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt)

- **HGNC (HUGO Gene Nomenclature Committee)**:
  <https://www.genenames.org/>

- **Ensembl Biomart**: <https://www.ensembl.org/biomart>
