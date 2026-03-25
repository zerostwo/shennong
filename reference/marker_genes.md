# Pan-Immune CellTypist Metadata

Hierarchical immune-cell annotations and curated marker genes from the
Pan-Immune CellTypist atlas.

## Usage

``` r
marker_genes
```

## Format

A data frame with 350 rows and 5 columns:

- `high_hierarchy_cell_type`:

  Broad cell type category.

- `low_hierarchy_cell_type`:

  More specific cell type label.

- `human`:

  Human marker gene.

- `mouse`:

  Mouse homolog marker gene.

- `high_hierarchy_marker_gene`:

  Whether the marker is considered high-hierarchy in the source atlas.

## Source

CellTypist Pan-Immune atlas table:
<https://github.com/Teichlab/celltypist_wiki/blob/main/atlases/Pan_Immune_CellTypist/v2/tables/Basic_celltype_information.xlsx>
