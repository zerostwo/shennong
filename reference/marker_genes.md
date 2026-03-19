# Pan-Immune CellTypist Metadata

This dataset contains hierarchical immune cell type metadata from the
Pan-Immune CellTypist atlas. It includes high- and low-hierarchy cell
type annotations, human and mouse marker genes, and curated marker gene
classification.

## Usage

``` r
data(marker_genes)
```

## Format

A data frame with 350 rows and 5 columns:

- high_hierarchy_cell_type:

  Character. Broad cell type category (e.g., "B cells").

- low_hierarchy_cell_type:

  Character. More specific cell type classification (e.g., "Follicular B
  cells").

- human:

  Character. Human marker gene associated with the cell type (e.g.,
  "CD79A").

- mouse:

  Character. Corresponding mouse homolog of the human marker gene (e.g.,
  "Cd79a").

- high_hierarchy_marker_gene:

  Logical. Indicates whether the marker gene is a high-hierarchy marker
  gene (TRUE/FALSE).

## Source

CellTypist Atlas: Pan-Immune dataset, Teichlab [GitHub
Repository](https://github.com/Teichlab/celltypist_wiki)

## Details

**Reference**: C. Dominguez Conde et al., Science, 2022. DOI:
[10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197)

**Original Data Source**: The dataset was obtained from the CellTypist
repository at: [GitHub:
Teichlab/celltypist_wiki](https://github.com/Teichlab/celltypist_wiki/blob/main/atlases/Pan_Immune_CellTypist/v2/tables/Basic_celltype_information.xlsx)

**Description**: A dataset providing hierarchical immune cell types with
their ontology identifiers, curated human and mouse marker genes, and a
classification flag for high-hierarchy marker genes.
