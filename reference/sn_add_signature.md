# Add a signature to the editable source registry

Add a signature to the editable source registry

## Usage

``` r
sn_add_signature(
  species = "human",
  path,
  genes,
  catalog_path = NULL,
  overwrite = FALSE,
  source = "custom"
)
```

## Arguments

- species:

  One of `"human"` or `"mouse"`.

- path:

  Slash-delimited signature path relative to the species root, for
  example `"Programs/MyProgram/MySignature"`.

- genes:

  Character vector of gene symbols stored at the leaf node.

- catalog_path:

  Optional path to a signature catalog `.rda` snapshot.

- overwrite:

  If `TRUE`, replace an existing signature at the same path.

- source:

  Optional source label recorded on the signature node.

## Value

Invisibly returns the normalized registry path.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_add_signature(
  species = "human",
  path = "Programs/custom/MySignature",
  genes = c("GENE1", "GENE2")
)
} # }
```
