# Update a signature in the editable source registry

Update a signature in the editable source registry

## Usage

``` r
sn_update_signature(
  species = "human",
  path,
  genes = NULL,
  registry_path = NULL,
  rename_to = NULL,
  source = NULL
)
```

## Arguments

- species:

  One of `"human"` or `"mouse"`.

- path:

  Slash-delimited signature path relative to the species root.

- genes:

  Optional replacement gene vector. If `NULL`, keep the existing genes.

- registry_path:

  Optional path to the editable registry JSON file.

- rename_to:

  Optional new terminal node name.

- source:

  Optional source label recorded on the signature node.

## Value

Invisibly returns the normalized registry path.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_update_signature(
  species = "human",
  path = "Programs/custom/MySignature",
  genes = c("GENE1", "GENE2", "GENE3")
)
} # }
```
