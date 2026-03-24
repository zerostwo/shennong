# Delete a signature from the editable source registry

Delete a signature from the editable source registry

## Usage

``` r
sn_delete_signature(species = "human", path, registry_path = NULL)
```

## Arguments

- species:

  One of `"human"` or `"mouse"`.

- path:

  Slash-delimited signature path relative to the species root.

- registry_path:

  Optional path to the editable registry JSON file.

## Value

Invisibly returns the normalized registry path.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_delete_signature(
  species = "human",
  path = "Programs/custom/MySignature"
)
} # }
```
