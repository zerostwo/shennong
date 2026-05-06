# List datasets available through Shennong

`sn_list_datasets()` returns a sample-level registry for datasets that
can be loaded with
[`sn_load_data`](https://songqi.org/shennong/dev/reference/sn_load_data.md).
By default it reads the current Shennong public Zenodo collection, whose
machine-readable index follows the `UPLOAD_RULES.md` layout stored in
the record.

## Usage

``` r
sn_list_datasets(
  source = c("public", "examples", "all"),
  record_id = NULL,
  save_dir = "~/.shennong/data",
  token = NULL,
  overwrite = FALSE,
  quiet = TRUE
)
```

## Arguments

- source:

  Which registry to list. `"public"` lists the current Shennong public
  Zenodo collection; `"examples"` lists the legacy PBMC examples;
  `"all"` returns both.

- record_id:

  Zenodo record ID for the public Shennong collection. Defaults to
  option `shennong.public_data_record` or `"20044788"`.

- save_dir:

  Local cache directory used to store `shennong_index.json`.

- token:

  Optional Zenodo access token for restricted/private records.

- overwrite:

  Logical; if `TRUE`, re-download the public index.

- quiet:

  Logical; if `TRUE`, suppress Zenodo download progress.

## Value

A data frame with one row per loadable sample-level dataset.

## Examples

``` r
if (FALSE) { # \dontrun{
datasets <- sn_list_datasets()
head(datasets[, c("dataset", "study_id", "organism", "technology")])
} # }
```
