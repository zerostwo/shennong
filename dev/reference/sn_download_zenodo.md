# Download files from a Zenodo record

`sn_download_zenodo()` downloads one or more files from a Zenodo record
into a local cache directory. Public Zenodo records do not require a
token. Supply `token` only for restricted files or private records that
your account can access.

## Usage

``` r
sn_download_zenodo(
  record_id,
  files = NULL,
  save_dir = "~/.shennong/data",
  token = NULL,
  sandbox = FALSE,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- record_id:

  Zenodo record ID.

- files:

  Character vector of file names to download. When `NULL`, the function
  queries the Zenodo record metadata and downloads all listed files.

- save_dir:

  Local directory where files are cached.

- token:

  Optional Zenodo access token for restricted/private downloads.

- sandbox:

  If `TRUE`, use the Zenodo sandbox service.

- overwrite:

  If `TRUE`, re-download files that already exist locally.

- quiet:

  If `TRUE`, suppress download progress messages.

## Value

A named character vector of local file paths.

## Examples

``` r
if (FALSE) { # \dontrun{
paths <- sn_download_zenodo(
  record_id = "14884845",
  files = "pbmc3k_filtered_feature_bc_matrix.h5"
)
paths
} # }
```
