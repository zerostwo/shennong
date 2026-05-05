# Upload reusable data files to Zenodo

`sn_upload_zenodo()` creates or updates a Zenodo draft through zen4R,
uploads user-selected files, and writes a Shennong manifest alongside
those files so downstream users can verify the exact dataset version and
file checksums they reused.

## Usage

``` r
sn_upload_zenodo(
  files,
  title,
  description = NULL,
  creators = NULL,
  version = NULL,
  record_id = NULL,
  new_version = FALSE,
  publish = FALSE,
  sandbox = FALSE,
  token = NULL,
  license = "cc-by-4.0",
  resource_type = "dataset",
  access = c("public", "restricted"),
  keywords = NULL,
  manifest = TRUE,
  manifest_dir = tempdir(),
  reserve_doi = TRUE,
  dry_run = FALSE,
  logger = NULL
)
```

## Arguments

- files:

  Character vector of local files to upload.

- title:

  Record title.

- description:

  Record description. When `NULL`, Shennong writes a short reusable-data
  description.

- creators:

  Creator names. A character vector such as `c("Duan, Songqi")` is the
  simplest form. A data frame or list with `name`, `firstname`,
  `lastname`, `orcid`, and `affiliation`/`affiliations` fields is also
  accepted.

- version:

  Dataset version string recorded in Zenodo metadata and the Shennong
  manifest.

- record_id:

  Optional Zenodo record ID. When supplied, `new_version = FALSE`
  updates that draft; `new_version = TRUE` creates a new version from
  the published record.

- new_version:

  If `TRUE`, create a new Zenodo version for `record_id`.

- publish:

  If `TRUE`, publish the record after upload. Defaults to `FALSE` so
  users can inspect the draft first.

- sandbox:

  If `TRUE`, use the Zenodo sandbox service.

- token:

  Zenodo access token. When `NULL`, Shennong checks
  `ZENODO_SANDBOX_TOKEN` for sandbox uploads and then `ZENODO_TOKEN` /
  `ZENODO_PAT`.

- license:

  License identifier stored in Zenodo rights metadata.

- resource_type:

  Zenodo resource type ID. Defaults to `"dataset"`.

- access:

  File and record access policy. Defaults to `"public"`.

- keywords:

  Optional subject keywords.

- manifest:

  If `TRUE`, upload a `shennong_zenodo_manifest.json` file with
  checksums and provenance.

- manifest_dir:

  Directory used to write the temporary manifest file.

- reserve_doi:

  Whether to reserve a DOI for new drafts.

- dry_run:

  If `TRUE`, build and return the planned record and manifest without
  contacting Zenodo.

- logger:

  Logger level passed to `zen4R::ZenodoManager$new()`.

## Value

A list containing the Zenodo identifiers, uploaded-file manifest,
manifest path, publication status, and the `ZenodoRecord` object.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- sn_upload_zenodo(
  files = c("data/reference_atlas.qs", "data/reference_markers.csv"),
  title = "Reusable PBMC reference atlas",
  creators = "Duan, Songqi",
  version = "2026.05.04",
  sandbox = TRUE
)
result$doi
result$manifest
} # }
```
