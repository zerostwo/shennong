# Build, validate, and export a Shennong Result Bundle v1

A Result Bundle is the package-owned, JSON-compatible handoff boundary
between a validated Shennong analysis result and an external
orchestrator. It contains no service calls or credentials. Input records
identify exact immutable revisions and SHA-256 digests; artifact records
describe candidate outputs that an authorized external service may
verify and promote.

## Usage

``` r
sn_build_result_bundle(
  result,
  inputs = list(),
  execution = list(),
  artifacts = list(),
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)

sn_validate_result_bundle(bundle, error = TRUE)

sn_export_result_bundle(bundle, path, pretty = TRUE, overwrite = FALSE)
```

## Arguments

- result:

  A canonical Shennong analysis result accepted by
  [`sn_validate_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_validate_result.md).

- inputs:

  A list of immutable input-reference records. Each record requires
  `role`, `revision`, at least one of `resource_id` or `artifact_id`,
  and
  `digest = list(algorithm = "sha256", value = "<64 hex characters>")`.
  Optional fields are `media_type`, `size_bytes`, and `metadata`.

- execution:

  A JSON-compatible list of execution provenance such as job, runtime,
  image, environment-lock, or source revision identifiers. Credentials
  are rejected.

- artifacts:

  A list of candidate output-artifact records. Each record requires
  `role` and a SHA-256 `digest`. Optional fields are `artifact_id`,
  bundle-relative `path`, `media_type`, `size_bytes`, and `metadata`.

- created_at:

  RFC 3339 UTC creation timestamp ending in `Z`. By default the current
  UTC time is used.

- bundle:

  A Result Bundle to validate or export.

- error:

  If `TRUE`, stop for an invalid bundle. If `FALSE`, return a validation
  report.

- path:

  Destination JSON file.

- pretty:

  Write indented JSON.

- overwrite:

  Replace an existing destination file.

## Value

`sn_build_result_bundle()` returns an `sn_result_bundle` list.
`sn_validate_result_bundle()` returns a report with `valid`, `errors`,
and `warnings`. `sn_export_result_bundle()` invisibly returns the
normalized file path, byte size, and SHA-256 digest of the exported
JSON.

## Examples

``` r
result <- list(
  schema_version = "2.0.0", analysis_type = "demo", result_id = "example",
  method = "mean", backend = "base", input = list(), parameters = list(),
  tables = list(primary = data.frame(feature = "gene1", value = 1)),
  embeddings = list(), graphs = list(), models = list(), diagnostics = list(),
  warnings = character(),
  provenance = list(
    package_versions = list(Shennong = "0.2.0"),
    random_seed = 1L,
    timestamp = "2026-01-01 UTC"
  )
)
bundle <- sn_build_result_bundle(result)
sn_validate_result_bundle(bundle, error = FALSE)
path <- tempfile(fileext = ".json")
sn_export_result_bundle(bundle, path)
unlink(path)
```
