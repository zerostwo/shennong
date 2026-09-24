# List local workflow usage records

List local workflow usage records

## Usage

``` r
sn_list_usage_runs(
  path = NULL,
  store = NULL,
  source = c("local", "remote"),
  workflow = NULL,
  mode = NULL,
  status = NULL,
  top_level = NULL,
  limit = 1000L
)
```

## Arguments

- path:

  Existing usage SQLite file. It may be omitted while tracking is
  enabled in the current process.

- store:

  Optional usage store. Use `source = "remote"` with a DBI store and
  reader-authorized connection factory.

- source:

  Read the local outbox or the managed remote table.

- workflow:

  Optional exact workflow name(s).

- mode:

  Optional execution mode(s).

- status:

  Optional run status values.

- top_level:

  If `TRUE`, keep only root calls; if `FALSE`, keep only nested calls;
  `NULL` keeps both.

- limit:

  Maximum number of newest rows returned.

## Value

A data frame ordered from newest to oldest. `invocation_number` is the
use number for the workflow; `parameter_set_invocation` is the use
number for the same sanitized parameter fingerprint.
