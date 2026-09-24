# Enable local workflow usage and timing records

Instruments Shennong's high-level computational workflows for the
current R process and writes one local SQLite row per call. Tracking is
always opt-in: loading or attaching Shennong never creates a database,
and no record is sent over the network. The parameter recorder stores
conservative scalar settings but redacts objects, matrices, identifiers,
paths, credentials, prompts, and free text.

## Usage

``` r
sn_enable_usage_tracking(
  path = NULL,
  store = NULL,
  consent = NULL,
  mode = c("development", "production", "test", "benchmark"),
  display = FALSE,
  nested = TRUE,
  functions = NULL,
  remote_flush = c("manual", "on_disable"),
  strict = FALSE
)
```

## Arguments

- path:

  Local SQLite file path. It is retained for compatibility and cannot be
  combined with `store`.

- store:

  Optional object from
  [`sn_create_usage_store()`](https://zerostwo.github.io/shennong/dev/reference/sn_create_usage_store.md).
  A remote DBI store still writes first to its local SQLite outbox.

- consent:

  Optional object from
  [`sn_confirm_usage_consent()`](https://zerostwo.github.io/shennong/dev/reference/sn_confirm_usage_consent.md).
  It is mandatory for a remote store.

- mode:

  Explicit execution context: `"development"`, `"production"`, `"test"`,
  or `"benchmark"`. Shennong never infers that a run is production.

- display:

  Show a compact elapsed-time message when each tracked call finishes.

- nested:

  Record nested high-level Shennong workflow calls. Summaries exclude
  nested calls by default to avoid double-counting elapsed time.

- functions:

  Optional exact character vector of exported functions to instrument.
  `NULL` records every public Shennong and rio adapter except the
  usage-tracking control/query functions themselves.

- remote_flush:

  For a DBI store, flush manually or when tracking is disabled. Remote
  failures on disable warn and retain the local outbox.

- strict:

  If `TRUE`, a database write failure stops the analysis. The default is
  fail-open: warn once and preserve the scientific call.

## Value

Invisibly, the result of
[`sn_check_usage_tracking()`](https://zerostwo.github.io/shennong/dev/reference/sn_check_usage_tracking.md).

## Examples

``` r
if (FALSE) { # \dontrun{
database <- file.path(tempdir(), "shennong-usage.sqlite")
sn_enable_usage_tracking(database, mode = "development")
object <- sn_run_cluster(object, integration_method = "unintegrated")
sn_disable_usage_tracking()
sn_summarize_usage(database)
} # }
```
