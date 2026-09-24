# Evaluate code with temporary local usage tracking

Evaluate code with temporary local usage tracking

## Usage

``` r
sn_with_usage_tracking(
  expr,
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

- expr:

  Code to evaluate.

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

The value of `expr`, with its visibility preserved.
