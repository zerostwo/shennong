# Summarize local workflow usage and runtime

Summarize local workflow usage and runtime

## Usage

``` r
sn_summarize_usage(
  path = NULL,
  store = NULL,
  source = c("local", "remote"),
  workflow = NULL,
  mode = NULL,
  top_level = TRUE,
  by_parameters = FALSE,
  sort_by = c("calls", "total_seconds", "median_seconds")
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

- top_level:

  If `TRUE`, keep only root calls; if `FALSE`, keep only nested calls;
  `NULL` keeps both.

- by_parameters:

  Group separate sanitized parameter fingerprints. This is useful for
  comparing calls such as different
  [`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md)
  method settings.

- sort_by:

  Rank summaries by call count, cumulative elapsed time, or median
  elapsed time.

## Value

A data frame ordered by call count and total elapsed time.
