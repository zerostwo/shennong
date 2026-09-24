# Create a local or managed-remote usage store

A remote DBI store always retains a local SQLite outbox at `path`.
Scientific calls write only to the outbox;
[`sn_flush_usage_tracking()`](https://zerostwo.github.io/shennong/dev/reference/sn_flush_usage_tracking.md)
performs the explicit, consent-gated remote delivery. The `connect`
closure and its credentials are never serialized into the usage
database.

## Usage

``` r
sn_create_usage_store(
  backend = c("sqlite", "dbi"),
  path,
  connect = NULL,
  table_prefix = "shennong_usage",
  allow_remote_ddl = FALSE
)
```

## Arguments

- backend:

  `"sqlite"` for local-only storage or `"dbi"` for a managed remote DBI
  destination backed by a local SQLite outbox.

- path:

  Local SQLite database/outbox path.

- connect:

  For `backend = "dbi"`, a zero-argument function returning a new valid
  `DBIConnection`. Do not return a connection opened before a fork.

- table_prefix:

  Remote table prefix.

- allow_remote_ddl:

  Permit an explicitly consented flush to create the two remote tables.
  Keep `FALSE` when migrations are administrator-owned.

## Value

An `sn_usage_store` object. It contains a connection factory in memory
and must not be serialized or committed.
