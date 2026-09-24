# Flush the local usage outbox to a managed remote DBI store

Delivery is explicit, consent-gated, and idempotent by `session_id` and
`run_id`. Only sessions recorded under the exact consent receipt are
sent; other receipts remain pending. Optional remote fields follow the
receipt's data categories. Remote error messages and process identifiers
are never uploaded. A failed flush leaves rows in the local outbox for
later retry.

## Usage

``` r
sn_flush_usage_tracking(store, consent, limit = 1000L, strict = TRUE)
```

## Arguments

- store:

  A `backend = "dbi"` store from
  [`sn_create_usage_store()`](https://zerostwo.github.io/shennong/dev/reference/sn_create_usage_store.md).

- consent:

  Active `remote_research` consent from
  [`sn_confirm_usage_consent()`](https://zerostwo.github.io/shennong/dev/reference/sn_confirm_usage_consent.md).
  It must exactly match the recorded receipt.

- limit:

  Maximum completed outbox rows delivered in one call.

- strict:

  Stop on a remote error. `FALSE` warns and retains the outbox.

## Value

A one-row delivery summary.
