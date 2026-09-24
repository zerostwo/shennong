# Record explicit consent for usage research

This constructor does not contact a server. It creates a bounded consent
receipt that is stored with the local session and required before any
remote DBI connection. The receipt binds its field categories and is
checked for modification before delivery. Shennong never creates a
participant or installation ID.

## Usage

``` r
sn_confirm_usage_consent(
  scope = c("local", "remote_research"),
  policy_id = NULL,
  policy_version = NULL,
  purposes = "performance_research",
  data_categories = c("api_name", "timing", "safe_parameters"),
  expires_at = NULL,
  study_id = NULL
)
```

## Arguments

- scope:

  `"local"` or `"remote_research"`.

- policy_id, policy_version:

  Research policy identifiers. Both are required for remote research.

- purposes:

  Non-empty research purposes.

- data_categories:

  Allowed remote fields. Supported values are `"api_name"`, `"timing"`,
  `"safe_parameters"`, `"package_versions"`, `"error_class"`, and
  `"acceleration"`.

- expires_at:

  Optional future `POSIXct` expiry.

- study_id:

  Optional short non-identifying study label. This is not a participant
  identifier.

## Value

An immutable-style `sn_usage_consent` list.
