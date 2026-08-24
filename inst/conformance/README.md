# Shennong backend conformance contracts

This directory contains machine-readable contracts for differential checks
between a Shennong entry point and the upstream implementation it wraps.
Contracts are installed with the package so that released validation claims can
be traced to an exact package source.

Contract status is intentionally conservative:

- `pilot` means that the declared unit-scale parameter envelope has an
  executable upstream comparison. It is not a release-conformance claim.
- `admitted` is reserved for a new registered method that has passed the full
  admission requirements in `docs/codex/BackendConformance.md`, including the
  required real/OOD evidence and zero unexpected skips in its CI profile.
- `stale` means that an upstream version, fixture, runner, parameter envelope,
  or comparator changed after the last accepted evidence.

The generated `public-api-parameters.json` inventory and files below
`matrices/` are static coverage artifacts. They classify current formals,
selectors, pairwise axes, and required high-risk cases so API drift cannot pass
silently. Their presence, freshness tests, and `runtime_classified` fields do
not mean that every listed case was executed or passed. Executed evidence comes
only from a contract runner and must report passed, skipped, unsupported, or
failed cells separately.

The repository-only immutable baseline at
`tests/conformance/legacy-methods.txt` prevents newly registered methods from
silently inheriting the historical backlog. The separate
`legacy-pending-methods.txt` is allowed only to shrink as existing methods are
admitted. Do not extend either inventory for a new method; add an `admitted`
contract and its evidence instead.
