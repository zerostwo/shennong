# Backend conformance test materials

`legacy-methods.txt` is the immutable baseline of methods that were already
marked `implemented: true` when the fail-closed conformance gate was introduced.
`legacy-pending-methods.txt` is the shrinking, currently unadmitted subset. Both
are migration accounting, not evidence that those methods conform. Promotion
removes a key from the pending file after its contract becomes `admitted`; the
baseline remains unchanged so post-gate additions stay detectable.

The static test requires every implemented `task::method` key to be either:

1. in the current legacy-pending inventory; or
2. referenced by an installed contract whose status is `admitted`.

Adding a new method to the legacy file is prohibited. A new method must include
an upstream reference runner, complete parameter mapping, semantic comparator,
required fixtures/evidence, and an admitted contract under
`inst/conformance/contracts/`.

Run the current static and micro-differential gates with:

```sh
Rscript -e 'testthat::test_local(filter = "backend-conformance", stop_on_failure = TRUE)'
```
