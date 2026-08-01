# Vendored SoupX AutoZyme patch

Shennong vendors the exact-scoped `soupx` patch so that
`sn_remove_ambient_contamination(method = "soupx")` can use AutoZyme even
when the official AutoZyme package does not ship this patch.

- Upstream project: <https://github.com/ElliotXie/autozyme>
- Upstream source revision: `9952189`
- Vendored files: installed `patches/soupx/patch.R`, its validated scope in
  `patches/soupx/SCOPE.md`, the successful equivalence-passing benchmark
  summary in `patches/soupx/speedups_finalized.tsv`, and package source
  `src/soupx.cpp`. The recorded speedups use cached baseline measurements from
  AutoZyme's verifier; failed or baseline-less runs are not finalized.
- License: MIT; see `LICENSE.autozyme.md`

Shennong compiles the native kernel into its own shared library. At runtime it
registers the vendored closures through AutoZyme's exported
`register_patch()` API, activates the patch only for the SoupX correction,
then restores the caller's previous AutoZyme state.
