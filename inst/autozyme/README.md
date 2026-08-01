# Vendored AutoZyme patches

Shennong vendors exact-scoped `scdblfinder` and `soupx` patches so its
doublet-detection and ambient-RNA workflows can use AutoZyme even when the
official AutoZyme package does not ship these patches.

- Upstream project: <https://github.com/ElliotXie/autozyme>
- Upstream source revisions: scDblFinder task candidate `c7a87c8` and
  Autozyme package commits `ab5a470`, `97747ea`; SoupX `9952189`
- Vendored files: installed `patches/scdblfinder/patch.R`, its validated scope
  and finalized benchmark summary under `patches/scdblfinder/`, plus
  `patches/soupx/patch.R`, its validated scope in
  `patches/soupx/SCOPE.md`, the successful equivalence-passing benchmark
  summary in `patches/soupx/speedups_finalized.tsv`, and package source
  `src/soupx.cpp`. The recorded speedups use cached baseline measurements from
  AutoZyme's verifier; failed or baseline-less runs are not finalized.
- License: MIT; see `LICENSE.autozyme.md`

Shennong compiles the SoupX native kernel into its own shared library. At
runtime it validates and registers the vendored closures through AutoZyme's
exported `register_patch()` API, activates each patch only around its compatible
workflow call, then restores the caller's previous AutoZyme state.
