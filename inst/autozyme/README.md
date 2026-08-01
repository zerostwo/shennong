# Vendored AutoZyme patches

Shennong vendors exact-scoped `scdblfinder`, `seurat_merge`,
`seurat_joinlayers`, and `soupx` patches so its single-cell workflows can use
them even when the official AutoZyme package does not ship these patches.

- Upstream project: <https://github.com/ElliotXie/autozyme>
- Upstream source revisions: scDblFinder task best
  `34b1fc7a9132b6e646512688018d4be353bba6a5`; Seurat merge
  `eeba1054ef0aa9def7010dbac6a1a559e6306039`; Seurat JoinLayers `55da460`;
  SoupX `9952189`
- Vendored files: installed `patches/scdblfinder/patch.R` and
  its validated scope and finalized benchmark summary under
  `patches/scdblfinder/`; `patches/soupx/patch.R`, plus SoupX's validated scope in
  `patches/soupx/SCOPE.md`, the successful equivalence-passing benchmark
  summary in `patches/soupx/speedups_finalized.tsv`; and the corresponding
  source, scope, manifest, changelog, and benchmark evidence under
  `patches/seurat_merge/` and `patches/seurat_joinlayers/`. The SoupX package
  source is `src/soupx.cpp`. The recorded speedups use cached baseline measurements from
  AutoZyme's verifier; failed or baseline-less runs are not finalized.
- License: MIT; see `LICENSE.autozyme.md`

Shennong compiles the SoupX native kernel into its own shared library. At
runtime it validates each vendored `patch.R` fingerprint and registers the
closures through AutoZyme's exported `register_patch()` API, activates each
patch only around its compatible workflow call, then restores the caller's
previous AutoZyme state.
