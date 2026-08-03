# AutoZyme fork integration

Shennong pins the `zerostwo/autozyme` fork at the revision declared in
`DESCRIPTION` and uses its registered patches for `scDblFinder`, UCell, LISI,
and SoupX. These are direct AutoZyme providers; Shennong does not copy their
patch sources into the package.

The ambient-correction path also requests the `decontx_standalone` patch
around the direct `decontX::decontX()` call. That patch is used automatically
when the installed fork build registers it. `decontX::decontPro()` is the
explicit CITE-seq/protein backend and currently has no Shennong-local native
replacement.

Only the exact-scoped `seurat_merge` and `seurat_joinlayers` patches remain
vendored under `patches/`; their source fingerprints and registration logic
are package-owned because the pinned fork does not provide those two entries.
There is no Shennong SoupX or scDblFinder native kernel to maintain.

- AutoZyme fork: <https://github.com/zerostwo/autozyme>
- Pinned revision: see `DESCRIPTION` and `R/acceleration.R`
- Runtime controls: `sn_check_autozyme()`, `sn_enable_autozyme()`,
  `sn_disable_autozyme()`, and `sn_with_autozyme()`
- License: MIT, following the upstream AutoZyme and patch licenses
