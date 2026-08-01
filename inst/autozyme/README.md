# Vendored SoupX AutoZyme patch

Shennong vendors the exact-scoped `soupx` patch so that
`sn_remove_ambient_contamination(method = "soupx")` can use AutoZyme even
when the official AutoZyme package does not ship this patch.

- Upstream project: <https://github.com/ElliotXie/autozyme>
- Upstream source revision: `9952189`
- Vendored files: installed `patches/soupx/patch.R` plus package source
  `src/soupx.cpp`
- License: MIT; see `LICENSE.autozyme.md`

Shennong compiles the native kernel into its own shared library. At runtime it
registers the vendored closures through AutoZyme's exported
`register_patch()` API, activates the patch only for the SoupX correction,
then restores the caller's previous AutoZyme state.
