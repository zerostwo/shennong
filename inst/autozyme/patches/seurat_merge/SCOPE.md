# seurat_merge release scope

## Verdict

This patch accelerates exactly one SeuratObject namespace method:
`merge.Assay5`. It does **not** register or rebind `merge.Seurat`,
`merge.StdAssay`, or `merge.Assay`. Top-level `merge.Seurat` can benefit only
indirectly when its normal S3 call dispatches to `merge.Assay5`; the package
does not claim a replacement for the whole Seurat merge workflow.

The task root had no `CAVEATS.md` when packaged. The absence was confirmed by
the parent agent; the release boundary below is therefore derived from
`memory/discoveries.md`, the accepted pipeline, and the strengthened task
evaluator.

## Provenance

- Task: `test_seurat_merge`
- Final task best: `eeba1054ef0aa9def7010dbac6a1a559e6306039`
- Audit-remediated parent carrying the packaged `merge.Assay5` implementation:
  `b1cf81af6250c160df88759b9cd2f319d6735915`
- Upstream source: `satijalab/seurat-object`
- Upstream commit: `b9a11fa46acc2142e142f08743b48c73a45b2742`
- Upstream package metadata: `SeuratObject 5.4.0.9001`

The accepted task also patched `merge.StdAssay` and `merge.Seurat`. Those
bindings are intentionally not lifted. The round-56 headline timing measured
that broader task pipeline, so it is evidence of convergence and exactness,
not a package-level speed claim for this narrower one-target release.

## Attested call

The public call represented by the direct Assay5 fixture is:

```r
merge(
  x = x_assay5,
  y = list(y_assay5, ...),
  labels = labels,
  add.cell.ids = NULL,
  collapse = FALSE,
  ...
)
```

The replacement mirrors the upstream formals and adds the package escape
hatch `zyme = TRUE` after `...`. `zyme = FALSE` calls the captured upstream
`merge.Assay5`. `autozyme::with_disabled()` is handled by the namespace
dispatcher and also calls the captured upstream method.

## Supported fast path

All of the following must hold:

- `x` and every member of `y` have the exact S4 class `Assay5` from the
  `SeuratObject` package. Inheritance alone is insufficient.
- `collapse` is exactly `FALSE`.
- `labels` is `NULL` or a character vector with one value per assay.
- `add.cell.ids` is `NULL` or a character vector with one value per assay.
- Generated layer names are unique.
- Layer backends carry only their canonical class/slot attributes. Arbitrary
  user-added backend attributes use upstream fallback because its canonical
  `LayerData<-` path may intentionally discard them.
- Feature-variable metadata names beginning with `vf_` have the canonical
  four-or-more-token shape used by upstream.
- If the private `.cell.names` dot is supplied, it is a list with one complete,
  non-missing character vector per assay. No other dot is evaluated.

The fast path creates the result's layer list and `cells`/`features` LogMaps
directly. Large maps are assembled in one logical matrix per map instead of
being grown repeatedly through `LayerData<-`. Layer objects themselves are
placed into the result without coercion, preserving raw sparse/backend class
and attributes. The canonical Assay5 metadata setter is retained, the default
layer is restored, and `methods::validObject()` is required before return.

## Fallback surface

The captured upstream `merge.Assay5` is used for:

- `zyme = FALSE`;
- `autozyme::with_disabled()` / global disable switches;
- v3 `Assay` inputs;
- `SCTAssay` inputs;
- `Assay5T`, custom `Assay5` subclasses, mixed classes, or malformed inputs;
- `collapse` values outside the exact supported call;
- malformed label, cell-id, private cell-map, layer-name, or variable-feature
  metadata shapes;
- layer backends with noncanonical user-added attributes;
- any unexpected error raised while assembling or validating the fast result.

`merge.StdAssay` and `merge.Assay` remain byte-for-byte upstream namespace
bindings. Their v3/SCT behavior is outside this patch, not emulated by it.

## Correctness contract

The package contract requires exact agreement with the upstream-disabled call
for:

- result class, dimensions, feature/cell order, layer names, default layer and
  default slot;
- raw layer class and attributes, including a real `dgTMatrix` backend;
- `cells` and `features` LogMaps;
- feature metadata and key;
- warnings, messages, and errors on fallback branches;
- caller-owned input serialization before and after the call;
- unrelated lazy promises in `...` remaining unforced.

The accepted task passed all 35 exact metrics on 15 fixtures, including v3,
unsupported subclass, raw backend, input immutability, and direct Assay5 lazy
dots. The package adds focused contract tests for these release boundaries.

## Performance and portability

No compiled kernel, thread pool, fork, external executable, network access, or
platform-specific path is introduced. Linux package-only attest passed three
repetitions per tier: 3.2x/2.5x/1.3x median speedups and
13.1%/19.2%/21.0% lower peak RSS for small/medium/large. Cross-platform speed
and correctness still require the normal CI matrix.
