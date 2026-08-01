# seurat_joinlayers release scope

## Verdict

This patch accelerates exactly `SeuratObject::JoinLayers.Assay5`. It does not
replace the generic or any Seurat-level method. The task root had no
`CAVEATS.md`; this boundary comes from the accepted implementation and the
strengthened 15-tier evaluator.

## Provenance

- Task: `test_seurat_joinlayers`
- Accepted task commit: `55da460`
- Upstream: `satijalab/seurat-object` at
  `b9a11fa46acc2142e142f08743b48c73a45b2742`
- Upstream metadata: `SeuratObject 5.4.0.9001`

## Supported fast path

The object must have the exact plain `Assay5` class from SeuratObject. The
public call must explicitly join `layers = "counts"` into `new = "counts"`
with no additional dots, and at least two selected layers must exist. Selected
layers must be uniformly exact `dgCMatrix` objects or public BPCells
`IterableMatrix` objects with valid feature/cell LogMap membership.

The dgC route performs a direct aligned cbind when membership is contiguous,
otherwise it delegates stitching to captured SeuratObject machinery. The
BPCells route uses only public sparse composition and `matrix_type()` checks;
it pads missing rows, restores selected cell order, and mirrors upstream
post-stitch preparation. Both routes update the layer list and LogMaps once,
restore variable features through the canonical setter, and require a valid
Assay5 before return.

## Fallback surface

Captured upstream behavior is retained for custom Assay5 subclasses, v3/SCT
objects, non-count layers, renamed output layers, nonempty dots, mixed or
unknown backends, overlapping cells, malformed maps, heterogeneous BPCells
matrix types, unsupported attributes, missing dependencies, and every failed
proof or assembly check. `autozyme::with_disabled()` and global disable flags
provide the package-level opt-out without changing the upstream method's
formals or forcing lazy dots.

## Correctness and performance evidence

All 15 task tiers passed exact comparisons for raw/accessor class, attributes,
dimnames and payload; selected and unselected layers; feature/cell order and
LogMaps; metadata; warnings; input immutability; fast-route expectation; and
object validity. Auto-structure checks were 12/12 on every tier.

On the accepted run, peak RSS changed from 672.7 to 563.9 MB on small,
745.2 to 605.8 MB on medium, and 1077.9 to 850.8 MB on large inputs. Runtime
changed from 2.224 to 2.030 seconds, 3.177 to 2.299 seconds, and 7.861 to 3.282
seconds respectively. These are task-host measurements, not cross-platform
guarantees. One small partial-feature edge case improved memory but can be
slower; it remains covered to preserve the broad correctness boundary.
