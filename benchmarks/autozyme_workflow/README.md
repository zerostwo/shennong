# AutoZyme Workflow Patch Benchmark

This benchmark measures actual Shennong call paths rather than isolated patch
functions. Each operation uses three fresh-process arms:

1. direct upstream public API with AutoZyme disabled;
2. the Shennong wrapper with automatic acceleration disabled;
3. the same Shennong wrapper with its lazy AutoZyme scope enabled.

The runner rotates arm order, records analytical elapsed time and whole-worker
peak RSS, canonicalizes the scientific output, and requires both direct-to-
wrapper and off-to-on equivalence. Patch eligibility, target intersection,
scoped activation, and internal fast-path hits remain separate fields. The
current AutoZyme API has no per-call hit counter, so `fast_path_hit` remains
unknown even when parity plus a repeated speed difference provides strong
empirical evidence.

## Fixture

`scripts/real-data/prepare-autozyme-pbmc30k.R` verifies and combines three
independent public PBMC captures:

- 10x PBMC 10k v3 source object SHA-256
  `c30f8ba7d4af1e136f5100cafd21d16f27090560646948a57e39044ba3b1cadb`;
- TENxPBMCData PBMC 33k source object SHA-256
  `263b81043cfd489861e72d17da9b62f283306ffe271d3391c462e37c5e513966`;
- 10x Fresh PBMC 68k donor A source object SHA-256
  `fd0b4cfb7737d9d5e48d02f8ce5a64979a0e7226d78bf5c9eaaeee532c03551d`.

At seed 717 it samples 10,000 cells per capture and retains 20,453 common gene
symbols. The ignored local `.qs2` artifact used below has SHA-256
`20b884fb6c8f20c932ae0eb3e56484e40089386208653156b447439032945013`.
Raw and derived matrices are intentionally not tracked by Git.

## Current pinned-source evidence

`results/pbmc30k-summary.csv` was generated on 2026-08-21 with R 4.6.1,
Seurat 5.5.1, SeuratObject 5.4.0, and exact AutoZyme revision
`8fc2e9c3a7f70302f97589aaa9b0395dcf86f9bc` on an Intel Core Ultra 9
285HX host. Every listed operation has three repetitions with rotating arm
order and a passing declared output comparator. The NormalizeData row comes
from its focused post-optimization rerun; the other main rows come from the
same pinned-library full matrix.

Seurat 5.5.1 and UCell 2.17.0 are explicit upstream-version-drift evidence,
not strict release admission for the manifest's older tested labels. The
source revision itself is exact. The separate GO-cache row is a one-repetition
focused test because it showed only a 1.04x difference and is not a general
enrichment-kernel acceleration.

## Automatic policy is narrower than this benchmark

The automatic Shennong subset is `cellchat`, `clusterprofiler`, `lisi`,
`nichenetr`, `scdblfinder`, `seurat`, `seurat_merge`, `soupx`, and `ucell`,
each behind its documented operation/input guard. Coralysis, standalone
decontX, broad Seurat targets beyond NormalizeData, JoinLayers, tradeSeq, and
WGCNA remain explicit-only. A row in `patch-intersection-audit.csv`, an eligible
manual patch, or the passing JoinLayers timing does not change that policy.

Likewise, the generated public-API parameter inventory plus clustering and
enrichment method matrices are static completeness/admission plans. They
enumerate selectors, dispatch cells and required cases; they are not
runtime-pass records. Only an executed direct / Shennong-off / Shennong-on
comparison with its comparator supports the evidence summarized here.

## Reproduce

```sh
Rscript scripts/real-data/prepare-autozyme-pbmc30k.R

Rscript scripts/real-data/benchmark-autozyme-workflow-patches.R \
  --input data-local/autozyme-benchmark/pbmc-captures-30k.qs2 \
  --output data-local/autozyme-benchmark/pbmc30k.json \
  --sample-by capture \
  --repetitions 3
```

Use an isolated library containing the exact AutoZyme revision declared in
`DESCRIPTION`. `sn_check_autozyme()` now always requires a pinned or exact
trusted patch source; `strict = FALSE` relaxes only the upstream version label.
