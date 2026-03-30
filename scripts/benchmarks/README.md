# Shennong Benchmark Assets

Purpose:
- Keep benchmark work versioned in the repository rather than as ad hoc notes.
- Start with lightweight workflow-surface benchmarking before adding heavier
  reproducibility, annotation, and external-tool benchmarks.

Files:
- `scripts/pbmc_shennong_workflow.R`: representative Shennong-first workflow
- `scripts/pbmc_seurat_baseline.R`: representative Seurat-first baseline
- `compare_workflow_surface.R`: static script audit that generates initial
  workflow-surface metrics
- `results/`: generated benchmark outputs

Usage:

```r
Rscript scripts/benchmarks/compare_workflow_surface.R
```

Current scope:
- Phase 0 static workflow-surface audit only

Planned additions:
- reproducibility replay benchmark
- AI-assisted annotation benchmark
- external-tool provenance benchmark
