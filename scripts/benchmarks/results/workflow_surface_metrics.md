# Workflow Surface Metrics

Generated: 2026-03-30 21:29:32 UTC

This phase-0 benchmark is a static script audit intended to quantify early
workflow-surface differences between a Shennong-first and a Seurat-first
implementation of a routine clustering plus marker-discovery workflow.

| workflow | file | executable_lines | assignment_steps | unique_function_calls | shennong_calls | seurat_calls | persistence_calls | conversion_calls |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| shennong | /home/duansq/personal/packages/shennong/scripts/benchmarks/scripts/pbmc_shennong_workflow.R | 35 |  6 | 10 | 5 | 0 | 2 | 0 |
| seurat_baseline | /home/duansq/personal/packages/shennong/scripts/benchmarks/scripts/pbmc_seurat_baseline.R | 50 | 10 | 13 | 0 | 8 | 1 | 0 |
