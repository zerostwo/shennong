# Structural Roadmap

Last updated: 2026-07-14

## Next safe refactors

- Split large durable domains only when a focused test boundary exists:
  `analysis_metrics.R`, `package_tools.R`, and `interpretation.R` are the next
  candidates.
- Add focused source-level tests where CodeGraph cannot infer dynamic
  `object@misc` dispatch or test reachability.
- Keep benchmark inputs and run logs outside git; retain only scripts, compact
  summaries, and reproducible metadata in `benchmarks/`.

## Deferred local-data decision

Ignored `dev/outputs/` and benchmark input/run caches occupy several gigabytes
and contain untracked work. They are intentionally not deleted by automated
repository cleanup. Review and migrate valuable scripts or source data before
removing those caches.

## Guardrails

- Preserve exported API behavior during file splits.
- Do not delete generated `man/` files independently of roxygen sources.
- Do not remove `data/*.rda`, `inst/pixi/`, `inst/codex/`, vignette sources, or
  CI workflows merely because CodeGraph does not index them.
- Re-sync CodeGraph after structural changes and verify both new and removed
  symbols explicitly.
