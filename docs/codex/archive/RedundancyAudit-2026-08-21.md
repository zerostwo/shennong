# Repository Redundancy Audit — 2026-08-21

The workspace occupied about 5.8 GB excluding `.git` at audit time. Tracked
package source accounted for about 7.7 MB; almost all disk use came from ignored
real-data fixtures, experimental outputs, benchmark inputs, and generated
build products. Exact SHA comparison found no duplicate tracked source files
and no duplicate R function definitions.

## Removed or safely reproducible

Use `Rscript scripts/clean-generated.R` for a dry run and add `--apply` only
after reviewing its exact whitelist. It selects R check/build products,
`README.html`, OmniPath logs, test plots, and Python `__pycache__` directories.
It never selects `.Rhistory`, `.Rproj.user`, `.Renviron`, `.codegraph`, Git
state, `data-local`, `dev`, or benchmark inputs. The ignored pkgdown `site/`
requires the separate `--include-site` flag because it can be useful for local
preview or publication.

All selected paths are reproducible from tracked source. The script refuses
symlinks, paths outside the repository, and any target containing a Git-tracked
file.

The whitelist was applied after validation on 2026-08-21 and removed 20.71 MiB:
`Shennong.Rcheck/`, the source tarball, `README.html`, two OmniPath log trees,
the generated test PDF, and seven Python bytecode caches. Git also pruned the
stale `shennong-coralysis-check` worktree metadata whose target no longer
existed. Scientific fixtures, research outputs, active worktrees, user history,
were retained. The stale 9.01 MiB pkgdown preview was subsequently removed and
replaced by a complete current core real-data build rather than left deleted.
After the final 0.3.0.9000 test/build/check pass on 2026-08-22, the whitelist
was applied again to the regenerated check directory, both obsolete/current
source tarballs, root OmniPath log, test plot, and seven bytecode caches. That
pass removed 22.92 MiB; a following dry run selected zero paths. The final
`site/dev` preview was intentionally retained.

## Large artifacts retained deliberately

| Area | Approximate size | Decision |
|---|---:|---|
| `dev/outputs/` | 2.50 GiB | Retain. These are staged ILC/PBMC research objects, not proven garbage. Archive binary outputs only after their scripts, summaries, manifests, and source data are governed. |
| `data-local/` | 1.30 GiB | Retain. It contains real pkgdown, scIB, and AutoZyme fixtures used for runtime evidence. |
| `benchmarks/coralysis_capacity/results*/inputs/` | 1.86 GiB | Retain pending archival. The inputs are generated, but their caller-supplied base object is not versioned, so exact regeneration is not yet guaranteed. |
| `.codegraph/` | 6.4 MiB | Retain while CodeGraph is used; it is a regenerable developer index. |
| `.Rproj.user/` | under 1 MiB on disk | Retain. Historical recovery copies include former R scripts and must be reviewed before any manual cleanup. |

The audit found about 68.7 MiB of exact duplicate content inside ignored
`data-local/`: repeated scIB input matrices, AutoZyme direct/off canonical
objects, and pkgdown runtime HTML/coverage copies. These paths should not be
deleted independently. Their runners should first migrate to one canonical or
content-addressed artifact, followed by an executable consumer check.

## Source-level consolidation candidates

Three pairs of internal helpers are structurally identical and are candidates
for small reviewed refactors, not immediate deletion:

- `.sn_default_scvi_run_dir()` and `.sn_default_python_run_dir()`;
- `.sn_validate_seurat_object()` and `.sn_validate_result_object()`;
- `.sn_resolve_spatial_result()` and `.sn_resolve_dynamics_result()`.

Standalone scripts also repeat `.absolute_path`, `.option`, `.flag`, and
argument parsers. Consolidating those helpers would reduce maintenance, but a
shared source must preserve the ability to execute each script directly.

The paired `README.Rmd`/`README.md`, roxygen source/`man`/`NAMESPACE`, and
`docs/codex`/`inst/codex` trees are intentional source/generated or
maintainer/installed layers. They are not redundant and must remain synchronized.

## Documentation overlap

No article is an exact duplicate. The older navigation makes some scopes look
duplicated because `annotation-pathways` also discusses programs, enrichment,
communication, and regulation; `clustering` includes integration diagnostics;
and composition/Milo overlaps abundance prioritization. The remedy is a
research-narrative landing page and grouped pkgdown navigation, not deleting
the executable module articles.

Legacy planning files under `dev/` should be moved to `docs/codex/archive/`
only after their still-relevant decisions are reconciled. In particular,
`dev/Shennong_Benchmark_Framework.md` references an obsolete
`scripts/benchmarks/` location. No `dev/` source or output was removed in this
audit.
