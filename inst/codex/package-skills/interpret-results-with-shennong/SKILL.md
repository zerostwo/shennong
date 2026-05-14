---
name: interpret-results-with-shennong
description: Use when interpreting stored Shennong analysis outputs, building prompt bundles from Seurat object results, or generating annotation notes, results prose, figure legends, and presentation summaries from packaged evidence helpers.
---

# interpret-results-with-shennong

## Purpose

Teach the agent how to interpret stored Shennong outputs through the package
interpretation layer.

## When To Use

- summarizing markers, enrichment, or results
- building prompt bundles from stored evidence
- writing annotation notes or result summaries

## Required Inputs

- Seurat object with stored results
- target interpretation task
- optional background context

## Required Outputs

- interpreted result
- prompt bundle or human-readable summary

## Rules

- use stored results rather than recomputing upstream analyses unnecessarily
- separate evidence from inference
- keep interpretation provider-agnostic at the package boundary
- use current public argument names such as `cluster_by`, not retired aliases
  such as `cluster_col`
- if the right interpretation helper is unclear, consult the shared package
  API map before composing raw prompt logic
- if the interpretation belongs to a governed project, update `memory/Status.md`
  or `memory/Decisions.md` when the result changes durable project knowledge

## References

- `../_shared/references/package_api_map.md`
- `../_shared/references/workflow_recipes.md`

## Procedure

1. Discover stored assets with `sn_list_results()`.
2. Retrieve them with the corresponding `sn_get_*_result()` helper.
3. Build prompt bundles with `sn_build_prompt()` when a prompt artifact is
   enough, or run the interpretation layer with `sn_interpret_annotation()`,
   `sn_interpret_de()`, `sn_interpret_enrichment()`, `sn_write_results()`,
   `sn_write_figure_legend()`, or `sn_write_presentation_summary()`.
4. Use `sn_prepare_annotation_evidence()`, `sn_prepare_de_evidence()`,
   `sn_prepare_enrichment_evidence()`, or `sn_prepare_results_evidence()` to
   keep evidence assembly explicit.
5. Use `sn_make_ellmer_provider()` and `sn_test_llm_provider()` when a live
   provider call is actually needed, then `sn_run_llm()` for direct message
   execution through the resolved provider.

## Common Mistakes

- attempting interpretation before storing results
- hiding the source evidence
- treating speculative interpretation as validated analysis

## Examples

- `sn_prepare_annotation_evidence()`
- `sn_prepare_results_evidence()`
- `sn_build_prompt()`
- `sn_interpret_annotation()`
- `sn_interpret_de()`
- `sn_interpret_enrichment()`
- `sn_make_ellmer_provider()`
- `sn_test_llm_provider()`
- `sn_run_llm()`
- `sn_write_results()`
- `sn_write_figure_legend()`
- `sn_write_presentation_summary()`
