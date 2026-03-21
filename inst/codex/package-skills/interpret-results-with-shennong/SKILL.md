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
- if the interpretation belongs to a governed project, update `memory/Status.md`
  or `memory/Decisions.md` when the result changes durable project knowledge

## Procedure

1. Discover stored assets with `sn_list_results()`.
2. Retrieve them with the corresponding `sn_get_*_result()` helper.
3. Build or run the interpretation layer with `sn_build_prompt()` or
   `sn_interpret_*()`.

## Common Mistakes

- attempting interpretation before storing results
- hiding the source evidence
- treating speculative interpretation as validated analysis

## Examples

- `sn_prepare_annotation_evidence()`
- `sn_build_prompt()`
- `sn_interpret_annotation()`
