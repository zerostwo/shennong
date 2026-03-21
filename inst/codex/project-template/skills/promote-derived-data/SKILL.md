# promote-derived-data

## Purpose

Standardize when run outputs may be promoted into `data/processed/`.

## When To Use

- a run produced a reusable derived asset
- a processed object should become a shared project input

## Required Inputs

- source run ID
- candidate derived asset
- evidence that the asset is reproducible

## Required Outputs

- promoted file in `data/processed/`
- provenance note linking back to the source run

## Rules

- promotion must only occur for reusable derived assets
- promoted assets must have stable names
- provenance must remain traceable to a run

## Procedure

1. Confirm the source run ID.
2. Confirm that the asset is reproducible and reusable.
3. Copy or materialize the asset into `data/processed/` with an explicit name.
4. Record the promotion in `memory/Status.md` or `memory/Decisions.md`.

## Common Mistakes

- promoting transient scratch files
- promoting outputs without run provenance
- using ambiguous filenames

## Examples

- promoting a curated Seurat object from `runs/<run_id>/outputs/` to
  `data/processed/pbmc_integrated_seurat.rds`
