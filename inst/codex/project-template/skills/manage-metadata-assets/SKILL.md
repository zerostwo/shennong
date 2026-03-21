# manage-metadata-assets

## Purpose

Standardize how metadata files are stored, curated, and referenced.

## When To Use

- ingesting new sample sheets
- curating annotations
- adding dictionaries or codebooks

## Required Inputs

- source metadata files
- intended metadata role
- provenance details

## Required Outputs

- metadata stored under `data/metadata/`
- clear separation between raw and curated metadata when relevant

## Rules

- metadata must live under `data/metadata/`
- raw metadata must remain immutable
- curated metadata must be traceable to raw inputs or scripts

## Procedure

1. Place source metadata in `data/metadata/raw/` when the raw/curated split is used.
2. Create curated derivatives in `data/metadata/curated/` with explicit names.
3. Store dictionaries or codebooks in `data/metadata/dictionaries/`.
4. Record important metadata assumptions in `memory/Decisions.md`.

## Common Mistakes

- scattering metadata across scripts or notebook folders
- overwriting source sample sheets
- omitting codebooks for encoded labels

## Examples

- `data/metadata/raw/sample_sheet_2026-03-21.csv`
- `data/metadata/curated/sample_sheet_curated_v1.tsv`
