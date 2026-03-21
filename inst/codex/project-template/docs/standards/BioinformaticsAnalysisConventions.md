# Bioinformatics Analysis Conventions

These conventions are mandatory for this initialized analysis project.

## 1. Top-Level Structure

The repository must use these top-level directories:

- `config/`
- `data/raw/`
- `data/processed/`
- `data/metadata/`
- `scripts/`
- `notebooks/`
- `runs/`
- `results/figures/`
- `results/tables/`
- `results/reports/`
- `memory/`
- `skills/`

The agent must not create alternative top-level analysis directories unless the
change is recorded in `memory/Decisions.md`.

## 2. Raw vs Processed vs Runs vs Results

- `data/raw/` must contain immutable source inputs only.
- `data/processed/` must contain reusable derived assets that have been promoted
  from a formal run.
- `runs/` must contain run-specific working artifacts, transient outputs,
  parameterized execution outputs, and per-run logs.
- `results/` must contain curated deliverables for review, reporting, or export.

The agent must not treat `runs/` as a substitute for `data/processed/`.
The agent must not promote a run output into `data/processed/` until the
promotion criteria are satisfied.

## 3. Metadata Placement

All metadata must live under `data/metadata/`.

Recommended substructure:

- `data/metadata/raw/`
- `data/metadata/curated/`
- `data/metadata/dictionaries/`

If the substructure is used, raw metadata must remain immutable and curated
metadata must be script-derived or explicitly documented.

## 4. Script Organization

Scripts must be organized by analysis stage, not by programming language.

Allowed examples:

- `scripts/01_ingest/`
- `scripts/02_qc/`
- `scripts/03_cluster/`
- `scripts/04_de/`
- `scripts/05_enrichment/`
- `scripts/06_reporting/`

The agent must not organize formal scripts primarily into folders such as
`python/`, `r/`, or `bash/` when those folders obscure the workflow stages.

## 5. Notebook Rules

- Notebooks may be used for exploration, review, and communication.
- Notebooks must not be the only record of formal analysis logic.
- Reusable notebook logic must be promoted into formal scripts.
- If a notebook materially changes project conclusions, the agent must update
  scripts or memory accordingly.

## 6. Naming Rules

- Names must be explicit and descriptive.
- The agent must not use ambiguous names such as `final`, `final2`, `new`, or
  `test`.
- File and directory names should encode analysis stage, cohort, or content.

## 7. Run ID Conventions

Formal run IDs must use the form:

`YYYY-MM-DD_HHMM_topic_vN`

Example:

`2026-03-21_1030_pbmc_qc_v1`

The agent must use one run directory per formal execution sequence.

## 8. Reproducibility and Traceability

Each formal run must have:

- a run ID
- a defined input set
- a recorded configuration
- scripts or commands that can be rerun
- outputs stored under the run directory

The agent must record promotion decisions and important run outcomes in
`memory/Status.md` or `memory/Decisions.md`.

## 9. Promotion Rules

Outputs may be promoted from `runs/` into `data/processed/` only when:

- the generating run is identifiable
- the derivation is reproducible
- the asset is expected to be reused
- the asset has a stable name and clear provenance

Outputs may be promoted into `results/` only when:

- they are curated for review, communication, or delivery
- their source run or source processed asset is traceable

## 10. Agent-Generated Artifacts

- Agent-generated working notes belong in `memory/` when they are durable.
- Agent-generated formal scripts belong in `scripts/`.
- Agent-generated run artifacts belong in `runs/`.
- Agent-generated curated deliverables belong in `results/`.

The agent must not leave important governance state only in transient chat.

## 11. Tooling and Environment Paths

Durable paths for tools, references, and environments must be stored in
`config/default.yaml`.

Examples include:

- `cellranger` executable path
- reference genome directory
- `mamba` executable path
- reusable environment names or prefixes

If the project relies on them operationally, the agent must also summarize them
in `memory/Decisions.md`.
