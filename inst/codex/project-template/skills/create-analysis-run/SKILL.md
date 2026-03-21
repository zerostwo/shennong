# create-analysis-run

## Purpose

Create a new formal analysis run under `runs/` using the project run-ID
convention and traceable inputs.

## When To Use

- starting a formal analysis execution
- rerunning a workflow with changed parameters
- creating a reproducible run record for review

## Required Inputs

- analysis topic
- intended scripts or commands
- input datasets or processed assets
- relevant configuration

## Required Outputs

- a new run directory under `runs/`
- a run ID using `YYYY-MM-DD_HHMM_topic_vN`
- a short run manifest or note describing inputs and intent

## Rules

- must use the run-ID convention
- must not reuse an old run directory for a materially different execution
- must record enough information to reproduce the run

## Procedure

1. Determine the analysis topic and next version number.
2. Create `runs/<run_id>/`.
3. Record inputs, scripts, and config references for the run.
4. Execute formal scripts against that run context.
5. Update `memory/Status.md` with the run ID and purpose.

## Common Mistakes

- using ambiguous run names
- mixing outputs from multiple runs
- failing to record which scripts produced the outputs

## Examples

- `runs/2026-03-21_1030_pbmc_qc_v1/`
- `runs/2026-03-21_1415_crc_de_v2/`
