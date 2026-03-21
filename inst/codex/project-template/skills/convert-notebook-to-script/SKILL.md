# convert-notebook-to-script

## Purpose

Standardize how exploratory notebook logic becomes formal scripts.

## When To Use

- notebook logic is being reused
- notebook logic affects formal project conclusions
- notebook steps need reproducible execution

## Required Inputs

- source notebook
- target workflow stage
- intended formal behavior

## Required Outputs

- a formal script under `scripts/`
- notebook trimmed to exploration or narrative use

## Rules

- reusable analysis logic must move into scripts
- notebooks must not remain the only source of formal execution logic

## Procedure

1. Identify the reusable notebook steps.
2. Place them into the correct stage-oriented script location.
3. Make inputs and outputs explicit.
4. Update the notebook to reference the formal script when appropriate.
5. Record the workflow change in `memory/Status.md`.

## Common Mistakes

- copying notebook cells without parameterizing paths
- leaving critical logic only in markdown or outputs
- using notebook execution order as hidden state

## Examples

- convert clustering notebook steps into `scripts/03_cluster/run_clustering.R`
