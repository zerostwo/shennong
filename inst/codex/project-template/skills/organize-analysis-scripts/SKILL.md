# organize-analysis-scripts

## Purpose

Keep formal scripts structured by analysis stage and maintain a readable,
traceable workflow layout.

## When To Use

- adding new formal scripts
- reorganizing stage folders
- promoting exploratory logic into maintained scripts

## Required Inputs

- workflow stages
- script responsibilities
- naming conventions from project standards

## Required Outputs

- scripts placed in stage-oriented directories
- explicit script names that describe their role

## Rules

- scripts must be organized by analysis stage, not by language alone
- script names must be explicit
- reusable logic must move out of notebooks into scripts

## Procedure

1. Identify the workflow stage.
2. Place or create the script in the matching stage area.
3. Use an explicit filename that reflects the action.
4. Update `memory/Status.md` if the script changes the formal workflow.

## Common Mistakes

- storing scripts in generic `misc` or `temp` folders
- naming scripts `run.R`, `test.py`, or `final.sh`
- keeping critical logic only in notebooks

## Examples

- `scripts/02_qc/filter_cells.R`
- `scripts/04_de/run_cluster_markers.R`
