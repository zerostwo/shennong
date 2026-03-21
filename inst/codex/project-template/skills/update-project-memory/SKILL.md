# update-project-memory

## Purpose

Standardize how and when the project memory files are updated.

## When To Use

- after meaningful progress
- after a durable decision
- after learning a durable environment location
- after changing planned scope or sequencing

## Required Inputs

- new durable knowledge
- updated project state
- scope or plan changes

## Required Outputs

- updated `memory/Decisions.md`, `memory/Plan.md`, `memory/Prompt.md`, or
  `memory/Status.md` as appropriate

## Rules

- durable decisions must go to `memory/Decisions.md`
- current state must go to `memory/Status.md`
- intended next work must go to `memory/Plan.md`
- durable operating context must go to `memory/Prompt.md`
- tool paths, reference paths, and environment names must also be kept in
  `config/default.yaml`

## Procedure

1. Classify the new information.
2. Update the correct memory file.
3. If the information is operationally reusable, also update `config/default.yaml`.
4. Keep entries concise, dated, and explicit.

## Common Mistakes

- storing durable knowledge only in transient chat
- putting current status into decisions
- forgetting to persist tool or reference locations after successful use

## Examples

- record a validated `cellranger` path in `config/default.yaml` and
  `memory/Decisions.md`
- record the active run and blocker in `memory/Status.md`
