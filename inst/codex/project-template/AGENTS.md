# {{project_name}} Analysis Governance

This repository is a formal bioinformatics analysis project initialized from the
Shennong package template. The agent must treat this repository as a governed
analysis project, not as an ad hoc scratch directory.

## Operating Priorities

1. Preserve reproducibility, traceability, and auditability.
2. Keep raw inputs immutable and derived artifacts attributable to a run.
3. Prefer formal scripts and recorded runs over notebook-only work.
4. Persist durable project knowledge so future sessions can reuse it safely.
5. Keep outputs organized so another analyst can re-run and review the work.

## Hard Constraints

- The agent must read the required governance files before starting formal work.
- The agent must use the directory conventions defined in
  `docs/standards/BioinformaticsAnalysisConventions.md`.
- The agent must record durable project decisions in `memory/Decisions.md`.
- The agent must record current work state in `memory/Status.md`.
- The agent must update `memory/Plan.md` when scope or sequencing changes.
- The agent must update `memory/Prompt.md` when the durable operating context
  changes.
- The agent must record durable environment details such as executable paths,
  reference locations, and reusable environment names in `config/default.yaml`
  and, when they become project assumptions, in `memory/Decisions.md`.
- The agent must use `skills/` for standard operating procedures inside this
  initialized project.

## The Agent Must Always Do

- Read governance in the required order before formal analysis work.
- Keep source data, derived data, runs, and exported results separated.
- Use explicit names, run identifiers, and dated artifacts.
- Promote reusable outputs only through documented promotion rules.
- Capture enough context so a later session can continue without guesswork.

## The Agent Must Never Do

- Must not overwrite `data/raw/` source files.
- Must not place formal outputs in ambiguous locations such as the repository
  root.
- Must not use ambiguous names such as `final`, `final2`, `new`, or `test`.
- Must not treat notebook output as the only record of formal analysis work.
- Must not silently change directory conventions or project rules.

## Reading Priority Order

1. `AGENTS.md`
2. `docs/standards/BioinformaticsAnalysisConventions.md`
3. `memory/Decisions.md`
4. `memory/Plan.md`
5. `memory/Status.md`
6. `memory/Prompt.md`

## Pointers

- Standards: `docs/standards/BioinformaticsAnalysisConventions.md`
- Memory: `memory/`
- Project skills: `skills/`
- Package usage skills can be installed separately from Shennong with
  `sn_install_codex_skill(type = "package_skills")` when API-level guidance is
  needed.

## Do / Do Not

Do:
- use explicit run IDs
- keep scripts stage-oriented
- persist durable knowledge after using it successfully

Do not:
- hide important paths only in transient chat
- export unchecked intermediates as final results
- invent new conventions without recording them
