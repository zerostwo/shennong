# {{project_name}} Analysis Governance

Keep raw inputs immutable, derived artifacts traceable to a run, and formal
analyses reproducible from maintained scripts. Preserve user-owned work.

## Task-specific context

- Use `docs/standards/BioinformaticsAnalysisConventions.md` when creating runs,
  placing outputs, organizing scripts, or promoting data.
- Read relevant `memory/Decisions.md` entries for existing scientific choices;
  use `memory/Status.md` and `memory/Plan.md` when resuming work.
- Consult `memory/Prompt.md` for durable operating context and
  `config/default.yaml` for validated environments and reference locations.
- Use the matching project skill under `skills/` when its workflow is needed.
  A small edit does not require reading every governance record.

## Execution and outputs

Use explicit run IDs and meaningful stage names. Keep source data, derived data,
run artifacts, and curated exports separate. Never overwrite `data/raw/` or
export unchecked intermediates as final results. Notebook output alone is not
the formal record; retain the script, inputs, parameters, and validation evidence.

Continue authorized analysis through the requested deliverables and relevant
validation. Reuse existing authorization for ordinary run bookkeeping and local
repairs. Ask for material scientific scope changes or operations requiring new
authority. Report incomplete execution and validation accurately.

Update only the project records whose durable information changed: decisions in
`memory/Decisions.md`, current state in `memory/Status.md`, scope/sequencing in
`memory/Plan.md`, and operating context in `memory/Prompt.md`. Store reusable
environment details in `config/default.yaml`. These are project files; any
separate agent-managed memory system follows its own write authorization.

## Shennong discovery and usage tracking

Install package guidance with `sn_install_codex_skill(type = "package_skills")`
when needed. `sn_get_mcp_server_config()` provides read-only method/help discovery;
keep analysis execution in explicit project scripts.

- If Shennong usage tracking is enabled, place its local SQLite database below
  ignored `runs/`, select development/production/test/benchmark mode
  explicitly, and never commit or upload the raw SQLite activity log. Enable
  before caching Shennong function references; the current namespace-binding
  tracker does not intercept references saved earlier. Remote DBI delivery must
  retain the local outbox and requires explicit versioned research consent and
  a sanitized flush.



### Composition plotting: dedicated string interfaces

Use `sn_plot_sankey(data, flow_by = c("level1", "level2"),
fill_by = "level1", style = "A")` for hierarchies. B is bracketed two-stage
annotation comparison (default final-label colors); C is vertical two-stage
composition with source-share pies (default first-stage colors). C pies use the
same retained weights within each target/facet; cross-facet widths are relative.
`show_pies = FALSE` omits pies. Factor levels control each stage's order.
Use `sn_plot_bar`, `sn_plot_sample_bar`, `sn_plot_sample_boxplot`, and
`sn_plot_histogram` for other composition types. Supply `x_by`, `y_by`,
`fill_by`, `facet_row_by`, and `facet_col_by` as column-name strings. Sample
summaries also require `sample_by`. All share `sn_plot_composition` preparation.
Legacy unquoted mappings remain accepted only by the unified function; never
mix a legacy mapping with its string replacement. See the composition-analysis
and visualization articles for runnable examples and styling controls.

Composition panel sizes use `panel_widths`/`panel_heights` in pt for all five
entry points. A/B support facet-size vectors; C requires scalars and preserves
circular pies. These are panel dimensions, not the full export canvas.

Composition plots default to 8 pt for labels, titles, axes, strips, and legends.
Use `sn_plot_sankey(..., style = "C", show_pies = FALSE)` to hide bottom pies;
`show_pies = TRUE` restores them without changing ribbon counts.
