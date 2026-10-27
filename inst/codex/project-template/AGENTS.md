# {{project_name}} Analysis Governance

Keep raw inputs immutable, derived artifacts traceable to a run, and formal
analyses reproducible from maintained scripts. Preserve user-owned work.


## Documentation entry points

- [Get started](https://zerostwo.github.io/shennong/dev/articles/get-started.html): bundled-data counts, clustering, markers, scoring, and result retrieval.
- [Parameters and results](https://zerostwo.github.io/shennong/dev/articles/parameters-and-results.html): common controls, return values, result IDs, safe reruns, group aggregation, and API migration.
- [Differential expression and enrichment](https://zerostwo.github.io/shennong/dev/articles/differential-expression.html): markers versus replicated contrasts, tested-gene backgrounds, ORA/GSEA formulas, and multiple databases.
- [Choose a backend](https://zerostwo.github.io/shennong/dev/articles/method-catalog.html): method catalog and optional runtime requirements.

Start with the bundled-data examples for basic usage; real-data articles require their stated input files. Do not assume displayed optional-backend recipes executed during an ordinary website build.

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



### Feature plotting: dense embeddings

For dense feature maps, `sn_plot_feature(..., pt_size = 0.05)` keeps high values
on top and uses zero stroke by default. Explicit `stroke_size` adds an outline;
point sizes remain consistent between raster and vector output.

For compact group-level heatmaps, use `sn_plot_heatmap(..., mode = "average")`.
The pheatmap and ComplexHeatmap backends support `cluster_rows` for groups and
`cluster_cols` for features; the default ggplot2 backend preserves input order.

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

## Styled embedding figures

For `sn_plot_dim()` / `sn_plot_feature(style = "nebula" | "glass")`, use an
existing 2D reduction with default dimensions for planar density contours, or
an existing three-dimensional reduction with explicit `dims = 1:3`. Preserve
original reductions. The browser viewer uses `interactive = TRUE`; download its camera
JSON or copy its R list and pass it as `camera` to a static plot.
`sn_get_plot_camera(viewer)` retrieves only the initial camera. Export static
ggplots with `raster_dpi = 600` and `ggsave(..., dpi = 600)`: surfaces and points
are rasterized together, labels and legends remain vector. Feature colors use
the specified assay/layer; `group_by` controls envelopes only. Envelopes are
visual density aids, not tissue boundaries. Never invent a third coordinate.

The 3D renderer uses WebGL with local shaders. Static 600 dpi capture requires
Chrome/Chromium plus chromote and png, and runs through the same shader as the
viewer. Do not fall back silently to the old flat polygon renderer. Density
envelopes use local-neighbor bandwidth estimates; do not duplicate cells to
imitate the apparent point density of a reference image.

## Verify statistical units and stored mappings

For program comparisons, use biological `sample_by` IDs. Complete donor pairs
are modeled as paired; partially paired profiles require explicit input
selection. Inspect the stored comparison's `paired` flag. For program scoring,
set the top-level `seed` for reproducibility and retrieve
`tables$metadata_columns` from `sn_get_result(object, "program_scoring", id)`
instead of guessing sanitized metadata names. Discover IDs with
`sn_list_results(object, type = "program_scoring")`.

Mean program scoring supports single cells and single aggregated groups.
Signature tables need both program and gene columns; alternatively supply a
named list of gene vectors. Single-entity scores are descriptive, not evidence
of biological replication.

Stored-DE ORA uses the finite backend-test background of each comparison,
recorded in `input$tested_features_by_comparison`; available
`input$candidate_features` are not an ORA background. Rerun older candidate-only
DE records or provide a justified explicit universe. Communication reads
matching split layers, Milo annotation filters use stored `annotation_by`,
and native RDS/RData/QS2 serialization preserves matrix classes.

DE subset labels must be complete and non-blank. Choose a `subset_by` name
that does not collide with backend columns such as `gene`, `comparison`, or
marker `cluster`; copy metadata to a distinct name if needed. Empty reported
subsets contribute zero rows without aborting other subsets. Retrieve empty
results normally, and inspect `input$tested_features_by_comparison` rather
than treating an empty marker table as proof that no genes were tested.

## Analysis API contracts

Use `batch_by` for clustering batch metadata, `backend_control` for optional
backend settings, and `n_workers` for the common parallel-worker control.
`sn_run_scvi`, `sn_run_scanvi`, and `sn_run_scpoli` are clustering shortcuts.
Core clustering `assay`, `layer`, `npcs`, `dims`, and `hvg_group_by` are explicit
parameters; advanced `...` controls must be named.

DE, enrichment, scoring, and Milo use `return_object = FALSE` for a unified
result. Read complete metadata with `sn_get_result()` and filtered tables with
typed getters. DE/enrichment getters accept `top_scope = "group"` or `"all"`,
`p_adjusted_cutoff`, and independent direction/effect filters for DE. An omitted
ID resolves only a unique result. Discover ambiguous choices with
`sn_list_results()`. Never assume a default write overwrites a prior analysis:
explicit storage replacement requires `overwrite = TRUE`.

Pass computational `seed` at the top level. Group program scoring requires
`aggregate = "expression"` (average expression, then score) or `"scores"`
(score cells, then average); record which question the analysis answers.
Milo's `keep_model` retains `models$milo` without changing its return shape.

## Styled embedding figures

For `sn_plot_dim()` / `sn_plot_feature(style = "nebula" | "glass")`, use an
existing 2D reduction with default dimensions for planar density contours, or
an existing three-dimensional reduction with explicit `dims = 1:3`. Preserve
original reductions. The browser viewer uses `interactive = TRUE`; download its camera
JSON or copy its R list and pass it as `camera` to a static plot.
`sn_get_plot_camera(viewer)` retrieves only the initial camera. Export static
ggplots with `raster_dpi = 600` and `ggsave(..., dpi = 600)`: surfaces and points
are rasterized together, labels and legends remain vector. Feature colors use
the specified assay/layer; `group_by` controls envelopes only. Envelopes are
visual density aids, not tissue boundaries. Never invent a third coordinate.

The 3D renderer uses WebGL with local shaders. Static 600 dpi capture requires
Chrome/Chromium plus chromote and png, and runs through the same shader as the
viewer. Do not fall back silently to the old flat polygon renderer. Density
envelopes use local-neighbor bandwidth estimates; do not duplicate cells to
imitate the apparent point density of a reference image.


## Spatial map presentation

Use `sn_plot_spatial()` / `sn_plot_spatial_feature()` with dimension-plot-style
`pt_size`, `cols`/`palette`, and visibility controls. Metadata coordinates need
an explicit section `split_by`; native Seurat image/FOV centroids are selected
with `images` and kept in separate panels. Preserve coordinate units and
orientation (`flip_y = FALSE` for Cartesian axes). Histology and segmentation
outlines are outside these coordinate maps. Select the exact expression assay
and layer; missing layers error and never become raw counts implicitly.
Numeric metadata maps need no expression layer. Shared `limits` affect colors
only, not exported source values. Discover stored spatial domains with
`sn_list_results(object, type = "spatial_domains")`, retrieve with
`sn_get_result(object, "spatial_domains", result_id)`, and plot with
`sn_plot_spatial_domain(object, result_id = result_id, label = TRUE)`.


### Publication-inspired atlas styles

Discover `sn_list_plot_styles()` and use opt-in `atlas`, `atlas_outline`, or
`atlas_module` only on supported plotting functions. See the
`atlas-figure-styles` article for paper provenance and the figure inventory.
Choose the statistical encoding before its appearance: expression percentage,
frequency, Ro/e, odds ratio, P value and model coefficients are distinct.
`sn_plot_matrix()` and `sn_plot_network()` render supplied results and do not
infer modules or communication. Preserve sample units, pairing, source values
and common color limits; inspect exported PDF/PNG and source-data readback.

### Structural atlas layouts

`sn_list_plot_styles(function_name = "sn_plot_matrix")` lists structural styles
and the exact paper panels/input requirements. Themes (`atlas`,
`atlas_outline`, `atlas_module`) are separate from the 23 geometry choices.
Use matrix `split_triangle`, `circular`, `annotated`, `dendrogram`,
`grouped_bubble`, `heatmap_bar`, `heatmap_density`; distribution `violin_box`,
`raincloud`, `paired_box`, `paired_slope`, `split_violin`, `stacked_violin`;
`sn_plot_effect()` for supplied coefficient/ratio forests, signed bars and
lollipops; `sn_plot_overlap()` for exclusive UpSet intersections; and
`sn_plot_network()` for bipartite, supplied clone-tree, weighted chord and
shared-coordinate rewiring structures. Existing `sn_plot_heatmap(mode="average")`
and `sn_plot_violin()` expose the matching matrix/distribution layouts.

Require explicit complete subject pairs (`pair_by`, optional specimen
`sample_by`), an explicit coefficient versus ratio scale, original interval
bounds, and supplied network/tree topology. No automatic P-values, BCR trees,
CoVarNet discovery or imputation. Reusable source tables survive PDF/PNG bundle
export. Consult `extdata/zhang-atlas-layouts.csv` and the rendered atlas article;
synthetic layout examples must not be described as original cancer results.
