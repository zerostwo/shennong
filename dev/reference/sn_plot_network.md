# Plot a supplied cell association or cellular-module network

Draws an already computed edge table. A circular layout is deterministic
and conveys no spatial distance. This renderer does not infer
co-occurrence, ligand-receptor communication or CoVarNet modules. Supply
those results explicitly, including the appropriate score label and
provenance.

## Usage

``` r
sn_plot_network(
  edges,
  nodes = NULL,
  group_by = NULL,
  label_by = NULL,
  layout = c("circle", "supplied"),
  directed = FALSE,
  style = "classic",
  palette = NULL,
  edge_palette = NULL,
  weight_label = "Weight",
  node_size = 5,
  label_size = 8,
  title = NULL,
  compare_edges = NULL,
  branch_length = NULL
)
```

## Arguments

- edges:

  Data frame containing `source`, `target` and finite numeric `weight`.
  Duplicate edges and self loops are rejected rather than aggregated.

- nodes:

  Optional data frame with unique `node` identifiers. Includes isolated
  nodes. Defaults to nodes found in the edge table.

- group_by:

  Optional node column mapped to fill (for example lineage).

- label_by:

  Optional node label column; defaults to `node`.

- layout:

  `circle` or `supplied`; supplied coordinates use node `x`, `y`.

- directed:

  Whether to draw arrows. Arrowheads are inset from target centers for
  visibility; supplied coordinates and edge values are retained.

- style:

  Presentation preset; see
  [`sn_list_plot_styles()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_plot_styles.md).

- palette:

  Optional node palette or named colors.

- edge_palette:

  Optional edge color gradient. Module style uses grey.

- weight_label:

  Legend title describing the supplied edge weight.

- node_size, label_size:

  Node size in mm and label size in points.

- title:

  Optional title.

- compare_edges:

  Second source/target/weight edge table for `rewiring`. Both panels
  share node coordinates and edge scales; absent edges stay absent.

- branch_length:

  Optional edge column containing non-negative branch lengths for
  `clone_tree`; otherwise depth uses unit edges. Topology is supplied,
  never inferred from BCR sequences or the `weight` score.

## Value

A ggplot or patchwork retaining supplied node/edge tables. Structural
styles: `bipartite` requires exactly two node groups and cross-group
edges; `clone_tree` requires a connected rooted tree with one parent per
node; `chord` draws ribbons whose endpoint widths encode non-negative
weights; `rewiring` draws two networks at identical coordinates.
Directed chords are unsupported; use a directed network to retain
arrows.

## Examples

``` r
edges <- data.frame(source = c("NK", "B"), target = c("DC", "DC"),
                    weight = c(0.4, 0.8))
sn_plot_network(edges, style = "atlas_module", weight_label = "Specificity")
```
