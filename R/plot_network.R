.sn_network_styles <- function() c("bipartite", "clone_tree", "chord", "rewiring")

.sn_plot_network_layout <- function(args) {
  style <- args$style; edges <- args$edges; nodes <- args$nodes
  base <- args; base$style <- "atlas"; base$compare_edges <- base$branch_length <- NULL
  if (style != "rewiring" && !is.null(args$compare_edges)) stop("compare_edges requires style = 'rewiring'.", call. = FALSE)
  if (style != "clone_tree" && !is.null(args$branch_length)) stop("branch_length requires style = 'clone_tree'.", call. = FALSE)
  if (style %in% c("bipartite", "clone_tree") && is.null(base$edge_palette) && all(edges$weight >= 0)) base$edge_palette <- c("#B5C0C6", "#344955")
  if (style == "bipartite") {
    if (is.null(args$group_by)) stop("bipartite requires group_by for the two node partitions.", call. = FALSE)
    g <- as.character(nodes[[args$group_by]]); groups <- unique(g)
    if (anyNA(g) || length(groups) != 2L) stop("bipartite requires exactly two non-missing node groups.", call. = FALSE)
    if (any(g[match(edges$source, nodes$node)] == g[match(edges$target, nodes$node)])) stop("Bipartite edges must cross the two groups.", call. = FALSE)
    nodes$x <- match(g, groups); nodes$y <- 0
    for (level in groups) nodes$y[g == level] <- seq(1, 0, length.out = sum(g == level))
    base$nodes <- nodes; base$layout <- "supplied"
    plot <- do.call(sn_plot_network, base)
  } else if (style == "clone_tree") {
    ids <- as.character(nodes$node); targets <- as.character(edges$target)
    roots <- setdiff(ids, targets)
    if (length(roots) != 1L || anyDuplicated(targets) || nrow(edges) != nrow(nodes) - 1L) stop("clone_tree needs one connected rooted tree and one parent per non-root node.", call. = FALSE)
    lengths <- rep(1, nrow(edges))
    if (!is.null(args$branch_length)) {
      if (length(args$branch_length) != 1L || !args$branch_length %in% names(edges)) stop("branch_length must name an edge column.", call. = FALSE)
      lengths <- edges[[args$branch_length]]
      if (!is.numeric(lengths) || any(!is.finite(lengths)) || any(lengths < 0)) stop("Branch lengths must be finite and non-negative.", call. = FALSE)
    }
    nodes$x <- nodes$y <- NA_real_; names_x <- stats::setNames(rep(NA_real_, length(ids)), ids)
    names_x[roots] <- 0; visited <- character(); leaf <- 0
    visit <- function(id) {
      if (id %in% visited) stop("Cycle in clone tree.", call. = FALSE)
      visited <<- c(visited, id); rows <- which(as.character(edges$source) == id)
      if (!length(rows)) { leaf <<- leaf + 1; nodes$y[match(id, ids)] <<- leaf; return(leaf) }
      yy <- vapply(rows, function(i) {
        child <- as.character(edges$target[i]); names_x[child] <<- names_x[id] + lengths[i]; visit(child)
      }, numeric(1))
      nodes$y[match(id, ids)] <<- mean(range(yy)); mean(range(yy))
    }
    visit(roots)
    if (length(visited) != length(ids)) stop("Disconnected node or cycle: input is not a rooted tree.", call. = FALSE)
    nodes$x <- unname(names_x[ids]); base$nodes <- nodes; base$layout <- "supplied"; base$directed <- FALSE
    plot <- do.call(sn_plot_network, base)
    links <- plot$layers[[1]]$data
    elbow <- links; elbow$xend <- elbow$x
    horizontal <- links; horizontal$y <- horizontal$yend
    plot$layers[[1]]$data <- rbind(elbow, horizontal)
    plot <- plot + ggplot2::labs(caption = if (is.null(args$branch_length)) "Branch x distance: unit edge depth" else paste("Branch x distance:", args$branch_length))
  } else if (style == "chord") {
    if (args$directed) stop("chord does not encode arrows; directed = TRUE is unsupported.", call. = FALSE)
    if (any(edges$weight < 0) || !any(edges$weight > 0)) stop("Chord ribbons need non-negative weights, with at least one positive edge.", call. = FALSE)
    ids <- as.character(nodes$node)
    total <- vapply(ids, function(id) sum(edges$weight[edges$source == id | edges$target == id]), numeric(1))
    if (any(total == 0)) stop("Chord nodes need positive incident weight; remove isolated/zero-weight nodes explicitly.", call. = FALSE)
    gap <- min(0.04, pi / length(ids)); available <- 2 * pi - gap * length(ids)
    span <- total / sum(total) * available
    start <- c(0, utils::head(cumsum(span + gap), -1)); cursor <- start
    arc <- function(a, b, r = 1) { theta <- seq(a, b, length.out = 24); cbind(x = r * cos(theta), y = r * sin(theta)) }
    bezier <- function(a, b) { t <- seq(0, 1, length.out = 32); cbind(x = (1 - t)^2 * cos(a) + t^2 * cos(b), y = (1 - t)^2 * sin(a) + t^2 * sin(b)) }
    ribbons <- list(); allocations <- list()
    for (i in which(edges$weight > 0)) {
      s <- match(edges$source[i], ids); t <- match(edges$target[i], ids)
      a <- cursor[s]; b <- cursor[t]; aw <- edges$weight[i] / total[s] * span[s]; bw <- edges$weight[i] / total[t] * span[t]
      xy <- rbind(arc(a, a + aw), bezier(a + aw, b), arc(b, b + bw), bezier(b + bw, a))
      ribbons[[length(ribbons) + 1L]] <- data.frame(xy, edge = i, category = if (is.null(args$group_by)) ids[s] else as.character(nodes[[args$group_by]][s]))
      allocations[[length(allocations) + 1L]] <- data.frame(edge = i, source_start = a, source_end = a + aw, target_start = b, target_end = b + bw, weight = edges$weight[i])
      cursor[s] <- a + aw; cursor[t] <- b + bw
    }
    sectors <- do.call(rbind, lapply(seq_along(ids), function(i) data.frame(rbind(arc(start[i], start[i] + span[i], 1.02), arc(start[i] + span[i], start[i], 1.1)),
      node = ids[i], category = if (is.null(args$group_by)) ids[i] else as.character(nodes[[args$group_by]][i]))))
    mid <- start + span / 2; nodes$x <- 1.23 * cos(mid); nodes$y <- 1.23 * sin(mid)
    nodes$.label <- as.character(nodes[[args$label_by %||% "node"]])
    ribbon_data <- do.call(rbind, ribbons)
    plot <- ggplot2::ggplot() + ggplot2::geom_polygon(data = ribbon_data, ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge, fill = .data$category), alpha = 0.55, colour = NA) +
      ggplot2::geom_polygon(data = sectors, ggplot2::aes(x = .data$x, y = .data$y, group = .data$node, fill = .data$category), colour = "white", linewidth = 0.2) +
      ggplot2::geom_text(data = nodes, ggplot2::aes(x = .data$x, y = .data$y, label = .data$.label), size = args$label_size / ggplot2::.pt) +
      ggplot2::coord_equal(xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45), clip = "off") + ggplot2::theme_void(base_size = 8) +
      ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8), plot.caption = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(title = args$title, fill = args$group_by %||% "Source node", caption = paste("Ribbon endpoint width:", args$weight_label))
    plot <- .sn_add_discrete_palette(plot, args$palette %||% .sn_atlas_palette("atlas"), length(unique(sectors$category)), aesthetic = "fill")
  } else {
    if (is.null(args$compare_edges)) stop("rewiring requires compare_edges.", call. = FALSE)
    first <- do.call(sn_plot_network, base)
    base$edges <- args$compare_edges; base$nodes <- merge(nodes, attr(first, "shennong_figure_data")$coordinates, by = "node", sort = FALSE, suffixes = c(".old", ""))
    base$nodes <- base$nodes[match(nodes$node, base$nodes$node), , drop = FALSE]; base$layout <- "supplied"
    second <- do.call(sn_plot_network, base)
    weights <- c(edges$weight, args$compare_edges$weight); magnitude <- max(abs(weights), 1e-10)
    fix_scale <- function(p) suppressMessages(p + ggplot2::scale_linewidth_continuous(limits = c(0, magnitude), range = c(0.2, 1), name = paste0("|", args$weight_label, "|")) +
      ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-magnitude, magnitude), name = args$weight_label))
    plot <- patchwork::wrap_plots(fix_scale(first) + ggplot2::labs(title = "Reference network"), fix_scale(second) + ggplot2::labs(title = "Comparison network"), nrow = 1, guides = "collect")
    nodes$x <- attr(first, "shennong_figure_data")$coordinates$x; nodes$y <- attr(first, "shennong_figure_data")$coordinates$y
  }
  plot <- .sn_style_plot(plot, "atlas", "network"); attr(plot, "shennong_plot_style")$name <- style
  sources <- list(edges = edges, nodes = args$nodes, coordinates = nodes[, c("node", "x", "y")])
  if (style == "chord") sources$ribbon_endpoints <- do.call(rbind, allocations)
  if (style == "rewiring") sources$compare_edges <- args$compare_edges
  .sn_attach_figure_spec(plot, "network", list(n_nodes = nrow(nodes), n_edges = nrow(edges), n_points = nrow(nodes)), source_data = sources)
}

#' Plot a supplied cell association or cellular-module network
#'
#' Draws an already computed edge table. A circular layout is deterministic and
#' conveys no spatial distance. This renderer does not infer co-occurrence,
#' ligand-receptor communication or CoVarNet modules. Supply those results
#' explicitly, including the appropriate score label and provenance.
#' @param edges Data frame containing `source`, `target` and finite numeric
#'   `weight`. Duplicate edges and self loops are rejected rather than aggregated.
#' @param nodes Optional data frame with unique `node` identifiers. Includes
#'   isolated nodes. Defaults to nodes found in the edge table.
#' @param group_by Optional node column mapped to fill (for example lineage).
#' @param label_by Optional node label column; defaults to `node`.
#' @param layout `circle` or `supplied`; supplied coordinates use node `x`, `y`.
#' @param directed Whether to draw arrows. Arrowheads are inset from target
#'   centers for visibility; supplied coordinates and edge values are retained.
#' @param style Presentation preset; see [sn_list_plot_styles()].
#' @param palette Optional node palette or named colors.
#' @param edge_palette Optional edge color gradient. Module style uses grey.
#' @param weight_label Legend title describing the supplied edge weight.
#' @param node_size,label_size Node size in mm and label size in points.
#' @param title Optional title.
#' @param compare_edges Second source/target/weight edge table for `rewiring`.
#'   Both panels share node coordinates and edge scales; absent edges stay absent.
#' @param branch_length Optional edge column containing non-negative branch
#'   lengths for `clone_tree`; otherwise depth uses unit edges. Topology is
#'   supplied, never inferred from BCR sequences or the `weight` score.
#' @return A ggplot or patchwork retaining supplied node/edge tables. Structural
#'   styles: `bipartite` requires exactly two node groups and cross-group edges;
#'   `clone_tree` requires a connected rooted tree with one parent per node;
#'   `chord` draws ribbons whose endpoint widths encode non-negative weights;
#'   `rewiring` draws two networks at identical coordinates. Directed chords are
#'   unsupported; use a directed network to retain arrows.
#' @examples
#' edges <- data.frame(source = c("NK", "B"), target = c("DC", "DC"),
#'                     weight = c(0.4, 0.8))
#' sn_plot_network(edges, style = "atlas_module", weight_label = "Specificity")
#' @export
sn_plot_network <- function(edges, nodes = NULL, group_by = NULL, label_by = NULL,
                            layout = c("circle", "supplied"), directed = FALSE,
                            style = "classic", palette = NULL, edge_palette = NULL,
                            weight_label = "Weight", node_size = 5, label_size = 8,
                            title = NULL, compare_edges = NULL, branch_length = NULL) {
  args <- as.list(environment())
  style <- .sn_match_plot_style(style, .sn_network_styles()); layout <- match.arg(layout)
  if (!is.data.frame(edges) || !all(c("source", "target", "weight") %in% names(edges))) stop("`edges` requires source, target and weight columns.", call. = FALSE)
  if (!is.numeric(edges$weight) || any(!is.finite(edges$weight))) stop("Edge weights must be finite numeric values.", call. = FALSE)
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) stop("`directed` must be TRUE or FALSE.", call. = FALSE)
  for (arg in c("node_size", "label_size")) {
    v <- get(arg); if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) stop("`", arg, "` must be positive.", call. = FALSE)
  }
  source <- as.character(edges$source); target <- as.character(edges$target)
  if (anyNA(source) || anyNA(target) || any(!nzchar(source)) || any(!nzchar(target))) stop("Edge identifiers cannot be missing or empty.", call. = FALSE)
  if (any(source == target)) stop("Self loops require a separate rendering; remove them explicitly.", call. = FALSE)
  pairs <- if (directed) data.frame(source, target) else data.frame(source = pmin(source, target), target = pmax(source, target))
  if (anyDuplicated(pairs)) stop("Duplicate edges must be resolved explicitly before plotting.", call. = FALSE)
  if (is.null(nodes)) nodes <- data.frame(node = sort(unique(c(source, target))))
  if (!is.data.frame(nodes) || !"node" %in% names(nodes) || !nrow(nodes) || anyNA(nodes$node) || any(!nzchar(as.character(nodes$node))) || anyDuplicated(as.character(nodes$node))) stop("`nodes` requires unique non-empty node identifiers.", call. = FALSE)
  if (length(setdiff(c(source, target), as.character(nodes$node)))) stop("Edge endpoints are missing from `nodes`.", call. = FALSE)
  if (length(setdiff(c(group_by, label_by), names(nodes)))) stop("Node grouping or label columns are missing.", call. = FALSE)
  if (style %in% .sn_network_styles()) {
    args$nodes <- nodes; args$style <- style
    return(.sn_plot_network_layout(args))
  }
  if (!is.null(compare_edges) || !is.null(branch_length)) stop("compare_edges and branch_length require their structural network styles.", call. = FALSE)
  data <- nodes
  if (layout == "circle") {
    theta <- seq(0, 2 * pi, length.out = nrow(nodes) + 1L)[seq_len(nrow(nodes))]
    data$x <- cos(theta); data$y <- sin(theta)
  } else if (!all(c("x", "y") %in% names(nodes)) || !is.numeric(nodes$x) || !is.numeric(nodes$y) || any(!is.finite(nodes$x)) || any(!is.finite(nodes$y))) stop("Supplied layout needs finite numeric node x and y columns.", call. = FALSE)
  data$.sn_group <- if (is.null(group_by)) "All" else as.character(nodes[[group_by]])
  data$.sn_label <- as.character(nodes[[label_by %||% "node"]])
  links <- edges
  links$x <- data$x[match(source, nodes$node)]; links$y <- data$y[match(source, nodes$node)]
  links$xend <- data$x[match(target, nodes$node)]; links$yend <- data$y[match(target, nodes$node)]
  if (directed && nrow(links)) {
    # Inset arrowheads from node centers so the node glyph cannot hide them.
    dx <- links$xend - links$x; dy <- links$yend - links$y
    distance <- sqrt(dx^2 + dy^2)
    span <- max(diff(range(data$x)), diff(range(data$y)))
    inset <- pmin(distance * 0.25, span * 0.04 * node_size / 5)
    fraction <- ifelse(distance > 0, inset / distance, 0)
    links$xend <- links$xend - dx * fraction
    links$yend <- links$yend - dy * fraction
  }
  if (is.null(palette)) palette <- .sn_atlas_palette(style)
  signed_default <- is.null(edge_palette) && style != "atlas_module" && any(edges$weight < 0)
  if (is.null(edge_palette)) edge_palette <- if (style == "atlas_module") c("#DDDDDD", "#333333") else c("#2166AC", "#F7F7F7", "#B2182B")
  plot <- ggplot2::ggplot(data) +
    ggplot2::geom_segment(data = links, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, colour = .data$weight, linewidth = abs(.data$weight)),
      arrow = if (directed) grid::arrow(length = grid::unit(2, "mm")) else NULL, alpha = 0.8) +
    ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y, fill = .data$.sn_group), shape = 21, size = node_size, stroke = 0.3, colour = "white") +
    ggplot2::geom_text(ggplot2::aes(x = .data$x, y = .data$y, label = .data$.sn_label), size = label_size / ggplot2::.pt, vjust = -1.1) +
    ggplot2::scale_linewidth_continuous(range = c(0.2, 1), name = paste0("|", weight_label, "|")) +
    ggplot2::coord_equal(clip = "off") + ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(plot.margin = ggplot2::margin(16, 24, 16, 24)) + ggplot2::labs(title = title, fill = group_by)
  edge_scale <- if (signed_default) ggplot2::scale_colour_gradient2(
    low = edge_palette[1], mid = edge_palette[2], high = edge_palette[3], midpoint = 0, name = weight_label) else
    ggplot2::scale_colour_gradientn(colours = edge_palette, name = weight_label)
  plot <- plot + edge_scale
  plot <- .sn_add_discrete_palette(plot, palette, length(unique(data$.sn_group)), aesthetic = "fill")
  if (is.null(group_by)) plot <- plot + ggplot2::guides(fill = "none")
  .sn_attach_figure_spec(.sn_style_plot(plot, style, "network"), "network",
    list(n_points = nrow(nodes), n_nodes = nrow(nodes), n_edges = nrow(edges),
      n_groups = length(unique(data$.sn_group)), labels = data$.sn_label),
    source_data = list(edges = edges, nodes = nodes, coordinates = data[, c("node", "x", "y")], directed = directed))
}
