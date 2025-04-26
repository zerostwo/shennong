#' @export
sn_plot_boxplot <- function(data, x, y, sort = FALSE) {
  ggplot(data = data, mapping = aes(x = {{ x }}, y = {{ y }})) +
    geom_boxplot(
      outlier.shape = NA,
      staplewidth = 0.2,
      fatten = 1
    )
}

#' @export
sn_plot_barplot <- function(data, x, y, fill, sort_by = NULL) {
  ggplot(
    data = data,
    mapping = aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})
  ) +
    geom_col()
}

#' Create a dimensionality reduction plot for categorical data
#'
#' This function creates a dimensionality reduction plot for categorical data using Seurat and ggplot2. It allows for the selection of the reduction method, grouping, and splitting variables, as well as the visualization of labels, rasterization, and color palette. The sn_plot_dim() function is intended to be used as a wrapper around Seurat's DimPlot() function.
#'
#' @param object A Seurat object containing categorical data.
#' @param dims The dimensions to plot. Default is c(1, 2).
#' @param cells The cells to plot. Default is NULL.
#' @param cols The columns to plot. Default is NULL.
#' @param pt_size The size of the points on the plot. Default is 1.
#' @param reduction The dimensionality reduction method. Default is NULL.
#' @param group_by The variable to group data by. Default is NULL.
#' @param split_by The variable to split data by. Default is NULL.
#' @param shape_by The variable to shape data by. Default is NULL.
#' @param order The order to plot the data in. Default is NULL.
#' @param shuffle Logical value indicating whether to shuffle the data before plotting. Default is FALSE.
#' @param seed The random seed to use for shuffling the data. Default is 1.
#' @param label Logical value indicating whether to show labels on the plot. Default is FALSE.
#' @param label_size The size of the labels on the plot. Default is 8 * 0.36.
#' @param label_color The color of the labels on the plot. Default is "black".
#' @param label_box Logical value indicating whether to show a box around the labels on the plot. Default is FALSE.
#' @param repel Logical value indicating whether to use point repulsion to avoid overlapping labels. Default is TRUE.
#' @param cells_highlight The cells to highlight on the plot. Default is NULL.
#' @param cols_highlight The columns to highlight on the plot. Default is NULL.
#' @param sizes_highlight The sizes to highlight on the plot. Default is NULL.
#' @param na_value The value to use for missing data. Default is "grey50".
#' @param ncol The number of columns to use for the plot. Default is NULL.
#' @param combine Logical value indicating whether to combine the plots into a single plot. Default is TRUE.
#' @param raster Logical value indicating whether to use rasterization for improved performance. Default is TRUE.
#' @param raster_dpi The DPI to use for rasterization. Default is c(512, 512).
#' @param show_legend Logical value indicating whether to show the legend on the plot. Default is TRUE.
#' @param show_axis Logical value indicating whether to show the axis on the plot. Default is TRUE.
#' @param show_border Logical value indicating whether to show the panel and axis borders on the plot. Default is TRUE.
#' @param title The title for the plot. Default is NULL.
#' @param palette The color palette to use for the plot. Default is "Paired".
#' @param ... Additional parameters to be passed to the DimPlot() function in Seurat.
#' @return A ggplot2 object containing the dimensionality reduction plot.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(ggplot2)
#' data("pbmc_small")
#' pbmc_small <- CreateSeuratObject(counts = pbmc_small)
#' Idents(pbmc_small) <- pbmc_small$percent.mt > 10
#' sn_plot_dim(object = pbmc_small, reduction = "tsne", group_by = "ident", palette = "Set1")
#' }
#'
#' @importFrom ggplot2 theme element_blank element_text margin scale_color_manual labs
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
sn_plot_dim <- function(
    object,
    dims = c(1, 2),
    cells = NULL,
    cols = NULL,
    pt_size = 2,
    reduction = NULL,
    group_by = NULL,
    split_by = NULL,
    shape_by = NULL,
    order = NULL,
    shuffle = FALSE,
    seed = 717,
    label = FALSE,
    label_size = 8 * 0.36,
    label_color = "black",
    label_box = FALSE,
    repel = FALSE,
    cells_highlight = NULL,
    cols_highlight = "#DE2D26",
    sizes_highlight = 1,
    na_value = "grey50",
    ncol = NULL,
    combine = TRUE,
    raster = TRUE,
    raster_dpi = c(512, 512),
    show_legend = TRUE,
    show_axis = TRUE,
    show_border = TRUE,
    title = NULL,
    palette = "Paired",
    ...) {
  p <- Seurat::DimPlot(
    object = object,
    dims = dims,
    cells = cells,
    cols = cols,
    pt.size = pt_size,
    reduction = reduction,
    group.by = group_by,
    split.by = split_by,
    shape.by = shape_by,
    order = order,
    shuffle = shuffle,
    seed = seed,
    label = label,
    label.size = label_size,
    label.color = label_color,
    label.box = label_box,
    repel = repel,
    cells.highlight = cells_highlight,
    cols.highlight = cols_highlight,
    sizes.highlight = sizes_highlight,
    na.value = na_value,
    ncol = ncol,
    combine = combine,
    raster = raster,
    raster.dpi = raster_dpi,
    ...
  )
  if (!show_legend) {
    p <- p + Seurat::NoLegend()
  }
  reduction <- reduction %||% SeuratObject::DefaultDimReduc(object = object)
  object[["ident"]] <- Seurat::Idents(object = object)
  group_by <- group_by %||% "ident"
  n <- length(table(object[[group_by]]))
  x <- ifelse(reduction == "tsne", "tSNE 1",
    paste0(stringr::str_to_upper(reduction), " 1")
  )
  y <- ifelse(reduction == "tsne", "tSNE 2",
    paste0(stringr::str_to_upper(reduction), " 2")
  )
  p <- p +
    catplot::theme_cat(aspect_ratio = 1) +
    theme(legend.margin = margin(l = -8)) +
    labs(
      x = x,
      y = y,
      title = title
    )

  p <- add_palette(p = p, palette = palette, n = n)


  if (!show_axis) {
    p <- p + theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  }
  if (!show_border) {
    p <- p + theme(
      panel.border = element_blank(),
      axis.line = element_blank()
    )
  }
  return(p)
}

#' Plot a violin plot with categorical groups
#'
#' This function plots a violin plot with categorical groups using the VlnPlot function from the Seurat package.
#'
#' @param object A Seurat object containing the data to plot.
#' @param features A character vector of feature names to plot.
#' @param pt_size The size of the points to plot.
#' @param sort Whether to sort the features by their mean expression or not.
#' @param group_by A character vector specifying the grouping variable. Defaults to "ident".
#' @param split_by A character vector specifying the splitting variable.
#' @param show_legend Whether to show the legend or not. Defaults to FALSE.
#' @param angle_x The angle of the x-axis labels. Defaults to 0.
#' @param aspect_ratio The aspect ratio of the plot. Defaults to 0.5.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 expansion guide_axis scale_fill_manual
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#'
#' @examples
#' \dontrun{
#' sn_plot_violin(object = mySeuratObject, features = c("CD3D", "CD8A", "CD4"))
#' }
#'
#' @export
sn_plot_violin <- function(object,
                           features,
                           pt_size = 0,
                           sort = FALSE,
                           group_by = NULL,
                           split_by = NULL,
                           show_legend = FALSE,
                           angle_x = 0,
                           palette = "Paired",
                           aspect_ratio = 0.5) {
  p <- Seurat::VlnPlot(
    object = object,
    features = features,
    pt.size = pt_size,
    sort = sort,
    group.by = group_by,
    split.by = split_by
  )
  if (!show_legend) {
    p <- p + Seurat::NoLegend()
  }
  object[["ident"]] <- Seurat::Idents(object = object)
  group_by <- group_by %||% "ident"
  n <- length(table(object[[group_by]]))
  if (n > 12) {
    values <- colorRampPalette(brewer.pal(12, palette))(n)
  } else {
    values <- brewer.pal(n, palette)
  }
  p <- p + catplot::theme_cat(aspect_ratio = aspect_ratio) +
    theme(
      axis.title.x = element_blank(),
      legend.margin = margin(l = -8)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    guides(x = guide_axis(angle = angle_x)) +
    scale_fill_manual(values = values)
  return(p)
}

#' Plot a dot plot with categorical groups
#'
#' This function plots a dot plot with categorical groups using the DotPlot function from the Seurat package.
#'
#' @param x A Seurat object containing the data to plot.
#' @param assay The assay to plot. Defaults to the default assay.
#' @param features A character vector of feature names to plot.
#' @param col_min The minimum value for the color scale. Defaults to -2.5.
#' @param col_max The maximum value for the color scale. Defaults to 2.5.
#' @param dot_min The minimum value for the dot size scale. Defaults to 0.
#' @param dot_scale The size of the dots to plot.
#' @param idents A character vector of identities to plot. Defaults to all identities.
#' @param group_by A character vector specifying the grouping variable. Defaults to "ident".
#' @param split_by A character vector specifying the splitting variable.
#' @param cluster_idents Whether to cluster the identities or not. Defaults to FALSE.
#' @param scale Whether to scale the dot size or not. Defaults to TRUE.
#' @param scale_by The variable to scale the dot size by. Defaults to "radius".
#' @param scale_min The minimum value for the dot size scale. Defaults to NA.
#' @param scale_max The maximum value for the dot size scale. Defaults to NA.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 coord_fixed guides theme element_text margin scale_color_distiller scale_y_discrete guide_axis guide_legend guide_colorbar
#'
#' @examples
#' \dontrun{
#' sn_plot_dot(x = mySeuratObject, features = c("CD3D", "CD8A", "CD4"))
#' }
#'
#' @export
sn_plot_dot <- function(x,
                        assay = NULL,
                        features,
                        col_min = -2.5,
                        col_max = 2.5,
                        dot_min = 0,
                        dot_scale = 4,
                        idents = NULL,
                        group_by = NULL,
                        split_by = NULL,
                        cluster_idents = FALSE,
                        scale = TRUE,
                        scale_by = "radius",
                        scale_min = NA,
                        scale_max = NA,
                        palette = "RdBu") {
  p <-
    Seurat::DotPlot(
      object = x,
      assay = assay,
      features = features,
      col.min = col_min,
      col.max = col_max,
      dot.min = dot_min,
      dot.scale = dot_scale,
      idents = idents,
      group.by = group_by,
      split.by = split_by,
      cluster.idents = cluster_idents,
      scale = TRUE,
      scale.by = scale_by,
      scale.min = scale_min,
      scale.max = scale_max
    ) + coord_fixed() +
    guides(
      x = guide_axis(angle = 90),
      size = guide_legend(title = "Percent (%)"),
      color = guide_colorbar(
        title = "Z score",
        frame.colour = "black",
        frame.linewidth = 0.2,
        ticks.colour = "black",
        ticks.linewidth = 0.2
      )
    ) +
    catplot::theme_cat(show_title = "none") +
    theme(
      axis.text.x = element_text(face = "italic"),
      axis.title = element_blank(),
      legend.margin = margin(l = -8),
      panel.grid = element_line(
        colour = "lightgrey",
        linewidth = 0.2 / 1.07
      )
    ) +
    scale_color_distiller(palette = palette) +
    scale_y_discrete(limits = rev)
  return(p)
}

#' Plot feature expression in reduced dimensions
#'
#' This function plots feature expression in reduced dimensions using the FeaturePlot function from the Seurat package.
#'
#' @param object A Seurat object containing the data to plot.
#' @param features A character vector of feature names to plot.
#' @param reduction A character string specifying the dimensionality reduction to use (e.g., "PCA", "UMAP", "tSNE"). Defaults to NULL.
#' @param label A character vector specifying the labels to use for each cell group. Defaults to label.
#' @param split_by A character vector specifying the cell groups to split the plot by. Defaults to NULL.
#' @param label_size A numeric value specifying the size of the labels. Defaults to 8 * 0.36.
#' @param pt_size A numeric value specifying the size of the points. Defaults to 1.
#' @param slot A character string specifying which slot in the Seurat object to use (e.g., "data", "scale.data", "integrated"). Defaults to "data".
#' @param max_cutoff A numeric value specifying the maximum expression cutoff. Defaults to NA.
#' @param raster A logical value specifying whether to use raster graphics. Defaults to TRUE.
#' @param seed An integer value specifying the random seed. Defaults to 717.
#' @param title A character string specifying the plot title. Defaults to NULL.
#' @param legend_title A character string specifying the legend title. Defaults to NULL.
#' @param show_legend A logical value specifying whether to show the legend. Defaults to TRUE.
#' @param show_axis A logical value specifying whether to show the plot axis. Defaults to TRUE.
#' @param show_border A logical value specifying whether to show the plot border. Defaults to TRUE.
#' @param palette A character string specifying the color palette to use. Defaults to "YlOrRd".
#' @param direction A numeric value specifying the direction of the color palette. Defaults to 1.
#' @param ... Additional parameters to pass to FeaturePlot.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 labs scale_color_distiller guides
#' @importFrom ggplot2 theme element_blank element_text margin
#'
#' @examples
#' \dontrun{
#' sn_plot_feature(x = mySeuratObject, features = c("CD3D", "CD8A", "CD4"), reduction = "UMAP")
#' }
#'
#' @export
sn_plot_feature <-
  function(object,
           features,
           reduction = NULL,
           label = label,
           split_by = NULL,
           label_size = 8 * 0.36,
           pt_size = 1,
           slot = "data",
           max_cutoff = NA,
           raster = TRUE,
           seed = 717,
           title = NULL,
           legend_title = NULL,
           show_legend = TRUE,
           show_axis = TRUE,
           show_border = TRUE,
           palette = "YlOrRd",
           direction = 1,
           ...) {
    p <- Seurat::FeaturePlot(
      object = object,
      features = features,
      reduction = reduction,
      split.by = split_by,
      label.size = 8 * 0.36,
      pt.size = pt_size,
      slot = slot,
      raster = raster,
      order = TRUE,
      max.cutoff = max_cutoff,
      ...
    )
    if (!show_legend) {
      p <- p + Seurat::NoLegend()
    }
    reduction <- reduction %||% SeuratObject::DefaultDimReduc(object = object)
    object[["ident"]] <- Seurat::Idents(object = object)
    title <- title %||% features
    p <- p +
      catplot::theme_cat(aspect_ratio = 1) +
      theme(legend.margin = margin(l = -8)) +
      labs(
        x = paste0(stringr::str_to_upper(reduction), " 1"),
        y = paste0(stringr::str_to_upper(reduction), " 2"),
        title = title
      )

    if (!is_null(palette)) {
      p <-
        p + scale_color_distiller(palette = palette, direction = direction)
    }
    if (!is_null(legend_title)) {
      p <- p + guides(color = guide_colorbar(title = legend_title))
    }
    if (!show_axis) {
      p <- p + theme(
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    }
    if (!show_border) {
      p <- p + theme(
        panel.border = element_blank(),
        axis.line = element_blank()
      )
    }
    return(p)
  }


add_palette <- function(p, palette, n) {
  if (length(palette) > 1) {
    values <- palette
  } else if (length(palette) == 1 &&
    palette %in% names(palette_db)) {
    values <- palette_db[[palette]]
  } else if (length(palette) == 1 &&
    palette %in% row.names(RColorBrewer::brewer.pal.info)) {
    values <- brewer.pal(12, palette)
  } else {
    stop(
      "Palette not found, Please check the name of the palette use `show_all_palettes()` to see all available palettes"
    )
  }

  if (n > length(values)) {
    values <- colorRampPalette(values)(n)
  } else {
    values <- values[1:n]
  }
  p <- p + scale_color_manual(values = values)
}


palette_db <- vector("list")

palette_db$ZhangJian2024 <- c(
  "#efcec9", "#ff8c72", "#23676e",
  "#fd70a9", "#aa96c0", "#4194d0",
  "#83c066", "#ffba64", "#3fa177",
  "#a6846a", "#49548f", "#34405c",
  "#f73a41", "#3284b8", "#8dd3c9",
  "#3aa08e", "#726f83"
)

palette_db$XuPan2024 <- c(
  "#c4b1a7", "#867790", "#af8e87", "#ecb772",
  "#af8e87", "#8ca287", "#f1be94", "#bc966d", "#a5beba", "#de9590", "#a3b8c5"
)


#' @export
show_all_palettes <- function() {
  # Collection from paper
  cat("Collection from paper\n")
  print(names(palette_db))
  # RColorBrewer
  cat("RColorBrewer\n")
  print(row.names(RColorBrewer::brewer.pal.info))
}

.all_palettes <- c(names(palette_db), row.names(RColorBrewer::brewer.pal.info))

#
# c("#efcec9", "#ff8c72", "#fd70a9", "#ffba64", "#f73a41", "#a6846a",
#   "#23676e", "#4194d0", "#3284b8", "#49548f", "#34405c", "#8dd3c9", "#3aa08e", "#83c066",
#   "#aa96c0", "#726f83")
#
#
# colorspace::swatchplot(
#     x = RColorBrewer::brewer.pal(n = 12, "Paired")
# )

# p <- iris |>
#     ggplot(aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
#     geom_point()
# p
# bb <- ggplot_build(p)
# bb$data[[1]] |> head()
