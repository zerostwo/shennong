# Reuse the derived PBMC fixture from embedding-styles.R; leave it unchanged.
devtools::load_all(quiet=TRUE)
o <- qs2::qs_read("dev/umap/style-demo/pbmc-3d.qs2")
o <- Seurat::RunUMAP(o, reduction="pca", dims=1:20, n.components=2L, reduction.name="umap", reduction.key="UMAP_", seed.use=717, verbose=FALSE)
out <- "dev/umap/style-demo-2d"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
qs2::qs_save(o, file.path(out,"pbmc-2d.qs2"))
plots <- list()
for (style in c("classic","nebula","glass")) {
 p <- sn_plot_dim(o, reduction="umap", group_by="seurat_clusters", style=style, label=TRUE, show_legend=FALSE, raster_dpi=c(600,600), title=paste("UMAP 1 / UMAP 2 -",style))
 plots[[style]] <- p
 ggplot2::ggsave(file.path(out,paste0(style,".png")), p, width=6,height=6,dpi=160)
 ggplot2::ggsave(file.path(out,paste0(style,".pdf")), p, width=6,height=6,dpi=600)
 if(style!="classic") {
  w <- sn_plot_dim(o, reduction="umap", group_by="seurat_clusters", style=style, label=TRUE, interactive=TRUE)
  htmlwidgets::saveWidget(w,file.path(out,paste0(style,".html")),selfcontained=TRUE)
 }
}
ggplot2::ggsave(file.path(out,"comparison.png"), patchwork::wrap_plots(plots,nrow=1), width=18,height=6,dpi=120)
p <- sn_plot_feature(o,"NKG7",reduction="umap",group_by="seurat_clusters",style="glass",raster_dpi=600)
ggplot2::ggsave(file.path(out,"NKG7.pdf"),p,width=6,height=5,dpi=600)
ggplot2::ggsave(file.path(out,"NKG7.png"),p,width=6,height=5,dpi=150)