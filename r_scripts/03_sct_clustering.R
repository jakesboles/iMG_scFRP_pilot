library(Seurat)
library(scCustomize)
library(tidyverse)
library(stringr)
library(glmGamPoi)
library(sctransform)
library(ggpubr)
library(stringr)
library(paletteer)
library(dittoSeq)
library(patchwork)

# Define directories and load in objects

proj_dir <- '/projects/p31535/boles/img_scfrp_pilot/'

files <- list.files(paste0(proj_dir, "data/02_doubletfinder/"))

obj_list <- list()

for (i in seq_along(files)){
  obj_list[[i]] <- readRDS(paste0(proj_dir, "data/02_doubletfinder/", files[i]))
  names(obj_list)[i] <- str_remove_all(files[i], ".rds")
}

# Collapse list into one Seurat object 
obj <- Merge_Seurat_List(obj_list,
                         add.cell.ids = names(obj_list)) %>%
  JoinLayers()

# Make some metadata columns into factors (just for ordering in plots)
obj@meta.data <- obj@meta.data %>%
  mutate(orig.ident = factor(orig.ident,
                             levels = c("KOLF", "GALC", "GBA", "LRRK2")),
         DF.adj = factor(DF.adj,
                         levels = c("Singlet", "Doublet")),
         DF.unadj = factor(DF.unadj,
                           levels = c("Singlet", "Doublet")))
Idents(obj) <- "orig.ident"

# Plotting tools
pal <- paletteer_d("tvthemes::kimPossible")
obj <- Store_Palette_Seurat(seurat_object = obj, 
                            palette = pal, palette_name = "sample_pal",
                            overwrite = T)

umap_labs <- list(
  ylab("UMAP 2"),
  xlab("UMAP 1")
)
umap_theme <- theme(
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
                                                 ends = "last"))
)

plots_dir <- paste0(proj_dir, "plots/03_sct_clustering/")
dir.create(plots_dir, 
           showWarnings = F, recursive = T)

# SCT normalization
obj <- SCTransform(obj, verbose = T)

# Run linear dimensionality reduction and determine number of PCs to carry forward
obj <- obj %>%
  RunPCA()

p <- ElbowPlot(obj)
# Inspecting this plot shows an inflection at around 5, but we'll keep 10 for now
ggsave(p,
       filename = paste0(plots_dir, "pca_elbow.png"),
       units = "in", dpi = 600,
       height = 3, width = 3)

pca_plot <- function(dims){
  DimPlot_scCustom(obj,
                   reduction = "pca",
                   label = F,
                   dims = dims,
                   colors_use = obj@misc$sample_pal) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10))
}

# See if samples distribute differently across PCs 
p <- ggarrange(
  pca_plot(c(1,2)), NULL, NULL,
  pca_plot(c(1,3)), pca_plot(c(2,3)), NULL,
  pca_plot(c(1,4)), pca_plot(c(2,4)), pca_plot(c(3,4)),
  ncol = 3, nrow = 3,
  common.legend = T, legend = "right",
  align = "hv")
ggsave(p,
       filename = paste0(plots_dir, "pca_grid.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

# Further dimensionality reduction with some inspection
obj <- obj %>%
  RunUMAP(dims = 1:10) %>%
  FindNeighbors(dims = 1:10)

# Visualize doublets
doublet_pal <- c("springgreen", "navy")

p1 <- DimPlot_scCustom(obj,
                       group.by = "DF.unadj",
                       label = F,
                       colors_use = doublet_pal) + 
  ggtitle("Unadjusted doublate rate") +
  umap_labs + 
  umap_theme

p2 <- DimPlot_scCustom(obj,
                       group.by = "DF.adj",
                       label = F,
                       colors_use = doublet_pal) + 
  ggtitle("Adjusted for homotypic proportion") +
  umap_labs + 
  umap_theme

p <- p1 + p2 +
  plot_layout(guides = "collect")

ggsave(p,
       filename = paste0(plots_dir, "umap_doublets.png"),
       units = "in", dpi = 600,
       height = 4, width = 9)
# Seems like doublets distribute unevenly across the UMAP space

p <- ggarrange(
  FeaturePlot_scCustom(obj,
                       features = "nCount_RNA",
                       na_cutoff = NULL) +
    umap_labs + 
    umap_theme,
  QC_Plots_UMIs(obj,
                group.by = "DF.adj",
                plot_boxplot = T),
  FeaturePlot_scCustom(obj,
                       features = "nFeature_RNA",
                       na_cutoff = NULL) +
    umap_labs + 
    umap_theme,
  QC_Plots_Genes(obj,
                 group.by = "DF.adj",
                 plot_boxplot = T),
  ncol = 2, nrow = 2, widths = c(1.5, 1, 1.5, 1),
  align = "v")
ggsave(p,
       filename = paste0(plots_dir, "doublet_qc.png"),
       units = "in", dpi = 600,
       height = 8, width = 8)

# Remove doublets predicted with homotypic proportion accounted for,
# re-SCT, then do some preliminary clustering

obj <- obj %>%
  subset(subset = DF.adj == "Singlet") %>%
  SCTransform() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.4)

p <- DimPlot_scCustom(obj,
                      colors_use = obj@misc$sample_pal,
                      group.by = "orig.ident",
                      label = F) +
  umap_labs + 
  umap_theme
ggsave(p,
       filename = paste0(plots_dir, "umap_sample.png"),
       units= "in", dpi = 600,
       height = 4, width = 5)

# Preliminary clustering & marker gene identification

p <- DimPlot_scCustom(obj,
                      group.by = "SCT_snn_res.0.4") + 
  umap_theme + 
  umap_labs
ggsave(p,
       filename = paste0(plots_dir, "umap_clusters.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

markers <- FindAllMarkers(obj)
top_markers <- Extract_Top_Markers(markers, num_genes = 5,
                                   make_unique = T,
                                   named_vector = F)

csv_dir <- paste0(proj_dir, "tab_data/03_sct_clustering/")
dir.create(csv_dir, showWarnings = F, recursive = T)
write.csv(markers,
          file = paste0(csv_dir, "cluster_markers.csv"),
          row.names = F)

# Heatmap of marker genes for each cluster 
pal <- DiscretePalette_scCustomize(num_colors = 11, palette = "polychrome")

png(filename = paste0(plots_dir, "clustering_0-4_marker_heatmap.png"),
    units = "in", res = 600,
    height = 8, width = 8)
dittoHeatmap(obj,
             genes = top_markers,
             scaled.to.max = T,
             cluster_cols = F,
             cluster_rows = F,
             order.by = "SCT_snn_res.0.4",
             annot.by = "SCT_snn_res.0.4",
             annot.colors = pal,
             treeheight_row = 0,
             heatmap.colors.max.scaled = paletteer_c("grDevices::Oslo", 30,
                                                     direction = 1))
dev.off()

# Looking for sample-specific clusters
obj@meta.data <- obj@meta.data %>%
  mutate(SCT_snn_res.0.4 = factor(SCT_snn_res.0.4,
                                  levels = c(0:10)))

p <- dittoBarPlot(obj,
             group.by = "SCT_snn_res.0.4",
             var = "orig.ident",
             color.panel = obj@misc$sample_pal,
             retain.factor.levels = T,
             main = "Fraction of cluster from each sample",
             ylab = "Proportion of cluster",
             xlab = "Cluster",
             x.labels.rotate = F) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )
ggsave(p,
       filename = paste0(plots_dir, "clusters_by_sample.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

# Cluster 3, while evenly distributed across all samples, shows an enrichment of 
# mitochondrial genes. A lot of these genes contributed to some of the first PCs,
# so show clusters by PC and visualize PC loadings 
pca_plot2 <- function(dims){
  DimPlot_scCustom(obj,
                   reduction = "pca",
                   label = F,
                   dims = dims) +
    # guides(color = guide_legend(position = "inside")) +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.justification = "left")
}

p <- ggarrange(
  pca_plot2(c(1,2)), NULL, NULL,
  pca_plot2(c(1,3)), pca_plot2(c(2,3)), NULL,
  pca_plot2(c(1,4)), pca_plot2(c(2,4)), pca_plot2(c(3,4)),
  ncol = 3, nrow = 3,
  common.legend = T, legend = "right",
  align = "hv")
ggsave(p,
       filename = paste0(plots_dir, "pca_grid_by_cluster.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

Iterate_PC_Loading_Plots(obj,
                         file_path = plots_dir,
                         file_name = "pca_loadings.pdf")

# Lastly, a few Feature Plots to highlight microglial genes 
feature_plot <- function(gene){
  p1 <- FeaturePlot_scCustom(obj,
                             features = gene,
                             reduction = "umap",
                             colors_use = viridis_inferno_light_high) + 
    theme(plot.title = element_text(face = "bold.italic")) +
    umap_theme + 
    umap_labs
  
  p2 <- VlnPlot_scCustom(obj,
                         features = gene,
                         group.by = "SCT_snn_res.0.4",
                         colors_use = pal,
                         slot = "counts") +
    xlab("Cluster") +
    ylab("Normalized counts") +
    theme(plot.title = element_blank(),
          legend.position = "none")
  
  p3 <- VlnPlot_scCustom(obj,
                         features = gene,
                         group.by = "orig.ident",
                         colors_use = obj@misc$sample_pal,
                         slot = "counts") +
    xlab("Sample") +
    ylab("Normalized counts") +
    theme(plot.title = element_blank(),
          legend.position = "none")
  
  design <- "
  AAAAAAAA
  BBBBBCCC
  "
  
  p1 + p2 + p3 +
    plot_layout(design = design)
}

p <- feature_plot("ITGAM")
ggsave(p,
       filename = paste0(plots_dir, "microglia_gene_itgam.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

p <- feature_plot("AIF1")
ggsave(p,
       filename = paste0(plots_dir, "microglia_gene_aif1.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

p <- feature_plot("CSF1R")
ggsave(p,
       filename = paste0(plots_dir, "microglia_gene_csf1r.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

p <- feature_plot("SPI1")
ggsave(p,
       filename = paste0(plots_dir, "microglia_gene_spi1.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

# Any differences in the genes that are supposed to be knocked out? 
p <- feature_plot("LRRK2")
ggsave(p,
       filename = paste0(plots_dir, "feature_lrrk2.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

p <- feature_plot("GBA1")
ggsave(p,
       filename = paste0(plots_dir, "feature_gba1.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

p <- feature_plot("GALC")
ggsave(p, filename = paste0(plots_dir, "feature_galc.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

# Save object for downstream DE, etc.
out_dir <- paste0(proj_dir, "data/03_sct_clustering/")
dir.create(out_dir,
           showWarnings = F,
           recursive = T)

saveRDS(obj,
        file = paste0(out_dir, "03_obj.rds"))
