# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggthemes")
  library("ggrepel")
  library("grid")
  library("DoubletFinder")
  library("doMC")
  library("xlsx")
  library("RColorBrewer")
  library(scCustomize)
  library(paletteer)
  library(stringr)
})

# Function to print clear log progress updates
message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

# Filter operator
`%notin%` <- Negate(`%in%`)

# Create directories 
message2("Creating directories")

proj_dir <- "/projects/p31535/boles/img_scfrp_pilot/"

data_in_dir <- paste0(proj_dir, "data/00_soupx/")

plots_dir <- paste0(proj_dir, "plots/01_qc1/")
dir.create(plots_dir, showWarnings = F,
           recursive = T)

csv_dir <- paste0(proj_dir, "tab_data/01_qc1/")
dir.create(csv_dir, showWarnings = F,
           recursive = T)

data_out_dir <- paste0(proj_dir, "data/01_qc1/")
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

# Create Seurat objects from SoupX-corrected integer matrices 
message2("Creating Seurat objects")

create_object <- function(sample){
  counts <- Read10X(paste0(data_in_dir, sample))
  return(CreateSeuratObject(counts = counts, project = sample))
}

samples <- list.dirs(data_in_dir,
                     recursive = F,
                     full.names = F)

list <- list()

for (i in seq_along(samples)){
  message(paste0("Creating object for ", samples[i]))
  list[i] <- create_object(samples[i])
}

message2("Merging Seurat objects")
obj <- Merge_Seurat_List(list,
                         add.cell.ids = samples)

message2("Joining layers")
obj <- JoinLayers(obj)

# Add QC metrics & organize object for plotting 
obj <- obj %>%
  Add_Cell_QC_Metrics(species = "human")

message2("Organizing object for QC plots")

obj@meta.data <- obj@meta.data %>%
  mutate(orig.ident = factor(orig.ident,
                             levels = c("KOLF", "GALC", "GBA", "LRRK2")))
Idents(obj) <- factor(Idents(obj),
                      levels = c("KOLF", "GALC", "GBA", "LRRK2"))

pal <- paletteer_d("tvthemes::kimPossible")
obj <- Store_Palette_Seurat(seurat_object = obj, 
                            palette = pal, palette_name = "sample_pal",
                            overwrite = T)

message2("Drawing pre-filter QC plots")

p1 <- QC_Plots_Genes(obj,
                     colors_use = obj@misc$sample_pal,
                     plot_boxplot = T) + 
  ggtitle("Genes per cell")
p2 <- QC_Plots_UMIs(obj,
                    colors_use = obj@misc$sample_pal,
                    y_axis_log = T,
                    plot_boxplot = T) + 
  ylab("log10(UMIs)") + 
  ggtitle("UMIs per cell")
p3 <- QC_Plots_Mito(obj,
                    colors_use = obj@misc$sample_pal,
                    plot_boxplot = T) + 
  ggtitle("Mitochondrial gene content per cell") +
  ylab("% mitochondrial\ngene counts")
p4 <- QC_Plots_Complexity(obj,
                          colors_use = obj@misc$sample_pal,
                          raster = F,
                          plot_boxplot = T) + 
  ggtitle("Cell complexity")

ggsave(plot = p1,
       filename = paste0(plots_dir, "nfeature.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)
ggsave(plot = p2,
       filename = paste0(plots_dir, "ncount.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)
ggsave(plot = p3, 
       filename = paste0(plots_dir, "mito.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)
ggsave(plot = p4,
       filename = paste0(plots_dir, "complexity.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

message2("Setting thresholds and drawing new plots")

mito_max <- 5
complex_min <- 0.8

feature_med <- median(obj@meta.data$nFeature_RNA)
feature_min <- min(obj@meta.data$nFeature_RNA)
feature_max <- max(obj@meta.data$nFeature_RNA)
feature_lb <- feature_med - 2*stats::mad(obj@meta.data$nFeature_RNA)
feature_ub <- feature_med + 5*stats::mad(obj@meta.data$nFeature_RNA)

umi_med <- median(obj@meta.data$nCount_RNA)
umi_lb <- umi_med - 3*stats::mad(obj@meta.data$nCount_RNA)
# this is a negative value, set to firm min
umi_lb <- 5000
umi_ub <- umi_med + 5*stats::mad(obj@meta.data$nCount_RNA)

cutoff <- function(intercept){
  geom_hline(yintercept = intercept,
             linetype = "longdash",
             color = "gray60",
             linewidth = 2)
}

ggsave(p1 + 
         cutoff(feature_lb) + 
         cutoff(feature_ub),
       filename = paste0(plots_dir, "nfeature_cutoffs.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

ggsave(p2 +
         cutoff(umi_lb) +
         cutoff(umi_ub),
       filename = paste0(plots_dir, "ncount_cutoffs.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

ggsave(p3 +
         cutoff(mito_max),
       filename = paste0(plots_dir, "mito_cutoffs.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

ggsave(p4 +
         cutoff(complex_min),
       filename = paste0(plots_dir, "complexity_cutoffs.png"),
       units = "in", dpi = 600,
       height = 3, width = 6)

message2("Computing stats")

median_stats <- Median_Stats(seurat_object = obj, 
                             group_by_var = "orig.ident")
median_stats

counts <- obj@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::summarize(Cell_count = n())
  
median_stats <- median_stats %>%
  left_join(counts, by = "orig.ident") 

median_stats[5, 6] <- sum(counts$Cell_count)

write.csv(median_stats,
          file = paste0(csv_dir, "median_stats_prefilter.csv"),
          row.names = F)

message2("Filtering object")

obj_filtered <- subset(x = obj,
                       subset = nCount_RNA > umi_lb &
                         nCount_RNA < umi_ub &
                         nFeature_RNA > feature_lb & 
                         nFeature_RNA < feature_ub &
                         percent_mito < mito_max & 
                         log10GenesPerUMI > complex_min)

# message2("Clustering to find any remaining low quality groups of cells")
# obj_filtered <- obj_filtered %>%
#   SCTransform(verbose = T) %>%
#   RunPCA() 
# 
# p <- ElbowPlot(obj_filtered)
# ggsave(p,
#        filename = paste0(plots_dir, "pca_elbow.png"),
#        units = "in", dpi = 600,
#        height = 3, width = 3)
# 
# obj_filtered <- obj_filtered %>%
#   RunUMAP(dims = 1:10) %>%
#   FindNeighbors(dims = 1:10) 
# 
# obj_filtered <- obj_filtered %>%
#   FindClusters(resolution = 0.4)
# 
# QC_Plots_Genes(obj_filtered,
#                     plot_boxplot = T) + 
#   ggtitle("Genes per cell")
# QC_Plots_UMIs(obj_filtered,
#                     y_axis_log = T,
#                     plot_boxplot = T) + 
#   ylab("log10(UMIs)") + 
#   ggtitle("UMIs per cell")
# QC_Plots_Mito(obj_filtered,
#                     plot_boxplot = T) + 
#   ggtitle("Mitochondrial gene content per cell") +
#   ylab("% mitochondrial\ngene counts")
# QC_Plots_Complexity(obj_filtered,
#                          raster = F,
#                           plot_boxplot = T) + 
#   ggtitle("Cell complexity")

median_stats_filtered <- Median_Stats(obj_filtered,
                                      group_by = "orig.ident")

counts_filtered <- obj_filtered@meta.data %>%
  group_by(orig.ident) %>%
  dplyr::summarize(Cell_count = n())

median_stats_filtered <- median_stats_filtered %>%
  left_join(counts_filtered, by = "orig.ident")

median_stats_filtered[5, 6] <- sum(counts_filtered$Cell_count)

write.csv(median_stats_filtered,
          file = paste0(csv_dir, "median_stats_postfilter.csv"))

message2("Creating before and after median stats figure")

colnames(median_stats) <- c("id", "nUMI", "nFeature", 
                            "percent_mito", "complexity", "cell_count")

colnames(median_stats_filtered) <- colnames(median_stats)

colnames(median_stats_filtered)[-1] <- paste0(colnames(median_stats_filtered)[-1],
                                              "_filtered")

stats <- median_stats %>%
  left_join(median_stats_filtered, by = "id")

stats_p <- stats %>%
  pivot_longer(cols = !id,
               names_to = "stat") %>%
  mutate(filter = if_else(str_detect(stat, "filtered"), "post", "pre")) %>%
  mutate(stat = str_remove_all(stat, "_filtered") %>%
           factor(levels = c("cell_count", "nFeature", "nUMI", "complexity", "percent_mito"),
                  labels = c("Cell count", "Median number\nof genes",
                             "Median number\nof UMIs", "Median log10(genes) /\n log10(UMIs)",
                             "Median %\nmitochondrial genes"))) %>%
  # dplyr::filter(stat == "cell_count") %>%
  ggplot(aes(x = factor(filter,
                        levels = c("pre", "post"),
                        labels = c("Pre-filtering", "Post-filtering")), y = value)) + 
  geom_point(aes(fill = id, shape = id),
             size = 5) +
  geom_line(aes(color = id, group = id),
            show.legend = F) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  scale_shape_manual(values = c(21:25)) +
  facet_wrap(. ~ stat,
             scales = "free",
             axes = "all_y",
             ncol = 5) + 
  theme_bw(base_size = 10) +
  theme(
    axis.title = element_blank(),
    legend.title = element_blank(),
    axis.text = element_text(color = "black"),
    strip.text = element_text(color = "white", face = "bold"),
    strip.background = element_rect(fill = "black")
  )

ggsave(stats_p,
       filename = paste0(plots_dir, "qc_change.png"),
       units = "in", dpi = 600,
       height = 3, width = 12)

message2("Saving filtered object to data/01_qc1/")

saveRDS(obj_filtered,
        file = paste0(data_out_dir, "01_filtered.rds"))