library(Seurat)
library(scCustomize)
library(ggplot2)
library(paletteer)
library(tidyverse)
library(DESeq2)
library(UpSetR)
library(ggrepel)

# Load project directory and data object ---------------------------------

proj_dir <- "/projects/p31535/boles/img_scfrp_pilot/"

obj <- readRDS(paste0(proj_dir, "data/03_sct_clustering/03_obj.rds"))

data_dir <- paste0(proj_dir, "data/04_de/")
dir.create(data_dir, 
           showWarnings = F, 
           recursive = T)

csv_dir <- paste0(proj_dir, "tab_data/04_de/")
dir.create(csv_dir,
           showWarnings = F,
           recursive = T)
# Plotting tools ----------------------------------------------------------

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

plots_dir <- paste0(proj_dir, "plots/04_de/")
dir.create(plots_dir, 
           showWarnings = F, recursive = T)

volcano_plot <- function(de_table, num_top, num_bottom){
  d <- de_table %>%
    mutate(de = case_when(p_val_adj < 0.005 &
                            avg_log2FC > 0.5 ~ "Upregulated",
                          p_val_adj < 0.005 & 
                            avg_log2FC < -0.5 ~ "Downregulated",
                          p_val_adj > 0.005 ~ NA)) %>%
    rownames_to_column(var = "gene") %>%
    mutate(rank = -log10(p_val_adj) * avg_log2FC)
  
  top_labs <- d %>%
    filter(de == "Upregulated") %>%
    dplyr::arrange(desc(rank)) %>%
    head(n = num_top)
  
  bottom_labs <- d %>%
    filter(de == "Downregulated") %>%
    dplyr::arrange(rank) %>%
    head(n = num_bottom)
  
  labs <- rbind(top_labs, bottom_labs) %>%
    mutate(p_val_adj = if_else(rank == Inf, 1e-280, if_else(rank == -Inf, 1e-280, p_val_adj)))
  
  d %>%
    ggplot(aes(x = avg_log2FC,
               y = -log10(p_val_adj))) +
    geom_point(aes(fill = de),
               shape = 21,
               alpha = 0.5,
               size = 1.5) + 
    geom_hline(yintercept = -log10(0.005),
               linetype = "dotted") + 
    geom_vline(xintercept = c(0.5, -0.5),
               linetype = "dotted") + 
    geom_label_repel(data= labs,
                    aes(x = avg_log2FC, 
                        y = -log10(p_val_adj),
                        label = gene,
                        color = de),
                    show.legend = F) +
    scale_fill_manual(breaks = c("Downregulated", "Upregulated"),
                      values = c("#00F100", "#F100F1")) +
    scale_color_manual(breaks = c("Downregulated", "Upregulated"),
                       values = c("darkgreen", "darkmagenta")) +
    labs(x = "log2(fold change)",
         y = "-log10(adjusted p-value)") +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    theme_bw(base_size = 12) +
    theme(axis.text = element_text(color = "black"),
          legend.title = element_blank())
}

# Perform DE with Wilcox test for each sample vs KOLF ---------------------

Idents(obj) <- "orig.ident"

DimPlot_scCustom(obj,
                 colors_use = obj@misc$sample_pal,
                 label = F) +
  umap_theme + 
  umap_labs

galc_de <- FindMarkers(obj,
                       ident.1 = "GALC",
                       ident.2 = "KOLF",
                       slot = "counts")
head(galc_de)
write.csv(galc_de,
          file = paste0(csv_dir, "galc_vs_kolf_wilcox.csv"),
          row.names = F)
p <- volcano_plot(galc_de, 30, 10)
ggsave(p,
       filename = paste0(plots_dir, "galc_volcano.png"),
       units = "in", dpi = 600,
       height = 6, width = 10)


gba_de <- FindMarkers(obj,
                      ident.1 = "GBA",
                      ident.2 = "KOLF")
head(gba_de)
write.csv(gba_de,
          file = paste0(csv_dir, "gba_vs_kolf_wilcox.csv"),
          row.names = F)

p <- volcano_plot(gba_de, 15, 10)
ggsave(p,
       filename = paste0(plots_dir, "gba_volcano.png"),
       units = "in", dpi = 600,
       height = 6, width = 10)

lrrk2_de <- FindMarkers(obj,
                        ident.1 = "LRRK2",
                        indent.2 = "KOLF")
head(lrrk2_de)
write.csv(lrrk2_de,
          file = paste0(csv_dir, "lrrk2_vs_kolf_wilcox.csv"),
          row.names = F)

p <- volcano_plot(lrrk2_de, 10, 15)
ggsave(p,
       filename = paste0(plots_dir, "lrrk2_volcano.png"),
       units = "in", dpi = 600,
       height = 6, width = 10)

# Upset plot --------------------------------------------------------------

de_tables <- list(galc_de, gba_de, lrrk2_de)

de_genes <- list()

for (i in seq_along(de_tables)){
  de_genes[[i]] <- de_tables[[i]] %>%
    filter(abs(avg_log2FC) > 0.5 &
             p_val_adj < 0.005) %>%
    rownames()
}

names(de_genes) <- c("GALC", "GBA", "LRRK2")

png(filename = paste0(plots_dir, "upset.png"),
    units = "in", res = 600,
    width = 8, height = 5)
UpSetR::upset(fromList(de_genes),
              sets.x.label = "# DEGs",
              order.by = "freq",
              queries = list(
                list(query = intersects,
                     params = list("GALC"),
                     color = obj@misc$sample_pal[2],
                     active = T),
                list(query = intersects,
                     params = list("GBA"),
                     color = obj@misc$sample_pal[3],
                     active = T),
                list(query = intersects,
                     params = list("LRRK2"),
                     color = obj@misc$sample_pal[4],
                     active = T)
              ),
              sets.bar.color = obj@misc$sample_pal[2:4])
dev.off()

# Rank genes for GSEA -----------------------------------------------------

ranked_de_genes <- list()

for (i in seq_along(de_tables)){
  d <- de_tables[[i]] %>%
    mutate(rank = -log10(p_val_adj) * avg_log2FC) %>%
    dplyr::arrange(dplyr::desc(rank)) %>%
    rownames_to_column(var = "gene") %>%
    dplyr::select(c("gene", "rank"))
  
  ranked_de_genes[[i]] <- d$rank
  names(ranked_de_genes[[i]]) <- d$gene
}


# Save data ---------------------------------------------------------------

saveRDS(ranked_de_genes,
        file = paste0(data_dir, "ranked_degs.rds"))
saveRDS(de_tables,
        file = paste0(data_dir, "wilcox_de_tables.rds"))
