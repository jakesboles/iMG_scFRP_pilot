library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)

proj_dir <- "/projects/p31535/boles/img_scfrp_pilot/"

in_dir <- paste0(proj_dir, "data/04_de/")

ranked_degs <- readRDS(paste0(in_dir, "ranked_degs.rds"))

bg_genes <- read.csv("/projects/p31535/Anne/cellranger_reffiles/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv",
                     skip = 5)
bg_genes <- bg_genes %>%
  dplyr::pull(gene_id)
  
bg_genes <- bitr(bg_genes, 
     fromType = "ENSEMBL",
       toType = 'SYMBOL',
       OrgDb = org.Hs.eg.db) %>%
  dplyr::pull(SYMBOL) %>%
  unique()
# Break "ties" at infinite and zero values by simply extending out the mins/maxes 
# p-values from Wilcox testing were sometimes == 0, meaning -log10(p) = Inf
# This is the product of an inappropriate analysis but there are few other 
# options with n = 1
for (i in seq_along(ranked_degs)){
  # Extend out non-Inf max to Inf values
  max <- max(ranked_degs[[i]][ranked_degs[[i]] != Inf])
  to_change_top <- ranked_degs[[i]][ranked_degs[[i]] == Inf]
  
  l <- length(to_change_top)
  
  for (j in seq_along(to_change_top)){
    to_change_top[j] <- max + (2*l - j)
  }
  
  # Extend out non- -Inf min to -Inf values 
  min <- min(ranked_degs[[i]][ranked_degs[[i]] != -Inf])
  to_change_bottom <- ranked_degs[[i]][ranked_degs[[i]] == -Inf]
  
  l <- length(to_change_bottom)
  
  for (j in seq_along(to_change_bottom)){
    to_change_bottom[j] <- min - (2*l + j)
  }
  
  # Break ties at zero 
  zero <- ranked_degs[[i]][ranked_degs[[i]] == 0]
  pos_min <- min(ranked_degs[[i]][ranked_degs[[i]] > 0])
  
  l <- length(zero)
  
  step <- pos_min / l
  
  for (j in seq_along(zero)){
    zero[j] <- pos_min - (step * j)
  }
  
  ranked_degs[[i]][ranked_degs[[i]] == Inf] <- to_change_top
  ranked_degs[[i]][ranked_degs[[i]] == -Inf] <- to_change_bottom
  ranked_degs[[i]][ranked_degs[[i]] == 0] <- zero
}

names(ranked_degs) <- c("GALC", "GBA", "LRRK2")

gsea <- compareCluster(ranked_degs,
                       fun = gseGO,
                       ont = "BP",
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       pvalueCutoff = 1)

tabs <- readRDS(paste0(in_dir, "wilcox_de_tables.rds"))

degs <- list()
for (i in seq_along(tabs)){
  degs[[i]] <- tabs[[i]] %>%
    dplyr::filter(abs(avg_log2FC) > 0.5 &
                    p_val_adj < 0.05) %>%
    rownames()
}
names(degs) <- c("GALC", "GBA", "LRRK2")

enrich <- compareCluster(degs,
                       fun = enrichGO,
                       ont = "BP",
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       pvalueCutoff = 1,
                       universe = bg_genes,
                       readable = T)

df <- enrich@compareClusterResult %>%
  as.data.frame()

df %>%
  filter(Cluster == "GALC" &
           p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>% 
  ggplot(aes(x = -log10(p.adjust),
             y = fct_inorder(Description))) +
  geom_col(aes(fill = Count)) +
  geom_vline(xintercept = -log10(0.05)) +
  labs(x = "-log10(BH-adjusted p-value)") +
  scale_y_discrete(limits = rev)
