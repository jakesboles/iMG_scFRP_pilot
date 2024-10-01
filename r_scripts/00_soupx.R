suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("Matrix")
  library("SoupX")
  library("DropletUtils")
  library("purrr")
  library("Hmisc")
})

# Code copied from SoupX SEA-AD script from Anne

run_soupx <- function(dir) {
  sample <- unlist(strsplit(dir, "/")) %>%
    tail(1)
  
  sample_dir <- paste0(data_out_dir, sample, "/")
  
  message(paste0("Processing sample: ", sample))
  
  # Initialize Seurat object
  seurat <- Read10X(paste0(dir, "/count/sample_filtered_feature_bc_matrix")) %>%
    CreateSeuratObject()
  
  # Run SCTransform and PCA 
  seurat <- SCTransform(seurat) %>% RunPCA()
  
  # Generate clusters using default parameters (SoupX is not sensitive to clustering quality - Seurat default should be fine)
  seurat <- seurat %>% FindNeighbors(dims = 1:10) %>% FindClusters() 
  
  # Load in 10X data for SoupX 
  toc = Seurat::Read10X(paste0(dir,"/count/sample_filtered_feature_bc_matrix"))
  tod = Seurat::Read10X(paste0(dir,"/count/sample_raw_feature_bc_matrix"))
  
  # Filter table of droplets (tod) so that it has same number of genes as table of counts (toc)
  tod <- tod[which(rownames(tod) %in% rownames(toc)),]
  
  # Create SoupX object (default calcSoupProfile = TRUE)
  soupx = SoupChannel(tod, toc)
  
  # Add cluster info to SoupX object
  soupx <- setClusters(soupx, setNames(seurat$seurat_clusters, colnames(seurat)))
  
  # Save automated estimate plot of contamination fraction 
  pdf(file = paste0(plots_dir, sample, "_rho_distribution.pdf"), width = 4, height = 4)
  soupx <- autoEstCont(soupx, forceAccept = TRUE) # Set forceAccept = TRUE to allow to proceed with high contamination fraction, can look into problematic samples later and decide to exclude or manually set the fraction more in line with other samples
  dev.off()
  
  # Adjust counts, rounded to integer
  message(paste0(sample, " Contamination Fraction: ", soupx$metaData$rho[1])) # Print contamination fraction
  rounded <- adjustCounts(soupx, roundToInt = TRUE)
  
  # Save corrected counts
  DropletUtils:::write10xCounts(sample_dir, rounded, 
                                overwrite = TRUE, version = "3")
  
  return(soupx$metaData$rho[1])
}

# Define directories

proj_dir <- "/projects/p31535/boles/img_scfrp_pilot/"

data_in_dir <- "/projects/b1042/Gate_Lab/boles/img_scfrp_pilot/cellranger/"

plots_dir <- paste0(proj_dir, "plots/00_soupx/")
dir.create(plots_dir, showWarnings = T,
           recursive = T)

csv_dir <- paste0(proj_dir, "tab_data/00_soupx/")
dir.create(csv_dir, showWarnings = T,
           recursive = T)

data_out_dir <- paste0(proj_dir, "data/00_soupx/")
dir.create(data_out_dir, showWarnings = T,
           recursive = T)

samples <- list.dirs(paste0(data_in_dir, "outs/per_sample_outs/"),
                     recursive = F, full.names = F)

cellranger_dir <- paste0(data_in_dir, "outs/per_sample_outs")

# Create general output folder
for (i in seq_along(samples)){
  output_dir <- paste0(data_out_dir, samples[i], "/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# Directories for sample-level cellranger output 
sample_dirs <- list.dirs(cellranger_dir, recursive = FALSE)
print(sample_dirs)

# Run SoupX over each Cell Ranger directory
rho <- lapply(sample_dirs, run_soupx)

names(rho) <- samples
rho <- unlist(rho)
print(rho)

# Export rho (contamination fractions per sample) as a csv
write.csv(rho, file = paste0(csv_dir, "/contamination_fraction.csv"))

