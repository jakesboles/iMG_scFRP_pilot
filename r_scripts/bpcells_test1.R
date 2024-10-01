library(scCustomize)
library(Seurat)
library(BPCells)

data <- "/projects/b1169/boles/als_multitissue_scfrp/data/00_soupx/AU-066_b/"
test <- "/projects/p31535/boles/als_multitissue_scfrp/test/"
dir.create(test)

counts <- Read10X(data.dir = data)
obj <- CreateSeuratObject(counts, project = "test")

write_matrix_dir(mat = obj[["RNA"]]$counts, dir = paste0(test, "bpcells_au066b"))


