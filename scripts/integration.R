# Integration

# Delete all variables to save RAM space
rm (list = ls())

# Restart R session to have a clean start
.rs.restartR()

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load split seurat that contains SCT normalization for both samples
split_seurat <- readRDS("data/split_seurat.rds")

# Remove scale.data layer
# Only needed this for PCA which we calculated
split_seurat$ctrl[['RNA']]$scale.data <- NULL
split_seurat$stim[['RNA']]$scale.data <- NULL


# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

head(integ_features)
length(integ_features)

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Allow for more memory
options(future.globals.maxSize = 4000 * 1024^2)

# Clear unused memory
gc() 

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Clear unused memory
gc() 

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
# Clear unused memory
gc() 

# Removing split_seurat object from the environment to save RAM 
rm(split_seurat)

# Rejoin the layers in the RNA assay that we split earlier
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])

seurat_integrated

seurat_integrated@assays

# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Clear unused memory
gc() 

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "sample",
        group.by = "sample")  


# Set seed
set.seed(123456)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
DimPlot(seurat_integrated, group.by = "sample")                             

# Clear unused memory
gc() 

# Save integrated seurat object
saveRDS(seurat_integrated, "results/integrated_seurat.rds")

rm(integ_anchors)


seurat_integrated <- readRDS("results/integrated_seurat.rds")



