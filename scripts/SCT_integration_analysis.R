# Single-cell RNA-seq - normalization

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)


# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
seurat_phase

view(seurat_phase)


# Load cell cycle markers
load("data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

view(seurat_phase@meta.data)

# View cell cycle scores and phases assigned to cells
head(seurat_phase@meta.data)                                


###SO WE ARE GOING IN WITH THE TOP 2000 GENES WITH THE MOST VARIABLILITY IN THE SAMPLES

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = T)



# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")



# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)



# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))


view(seurat_phase@meta.data)

#######

DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "mitoFr",
        split.by = "mitoFr")






########################################################################

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- readRDS("data/split_seurat.rds")
ctrl_reps <- split_seurat[c("ctrl_1", "ctrl_2")]
options(future.globals.maxSize = 4000 * 1024^2)


for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                   vars.to.regress = c("mitoRatio"),
                                   vst.flavor = "v2")
}


# Check which assays are stored in objects

split_seurat$ctrl@assays

split_seurat$stim@assays

split_seurat$stim@assays[["SCT"]]
split_seurat$ctrl@assays[["SCT"]]
########################################################################
# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")
########################################################################
# Load the split seurat object into the environment
split_seurat <- readRDS("data/split_seurat.rds")
