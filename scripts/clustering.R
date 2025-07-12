# Single-cell RNA-seq - clustering

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

seurat_integrated <- readRDS("results/integrated_seurat.rds")


# Explore heatmap of PCs
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)


# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)


# Plot the elbow plot
ElbowPlot(object = seurat_integrated, 
          ndims = 40)



# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)


# Determine the clusters for various resolutions                              
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))


seurat_integrated@meta.data %>% view()


Idents(object = seurat_integrated) %>% head()


# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"




# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)



# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


# Download the integrated object
download.file("https://d2lye0i2kmx8r6.cloudfront.net/cell-to-insight/seurat_integrated.RData.bz2", destfile="data/seurat_integrated.RData.bz2",headers =  c("Referer"="https://posit.cloud"))

# Load the seurat object
load(bzfile("data/seurat_integrated.RData.bz2"))


# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample)

# Barplot of number of cells per cluster by sample
ggplot(n_cells, aes(x=ident, y=n, fill=sample)) +
  geom_bar(position=position_dodge(), stat="identity") +
  theme_classic() +
  geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))



# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()



# Barplot of proportion of cells in each cluster by sample
ggplot(seurat_integrated@meta.data) +
  geom_bar(aes(x=integrated_snn_res.0.8, fill=sample), 
           position=position_fill())  +
  theme_classic()



# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()


# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)




# Boxplot of nGene per cluster
ggplot(seurat_integrated@meta.data) +
  geom_boxplot(aes(x=integrated_snn_res.0.8, y=nGene, 
                   fill=integrated_snn_res.0.8)) +
  theme_classic() +
  NoLegend()



# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)
head(pc_data)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  dplyr::summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)




# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)



DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()


DimPlot(object = seurat_integrated, 
        reduction = "umap", 
        label = TRUE) + NoLegend()



# Select the RNA counts slot to be the default assay
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
seurat_integrated




FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)



FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCGR3A", "MS4A7"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)



FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("MARCO", "ITGAM", "ADGRE1"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)




FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("FCER1A", "CST3"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)




FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


##lets make violin 

VlnPlot(seurat_integrated, 
             
            features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"))




# List of known celltype markers
markers <- list()
markers[["CD14+ monocytes"]] <- c("CD14", "LYZ")
markers[["FCGR3A+ monocyte"]] <- c("FCGR3A", "MS4A7")
markers[["Macrophages"]] <- c("MARCO", "ITGAM", "ADGRE1")
markers[["Conventional dendritic"]] <- c("FCER1A", "CST3")
markers[["Plasmacytoid dendritic"]] <- c("IL3RA", "GZMB", "SERPINF1", "ITM2C")

# Create dotplot based on RNA expression
DotPlot(seurat_integrated, markers, assay="RNA")



#homework
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("GNLY", "NKG7"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)


##nkcells, clusters 8,12 and maybe 5



DefaultAssay(seurat_integrated) <- "RNA"

Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"

cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var = "sample",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

cluster0_conserved_markers %>%  view()


annotations <- read.csv("data/annotation.csv")



# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

head(cluster0_ann_markers)


cluster0_ann_markers %>%  view()



##homework


# FindConservedMarkers for cluster 10
cluster10_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                    ident.1 = 10,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)

# Combine markers with gene descriptions 
cluster10_ann_markers <- cluster10_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# View the top 6 marker genes
head(cluster10_ann_markers)

################################################


# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}



# Iterate function across desired clusters
conserved_markers <- map_dfr(c(4,0,6,2), get_conserved)
head(conserved_markers) %>% view()


############################
cluster10_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                    ident.1 = 10,
                                                    grouping.var = "sample",
                                                    only.pos = TRUE,
                                                    logfc.threshold = 0.25)


cluster10_ann_markers <- cluster10_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# View the top 6 marker genes
head(cluster10_ann_markers)


