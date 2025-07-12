#quality control
#HBC single cellranger
#july 9 2025



## Cellranger metrics evaluation

# Load libraries
library(tidyverse)
library(reshape2)
library(plyr)

# Names of samples (same name as folders stored in data)
samples <- c("ctrl", "stim")

# Loop over each sample and read the metrics summary in
metrics <- list()
for (sample in samples) {
  path_csv <- paste0("data/", sample, "_metrics_summary.csv")
  df <- read.csv(path_csv)
  rownames(df) <- sample
  metrics[[sample]] <- df
}
# Concatenate each sample metrics together
metrics <- ldply(metrics, rbind)

# Remove periods and percentags to make the values numeric
metrics <- metrics %>%
  column_to_rownames(".id") %>%
  mutate_all(funs(parse_number(str_replace(., ",", "")))) %>%
  mutate_all(funs(parse_number(str_replace(., "%", ""))))
metrics$sample <- rownames(metrics)


# Columns of interest
cols <- c("Reads.Mapped.Confidently.to.Intergenic.Regions",
          "Reads.Mapped.Confidently.to.Intronic.Regions",
          "Reads.Mapped.Confidently.to.Exonic.Regions",
          "sample")

# Data wrangling to sculpt dataframe in a ggplot friendly manner
df <- metrics %>%
  select(cols) %>%
  melt() %>%
  mutate(variable = str_replace_all(variable, "Reads.Mapped.Confidently.to.", "")) %>%
  mutate(variable = str_replace_all(variable, ".Regions", ""))

# ggplot code to make a barplot
df %>% ggplot() +
  geom_bar(
    aes(x = sample, y = value, fill = variable),
    position = "stack",
    stat = "identity") +
  coord_flip() +
  labs(
    x = "Sample",
    y = "Percentage of Reads",
    title = "Percent of Reads Mapped to Each Region",
    fill = "Region")


########################################################################
#Loading single cell RNA-seq data into Seurat
########################################################################



## QC setup - Loading in count data

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)





# How to read in 10X data for a single sample (output is a sparse matrix)
ctrl_counts <- Read10X(data.dir = "data/ctrl_raw_feature_bc_matrix")

# Turn count matrix into a Seurat object (output is a Seurat object)
ctrl <- CreateSeuratObject(counts = ctrl_counts,
                           min.features = 100)
ctrl

# Explore the metadata
head(ctrl@meta.data)


#view(ctrl@meta.data)





###Reading in multiple samples


sample_names <- c("ctrl", "stim")

# Empty list to populate seurat object for each sample
list_seurat <- list()

for (sample in sample_names) {
  # Path to data directory
  data_dir <- paste0("data/", sample, "_raw_feature_bc_matrix")
  
  # Create a Seurat object for each sample
  seurat_data <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,
                                   project = sample)
  
  # Save seurat object to list
  list_seurat[[sample]] <- seurat_obj
}


list_seurat


#merging the control and the stimulated ones

# Create a merged Seurat object
merged_seurat <- merge(x = list_seurat[["ctrl"]], 
                       y = list_seurat[["stim"]], 
                       add.cell.id = c("ctrl", "stim"))

# Concatenate the count matrices of both samples together
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat


# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)


view(merged_seurat@meta.data)


tail(merged_seurat@meta.data)

view(merged_seurat@meta.data)


# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)



# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100


# Create metadata dataframe, just to avoid making change to the surat objet itself 
##kinnda like check point if the code below fails 
#we can have everything upto this point are preserved

metadata <- merged_seurat@meta.data

view(metadata)

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


view(metadata)




# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")




# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")




# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)


# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)



# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)



# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)



# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)


merged_seurat

# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))
filtered_seurat


##with this number of gene and mito ratio and log of gene per UMI we remove about 2000 cells 

# Extract counts
counts <- GetAssayData(object = filtered_seurat, layer = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0


view(counts)


# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

view(filtered_counts)


# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
filtered_seurat


# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data


# Visualize the number of cell counts per sample
metadata_clean %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


# Visualize the number UMIs/transcripts per cell
metadata_clean %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)


# Visualize the distribution of genes detected per cell via histogram
metadata_clean %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)



# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered.RData")



