# Rscript parithera.R ../files/a649e326-8ef1-4eb4-b185-bcb369c7e0b1/60adc5e8-cfc0-42e0-ba83-208fa3936f15 patient2_tube3_filtered_feature_bc_matrix.h5
# Load required libraries
library(Seurat)
library(hdf5r)
library(tidyverse)
library(devtools)
library(ggplot2)
library(rjson)

# Install 'presto' package from GitHub
devtools::install_github("immunogenomics/presto")

args <- commandArgs(trailingOnly = TRUE)

########## Create a Seurat object from H5 file ##########

# Set working directory
setwd(args[1])

# Load data for a single tube
seurat_data <- Read10X_h5(args[2])

# Create a Seurat object
seurat_object <- CreateSeuratObject(
    counts = seurat_data,
    project = "patient3_Tube1"
)

# Check objects
head(seurat_object)

########## Visualize reads, counts, and mitochondrial ##########
seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")

# Make violin plot, features = genes detected, counts = total molecules detected
plot_Vln <- VlnPlot(seurat_object, c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)
# ggsave("plot_violin.svg", plot = plot_Vln, device = "svg")
ggsave("plot_violin.png", plot = plot_Vln)


########## Normalize the data and find variable features ##########

final_seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
final_seurat_object <- FindVariableFeatures(final_seurat_object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(final_seurat_object), 10)

# Plot variable features with and without labels
plot_VFP <- VariableFeaturePlot(final_seurat_object)
plot_LP <- LabelPoints(plot = plot_VFP, points = top10, repel = TRUE)
plot_VFP + plot_LP

# ggsave("plot_variable_features.svg", plot = plot_VFP + plot_LP, device = "svg")
ggsave("plot_variable_features.png", plot = plot_VFP + plot_LP)

########## Scale data ##########

final_seurat_object <- ScaleData(final_seurat_object)

# Regress out some parameters, like in this case the percentage of MT reads
final_seurat_object <- ScaleData(final_seurat_object, vars.to.regress = "percent.mt")

########## Run PCA ##########

final_seurat_object <- RunPCA(final_seurat_object, features = VariableFeatures(object = final_seurat_object))

# Determine the dimensionality of the dataset
plot_elbow <- ElbowPlot(final_seurat_object)
# ggsave("plot_elbow.svg", plot = plot_elbow, device = "svg")
ggsave("plot_elbow.png", plot = plot_elbow)

####### CLUSER ANALYSIS
####### Find neighbors and clusters ##########
# Using 16 PCs for clustering, as it seems to be the inflection point
final_seurat_object <- FindNeighbors(final_seurat_object, dims = 1:16)
final_seurat_object <- FindClusters(final_seurat_object, resolution = 0.5)

###### Run UMAP for visualization
final_seurat_object <- RunUMAP(final_seurat_object, dims = 1:16)
plot_umap <- DimPlot(final_seurat_object, reduction = "umap", label = TRUE)
# ggsave("plot_umap.svg", plot = plot_umap, device = "svg")
ggsave("plot_umap.png", plot = plot_umap)


###### Describe number of cells in each cluster
# Get the cluster assignments from the Seurat object
cluster_assignments <- Idents(final_seurat_object)

# Create a table with the number of cells per cluster
cluster_counts <- table(cluster_assignments)

# Convert the table to a data frame
cluster_counts_df <- as.data.frame(cluster_counts)

# Rename the columns
colnames(cluster_counts_df) <- c("Cluster", "Number_of_Cells")

# Calculate the total number of cells
total_cells <- sum(cluster_counts_df$Number_of_Cells)

# Print the table
print(cluster_counts_df)
cat("Total number of cells:", total_cells, "\n")

# Save final_seurat_object to a file
saveRDS(final_seurat_object, file = "final_seurat_object.rds")
