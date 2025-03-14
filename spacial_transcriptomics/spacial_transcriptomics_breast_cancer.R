options(java.parameters = "-Xmx8g")
library(VoltRon)

####
# Visium Analysis ####
####

####
## Import Data in VoltRon ####
####

# Dependencies
if(!requireNamespace("rhdf5"))
  BiocManager::install("rhdf5")
library(rhdf5)

# import Breast Cancer Visium data
bc_visium <- importVisium("/home/rolv-user/workshop/data/BreastCancer/Visium/",
                          sample_name = "bc_visium")

# sample metadata
SampleMetadata(bc_visium)

# metadata
Metadata(bc_visium)

# features 
vrFeatures(bc_visium)

####
# Images ####
####

# images and channels
vrImageChannelNames(bc_visium)

# get images
vrImages(bc_visium)
vrImages(bc_visium, scale.perc = 20)
vrImages(bc_visium, channel = "H&E", scale.perc = 20)

####
## Omic Profile Clustering ####
####

####
### Processing and Filtering ####
####

# features
head(vrFeatures(bc_visium))
length(vrFeatures(bc_visium))

# normalize and select the top 3000 highly variable features
bc_visium <- normalizeData(bc_visium)
bc_visium <- getFeatures(bc_visium, n = 3000)  

# selected features
head(vrFeatureData(bc_visium))
selected_features <- getVariableFeatures(bc_visium)
head(selected_features, 20)

####
### Dimensionality Reduction ####
####

# embedding
bc_visium <- getPCA(bc_visium, features = selected_features, dims = 30)
bc_visium <- getUMAP(bc_visium, dims = 1:30)
vrEmbeddingNames(bc_visium)

# embedding visualization
vrEmbeddingPlot(bc_visium, embedding = "umap")

####
### Clustering ####
####

# Clustering of the Visium spots

# graph for neighbors
bc_visium <- getProfileNeighbors(bc_visium, dims = 1:30, k = 10, method = "SNN")
vrGraphNames(bc_visium)

# clustering
bc_visium <- getClusters(bc_visium, resolution = 0.5, label = "Clusters", graph = "SNN")

####
### Visualization ####
####

# embedding
vrEmbeddingPlot(bc_visium, embedding = "umap", group.by = "Clusters")
vrSpatialPlot(bc_visium, group.by = "Clusters", plot.segments = TRUE)
vrSpatialPlot(bc_visium, group.by = "Clusters", plot.segments = TRUE, alpha = 0.5)

####
## Import Data in Seurat ####
####

# Load a 10x Genomics Visium Spatial Experiment into a Seurat object
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(dplyr)

bc_visium_seu2 <- Load10X_Spatial(data.dir = "/home/rolv-user/workshop/data/BreastCancer/Visium/",
                                 filename = "CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5")

plot1 <- VlnPlot(bc_visium_seu2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(bc_visium_seu2, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

#  increasing the future.globals.maxSize
options(future.globals.maxSize = 1000 * 1024^2)  # Increase limit to 1 GB
bc_visium_seu2 <- SCTransform(bc_visium_seu2, assay = "Spatial", verbose = FALSE)

# overlay molecular data on top of tissue histology
SpatialFeaturePlot(bc_visium_seu2, features = c("COL1A1","ERBB2","MIEN1"), slot = "counts")

# Dimensionality reduction, clustering, and visualization
bc_visium_seu2 <- RunPCA(bc_visium_seu2, assay = "SCT", verbose = FALSE)

ElbowPlot(bc_visium_seu2)

bc_visium_seu2 <- FindNeighbors(bc_visium_seu2, reduction = "pca", dims = 1:15)
bc_visium_seu2 <- FindClusters(bc_visium_seu2, verbose = FALSE)
bc_visium_seu2 <- RunUMAP(bc_visium_seu2, reduction = "pca", dims = 1:15)

p1 <- DimPlot(bc_visium_seu2, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(bc_visium_seu2, label = TRUE, label.size = 3)
p1 + p2

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(bc_visium_seu2, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5
DoHeatmap(bc_visium_seu2, features = top5$gene) + NoLegend()


#  highlight particular cells of interest, for distinguishing the spatial localization of individual clusters
SpatialDimPlot(bc_visium_seu2, cells.highlight = CellsByIdentities(object = bc_visium_seu2, idents = c(0, 2, 6, 3, 5, 14)), facet.highlight = TRUE, ncol = 3)

####
### Deconvolution ####
####

# This part would have been more interesting if we had annotated the clusters according to cell types
# to understand the proportion of cellular mixtures within each spot
# for now we do this on the proportion of clusters

# visualize
Idents(bc_visium_seu2) <- "seurat_clusters"
gsubclass <- DimPlot(bc_visium_seu2, reduction = "umap", label = T) + NoLegend()

# Deconvolute anterior and posterior Visium spots
if(!requireNamespace("spacexr"))
  devtools::install_github("dmcable/spacexr")
library(spacexr)

# deconvolute spots
bc_visium <- getDeconvolution(bc_visium, sc.object = bc_visium_seu2, sc.assay = "Spatial",
                               sc.cluster = "seurat_clusters", max_cores = 2)

# Visualize
vrMainFeatureType(bc_visium) <- "Decon"
vrFeatures(bc_visium)
vrSpatialFeaturePlot(bc_visium, features = c("10", "6", "14", "11"),
                     crop = TRUE, ncol = 2, alpha = 1, keep.scale = "all")

####
### Processing ####
####

vrMainFeatureType(bc_visium) <- "Decon"
bc_visium <- normalizeData(bc_visium, method = "CLR")

# embedding visualization
bc_visium <- getUMAP(bc_visium, data.type = "norm", umap.key = "umap_niche")
vrEmbeddingPlot(bc_visium, embedding = "umap_niche", group.by = "Sample")

####
### Niche Clustering ####
####

# clustering
bc_visium <- getProfileNeighbors(bc_visium, data.type = "norm", method = "SNN", graph.key = "SNN_niche")
bc_visium <- getClusters(bc_visium, resolution = 0.4, graph = "SNN_niche", label = "Niche_Clusters")

# visualize clustering
g1 <- vrEmbeddingPlot(bc_visium, embedding = "umap", group.by = "Sample")
g2 <- vrEmbeddingPlot(bc_visium, embedding = "umap", group.by = "Niche_Clusters", label = TRUE)
g1 | g2

####
### Visualization ####
####

# spatial clustering plot
vrSpatialPlot(bc_visium, group.by = "Niche_Clusters", crop = TRUE, alpha = 1)

# heatmap plot
if(!requireNamespace("ComplexHeatmap"))
  BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
vrHeatmapPlot(bc_visium, features = vrFeatures(bc_visium), group.by = "Niche_Clusters",
              show_row_names = T, show_heatmap_legend = T)




####
# Xenium Analysis ####
####

# Dependencies
if(!requireNamespace("rhdf5"))
  BiocManager::install("rhdf5")
if(!requireNamespace("RBioFormats"))
  BiocManager::install("RBioFormats")
library(rhdf5)

# import data
Xen_R1 <- importXenium("/home/rolv-user/workshop/data/BreastCancer/Xenium_R1/outs/", 
                       sample_name = "XeniumR1", 
                       overwrite_resolution = TRUE, #regenerate tissue image from the pyramide image of the Xenium data
                       resolution_level = 3)

####
## Filtering ####
####

# Metadata
Xen_R1@metadata
head(Metadata(Xen_R1))

# filter out counts - keep cells with more than 5 transcripts
Xen_R1 <- subset(Xen_R1, Count > 5)

####
## Disk Representation ####
####

# Representation to store large datasets

# Dependencies
if(!requireNamespace("rhdf5"))
  BiocManager::install("rhdf5")
if(!requireNamespace("HDF5Array"))
  BiocManager::install("HDF5Array")
if(!requireNamespace("ImageArray"))
  devtools::install_github("BIMSBbioinfo/ImageArray")
if(!requireNamespace("BPCells"))
  devtools::install_github("bnprks/BPCells/r")
library(rhdf5)
library(HDF5Array)
library(ImageArray)
library(BPCells)

# save voltron on disk
Xen_R1_ondisk <- saveVoltRon(Xen_R1, 
                             format = "HDF5VoltRon", 
                             output = "/home/rolv-user/workshop/data/ondisk/Xen_R1", 
                             replace = TRUE)
# Xen_R1 takes up 1.2 GB of memory, has been reduced to 642 MB

# load voltron from disk
Xen_R1_ondisk <- loadVoltRon("/home/rolv-user/workshop/data/ondisk/Xen_R1/")

####
## Omic Profile Clustering ####
####

# normalize
Xen_R1_ondisk <- normalizeData(Xen_R1_ondisk, sizefactor = 1000)

# PCA reduction
Xen_R1_ondisk <- getPCA(Xen_R1_ondisk, dims = 20, overwrite = TRUE)
Xen_R1_ondisk <- getUMAP(Xen_R1_ondisk, dims = 1:20, overwrite = TRUE)
vrEmbeddingNames(Xen_R1_ondisk)

vrEmbeddingPlot(Xen_R1_ondisk, embedding = "pca")
vrEmbeddingPlot(Xen_R1_ondisk, embedding = "umap")

# neighbors
Xen_R1_ondisk <- getProfileNeighbors(Xen_R1_ondisk, dims = 1:20, method = "SNN")
vrGraphNames(Xen_R1_ondisk)

# clustering
Xen_R1_ondisk <- getClusters(Xen_R1_ondisk, resolution = 1.3, label = "Clusters", graph = "SNN")

# visualization
vrEmbeddingPlot(Xen_R1_ondisk, group.by = "Clusters", embedding = "umap", 
                pt.size = 0.4, label = TRUE)

# spatial plot
vrSpatialPlot(Xen_R1_ondisk, group.by = "Clusters", pt.size = 0.18)

####
### Marker Analysis ####
####
# find genes expressed more in each cluster to understand cell types

# convert the voltron object into a Seurat object to do stat analysis:
# what genes are diferentially expressed as oposed to other clusters
Xen_R1$Clusters <- Xen_R1_ondisk$Clusters
Xen_R1_seu <- VoltRon::as.Seurat(Xen_R1, cell.assay = "Xenium", type = "image")
Idents(Xen_R1_seu) <- "Clusters"
Xen_R1_seu <- NormalizeData(Xen_R1_seu, scale.factor = 1000)
markers <- FindAllMarkers(Xen_R1_seu) # compares the expression of each individual gene of of each individual cluster against all clusters

####
### Annotation ####
####

# get predefined annotations
annotations <- read.table("/home/rolv-user/workshop/data/BreastCancer/Xenium_R1/annotation.txt")[,1]

# annotate clusters
CellType <- factor(Xen_R1_ondisk$Clusters)
levels(CellType) <- annotations
Xen_R1_ondisk$CellType <- as.character(CellType)

# visualization
vrSpatialPlot(Xen_R1_ondisk, group.by = "CellType", pt.size = 0.18, alpha = 1)
vrEmbeddingPlot(Xen_R1_ondisk, group.by = "CellType", embedding = "umap",
                pt.size = 0.4, label = TRUE)

####
## Spatially aware clustering ####
####

####
### Niche Clustering ####
####

# spatial neighbors
Xen_R1_ondisk <- getSpatialNeighbors(Xen_R1_ondisk, radius = 30, method = "radius")
vrGraphNames(Xen_R1_ondisk)

# get niche assay
Xen_R1_ondisk <- getNicheAssay(Xen_R1_ondisk, label = "CellType", graph.type = "radius")
Xen_R1_ondisk # RNA is the main feature

# normalizing niche assay 
vrMainFeatureType(Xen_R1_ondisk) <- "Niche"  # now Niche is the main feature
vrFeatures(Xen_R1_ondisk) # the cell types we set before
Xen_R1_ondisk <- normalizeData(Xen_R1_ondisk, method = "CLR")

# clustering niches
Xen_R1_ondisk <- getClusters(Xen_R1_ondisk, nclus = 9, method = "kmeans", label = "niche_clusters")

# visualization
vrSpatialPlot(Xen_R1_ondisk, group.by = "niche_clusters", alpha = 1)
library(ComplexHeatmap)
vrHeatmapPlot(Xen_R1_ondisk, features = vrFeatures(Xen_R1_ondisk), group.by = "niche_clusters")

# visualization of specific cell type
vrSpatialPlot(Xen_R1_ondisk, group.by = "CellType", pt.size = 0.18, alpha = 1, group.ids = c("ACTA2_myoepithelial", "KRT15_myoepithelial"))
vrSpatialPlot(Xen_R1_ondisk, group.by = "CellType", pt.size = 1, alpha = 1, group.ids = c("CD4_TCells", "CD8_TCells", "BCells"), n.tile = 400)

####
### Hot Spot Analysis ####
####

# get spatial neighbor plot
Xen_R1_ondisk <- getSpatialNeighbors(Xen_R1_ondisk, method = "radius", radius = 15, graph.key = "radius_hot")

# visualize 
vrMainFeatureType(Xen_R1_ondisk) <- "RNA"
vrSpatialFeaturePlot(Xen_R1_ondisk, features = "PGR", alpha = 1, background.color = "black", n.tile = 300)

# analysis, find hot spots around a selected feature, search for positions around that gene being expressed
Xen_R1_ondisk <- getHotSpotAnalysis(Xen_R1_ondisk, features = "PGR", graph.type = "radius_hot", alpha.value = 0.001)

# visualize
vrSpatialFeaturePlot(Xen_R1_ondisk, features = "PGR_hotspot_stat", alpha = 1, background.color = "black", n.tile = 400)
vrSpatialPlot(Xen_R1_ondisk, group.by = "PGR_hotspot_flag", alpha = 1, background.color = "black", n.tile = 400)








# (OPTIONAL) import H&E data, and save to disk
Xen_R1_image <- importImageData("/home/rolv-user/workshop/data/BreastCancer/Xenium_R1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_highres.tif",
                                sample_name = "HEimage",
                                channel_names = "H&E", tile.size = 100)

# Store the image as a pyramide image on the disk
Xen_R1_image_disk <- saveVoltRon(Xen_R1_image, 
                                 format = "HDF5VoltRon", 
                                 output = "/home/rolv-user/workshop/data/ondisk/Xen_R1_image/", 
                                 replace = TRUE)

# load H&E data from disk
Xen_R1_image_disk <- loadVoltRon("/home/rolv-user/workshop/data/ondisk/Xen_R1_image/")






# Align 

# alignment
xen_reg <- registerSpatialData(object_list = list(bc_visium, Xen_R1_ondisk))

# transfer image
Visium_reg <- xen_reg$registered_spat[[1]]
vrImages(Visium_reg)
Xenium_reg <- xen_reg$registered_spat[[2]]
vrImages(Xenium_reg)









