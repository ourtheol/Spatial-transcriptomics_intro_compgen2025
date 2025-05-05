# Spatial-transcriptomics_intro_compgen2025
Intro to spatial omics analysis using the VoltRon and Seurat R packages on Visium and Xenium breast cancer tumor microenvironment datasets during the compgen2025 2nd module, taught by Artur Manukyan.

Datasets: 

https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast

https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-breast-cancer-fresh-frozen

This training took place in a rolv.io compute R session (R 4.4.2)

Prerequisits before installing VoltRon:

rhdf5: BiocManager::install("rhdf5") 

RBioFormats: BiocManager::install("RBioFormats") Note: RBioFormats might be hard to install due to java dependency. Check https://github.com/BIMSBbioinfo/VoltRon/tree/dev?tab=readme-ov-file#dependencies for more information. 

Seurat: install.packages("Seurat") 

spacexr: devtools::install_github("dmcable/spacexr") 

ComplexHeatmap: BiocManager::install("ComplexHeatmap") 

HDF5Array: BiocManager::install("HDF5Array") 

HDF5DataFrame: devtools::install_github("BIMSBbioinfo/HDF5DataFrame") 

ImageArray: devtools::install_github("BIMSBbioinfo/ImageArray") 

BPCells: devtools::install_github("bnprks/BPCells/r@v0.3.0") Note: You may need to install BPCells before installing VoltRon 

vitessceR: devtools::install_github("vitessce/vitessceR") 

basilisk: BiocManager::install("basilisk") 

ggnewscale: install.packages("ggnewscale") 

presto: devtools::install_github('immunogenomics/presto')


The following plots have been created by running the R script:

1. Visualization of the Visium microscopy image and the spots.
![Alt text](spacial_transcriptomics/rstudio-export/0.png)
2. The spaceranger pipeline output has been read into a Seurat object and I plot the UMAP clustering (with default parameters) and then both the spot-level expression data (number of UMIs) along with the associated image of the tissue slice. 
![Alt text](spacial_transcriptomics/rstudio-export/1.png)
![Alt text](spacial_transcriptomics/rstudio-export/2.png)
3. Molecular data on top of tissue histology. Expression of selected genes.
![Alt text](spacial_transcriptomics/rstudio-export/3.png)
4. Clustering of the seurat object in UMAP plot (after selecting the number of pcs with an elbow plot).
![Alt text](spacial_transcriptomics/rstudio-export/4.png)
5. Heatmap of the top 5 genes per cluster.
![Alt text](spacial_transcriptomics/rstudio-export/5.png)
6. Distinguishing the spatial localization of individual clusters.
![Alt text](spacial_transcriptomics/rstudio-export/6.png)
7. After deconvolution I plot selected features (clusters in this case, because I have not annotated the cell types) on the spacial image.
![Alt text](spacial_transcriptomics/rstudio-export/7.png)
8. Niche clustering.
![Alt text](spacial_transcriptomics/rstudio-export/8.png)
9. Spatial clustering plot of the niche clusters.
![Alt text](spacial_transcriptomics/rstudio-export/9.png)
10. Heatmap of the composition of seurat_clusters (should have been cell types here) in each niche cluster.
![Alt text](spacial_transcriptomics/rstudio-export/10.png)
11. UMAP Visualization of the clustering on the Xenium dataset.
![Alt text](spacial_transcriptomics/rstudio-export/11.png)
12. UMAP Visualization of the cell type annotaion on the Xenium dataset.
![Alt text](spacial_transcriptomics/rstudio-export/12.png)
13. Spatial Plot of the clusters.
![Alt text](spacial_transcriptomics/rstudio-export/13.png)
14. Spatial Plot of the cell type annotation.
![Alt text](spacial_transcriptomics/rstudio-export/14.png)
15. Heatmap of the cell types per niche cluster to examine their distance.
![Alt text](spacial_transcriptomics/rstudio-export/15.png)
16. Spacial plot of selected cell types.
![Alt text](spacial_transcriptomics/rstudio-export/16.png)
17. Hot spot analysis.
![Alt text](spacial_transcriptomics/rstudio-export/17.png)
18. Hot spot analysis around selected feature: PGR (PGR_hotspot_flag).
![Alt text](spacial_transcriptomics/rstudio-export/18.png)
19. Image alignment of the Visium and the Xenium datasets. The Visium data is the reference and the Xenium data is he query data. 
![Alt text](spacial_transcriptomics/rstudio-export/19.png)

