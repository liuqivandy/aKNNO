# aKNNO

==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Function](#Function)
* [Incorporate with Seurat](#incorporate_with_Seurat)

<a name="introduction"/>

# Introduction

aKNNO is a R package to build an optimized adaptive k-Nearst-Neighbors graph for single-cell and spatial transcriptomics clustering.

<a name="installation"/>

# Installation

```
devtools::install_github("liuqivandy/aKNNO")
```

<a name="Function"/>

# Function

aKNN uses the PCA matrix in Seurat object as input and return a seurat object with an adaptive nearest neighbor graph (name as aKNN) stored in the respective slot

```R
obj <- FindNeighbors_aKNN(obj,
                          reduction = "pca",
                          kmax=20,
                          prune.SNN=1/15,
                          delta=-0.5,
                          dims=1:50)
```

`obj` A Seurat object

`redunction` Reduction to calculate the distance. Default is 'pca'

`kmax` Defines the largest k. Default is 20

`prune.SNN` The cutoff for adjusted Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Default is 1/15

`delta` The similarity weight to adjust the degree of changes in ascending distance distribution. Default is -0.5

`dims` Dimensions of reduction to use. Default is 1:50

<a name="incorporate_with_Seurat"/>

# Incorporate with Seurat

```R
library(akNN)
load(system.file("data", "PBMC_toy.Rdata", package = "akNN"))
obj <- FindNeighbors_aKNN(obj)
obj <- FindClusters(obj,graph.name = "aKNN",resolution = 0.1,verbose = F)

#corrsponding kNN construction method
obj <- FindNeighbors(obj,nn.method = "rann",dims = 1:50,verbose = F)
obj <- FindClusters(obj,graph.name = "RNA_snn",verbose = F,resolution = 0.1)

#plot
library(ggplot2)
p1 <- DimPlot(obj,label = T,repel=T,reduction = "pca",group.by="Ref")+ggtitle("Ref")+NoLegend()
p2 <- DimPlot(obj,label = T,repel=T,reduction = "pca",group.by="akNN_res.0.1")+ggtitle("akNN")+NoLegend()
p3 <- DimPlot(obj,label = T,repel=T,reduction = "pca",group.by="RNA_snn_res.0.1")+ggtitle("kNN")+NoLegend()
p1+p2+p3
```

<p align="center">
  <img width="800"  src="https://github.com/JiaLiVUMC/akNN/blob/main/PBMC_toy.png">
</p>


