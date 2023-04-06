# akNN

==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Function](#Function)
* [Incorporate with Seurat](#incorporate_with_Seurat)

<a name="introduction"/>

# Introduction

akNN is a R package to perform an adjusted k-Nearst-Neighbors graph construction for the clustering of scRNAseq Data. 

<a name="installation"/>

# Installation

```
devtools::install_github("JiaLiVUMC/akNN")
```

<a name="Function"/>

# Function

akNN uses the PCA matrix in Seurat object as input and return a graph object name as 'akNN'.

```R
obj <- FindNeighbors_akNN(obj,
                          reduction = "pca",
                          knn=20,
                          prune=1/15,
                          delta=-0.5,
                          dims=50)
```

`obj` An Seurat object

`redunction` Name of reduction to pull cell embeddings for. Default is 'pca'

`knn` Defines the largest k. Default is 20

`prune` The cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Default is 1/15

`delta` Default is -0.5

`dims` Dimensions of reduction to use. Default is 50

<a name="incorporate_with_Seurat"/>

# Incorporate with Seurat

```R
library(akNN)
load(system.file("data", "PBMC_toy.Rdata", package = "akNN"))
obj <- FindNeighbors_akNN(obj)
obj <- FindClusters(obj,graph.name = "akNN",resolution = 0.1,verbose = F)

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


