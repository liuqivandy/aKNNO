# akNN

==========
* [Introduction](#introduction)
* [Installation](#installation)
* [Example](#example)
* [Incorporate with Seurat](#incorporate with Seurat)

<a name="introduction"/>

# Introduction

akNN is a R package to perform an adjusted k-Nearst-Neighbors graph construction for the clustering of scRNAseq Data. 

<a name="installation"/>

# Installation

```
devtools::install_github("JiaLiVUMC/akNN")
```

<a name="Example"/>

# Example

akNN uses the PCA matrix in Seurat object as input and return a graph object name as 'akNN'.

```R
library(akNN)
load(system.file("data", "PBMC_toy.Rdata", package = "akNN"))
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

<a name="Incorporate with Seurat"/>

# Incorporate with Seurat

```R
load(system.file("data", "PBMC_toy.Rdata", package = "akNN"))
obj <- FindNeighbors_aKNN(obj)
obj <- FindClusters(obj,graph.name = "akNN",resolution = 0.1,verbose = F)
```

