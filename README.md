# aKNNO

==========
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Function](#Function)
* [Tutorial](#Tutorial)
* [Citation] (#Citation)

<a name="introduction"/>

# Introduction

aKNNO is a R package to build an optimized adaptive k-Nearst-Neighbor graph for single-cell and spatial transcriptomics clustering.

<a name="installation"/>

# Installation

```
devtools::install_github("liuqivandy/aKNNO")
```

<a name="Function"/>

# Function

aKNNO uses the PCA in a Seurat object to calculate distances and return a Seurat object with an optimized adaptive nearest neighbor graph (name as aKNN_O) stored in the respective slot

```R
obj <- FindNeighbors_aKNN(obj)
```


<a name="Tutorial"/>

# Tutorial


<p align="center">
  <img width="800"  src="https://github.com/JiaLiVUMC/akNN/blob/main/PBMC_toy.png">
</p>


