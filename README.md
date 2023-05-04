# aKNNO

==========
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Function](#Function)
* [Tutorial](#Tutorial)
* [Citation](#Citation)

<a name="Introduction"/>

# Introduction

aKNNO is a R package to build an optimized adaptive k-nearst neighbor graph for single-cell and spatial transcriptomics clustering.

<a name="Installation"/>

# Installation

```
devtools::install_github("liuqivandy/aKNNO")
```

<a name="Function"/>

# Function

aKNNO uses the PCA in a Seurat object to calculate distances and return a Seurat object with an optimized adaptive k-nearest neighbor graph (named as aKNN_O) stored in the respective slot

```R
obj <- FindNeighbors_aKNNO(obj)
```


<a name="Tutorial"/>

# Tutorial

- [Analyze Mouse Instestinal Epithelium scRNAseq data](https://htmlpreview.github.io/?https://github.com/liuqivandy/aKNNO/blob/master/Tutorial/mouseInstestine.html)
- [Analyze 10x Visium data from mouse brain sagittal posterior](https://htmlpreview.github.io/?https://github.com/liuqivandy/aKNNO/blob/master/Tutorial/mousebrain_SagittalPosterior.html)
- [Analyze Human Pancreas scRNAseq data](https://htmlpreview.github.io/?https://github.com/liuqivandy/aKNNO/blob/master/Tutorial/humanPancreas.html)
- [Analyze Mouse Intestinal Organoid scRNAseq data](https://htmlpreview.github.io/?https://github.com/liuqivandy/aKNNO/blob/master/Tutorial/mouseInstestineOrganoids.html)
- [Analyze Mouse Habenula scRNAseq data](https://htmlpreview.github.io/?https://github.com/liuqivandy/aKNNO/blob/master/Tutorial/mouseHabenula.html.html)


<a name="Citation"/>

# Citation
Jia Li, Yu Shyr, Qi Liu. Clustering of Single-cell and Spatial Transcriptomics with an Optimized Adaptive K-nearest Neighbor Graph
