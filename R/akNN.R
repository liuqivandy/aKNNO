#'akNN

#' @description An Adjusted K-Nearst-Neighbors Graph Construction Method for the Clustering of Single Cell RNA Sequencing Data
#'
#' @param obj An Seurat object.
#' @param reduction Name of reduction to pull cell embeddings for. Default is 'pca'
#' @param knn Defines the largest k. Default is 20
#' @param prune The cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Default is 1/15
#' @param delta Default is -0.5
#' @param dims Dimensions of reduction to use. Default is 50
#'
#' @return This function will return a graph object name as akNN.
#'
#' @examples
#' library(akNN)
#' load(system.file("data", "PBMC_toy.Rdata", package = "akNN"))
#' obj <- FindNeighbors_akNN(obj)
#' @export

FindNeighbors_akNN<-function(obj,reduction="pca",knn=20,prune=1/15,delta=-0.5,dims=50){

  pcatmp<-Embeddings(obj,reduction=reduction)
  pcatmp <- pcatmp[,1:max(dims,dim(pcatmp)[2])]
  nnmtx<-nn2(pcatmp,k = knn+1)
  ind<-t(apply(nnmtx$nn.dists,1,function(x){x<(sum(sqrt(x))/(knn-1-delta))^2}))
  ind<-ind[,-ncol(ind)]
  nnmtx$nn.idx<-nnmtx$nn.idx[,-(knn+1)]

  rowind<-rep(nnmtx$nn.idx[,1],apply(ind,1,sum))
  colind<-t(nnmtx$nn.idx)[t(ind)]
  rownum<-nrow(nnmtx$nn.idx)

  tmpsparse<- sparseMatrix(
    i = rowind,
    j = colind,
    x = rep(1,length(rowind)),
    dims = c(rownum,rownum),
    dimnames = list(colnames(obj),colnames(obj)))

  snnmatrix<-tmpsparse %*% Matrix::t(tmpsparse)
  snnmatrix@x<-snnmatrix@x/(knn*2-snnmatrix@x)
  snnmatrix<-as(snnmatrix,"TsparseMatrix")
  keep<-snnmatrix@x>prune
  snnmatrix@i<-snnmatrix@i[keep]
  snnmatrix@j<-snnmatrix@j[keep]
  snnmatrix@x<-snnmatrix@x[keep]

  obj@graphs$akNN<-as.Graph(snnmatrix)
  #obj$klist <- rowSums(ind[,1:knn])
  return(obj)
}
