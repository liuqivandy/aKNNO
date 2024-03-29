#'aKNN

#' @description An Adaptive K-Nearest-Neighbors Graph (aKNN) Construction Method for  Clustering of Single Cell and Spatial transcriptomics Data
#'
#' @param obj A Seurat object
#' @param reduction Reduction to use as input to build an aKNN graph. Default is 'pca'
#' @param delta The parameter controlling the sensitivity to local distances change (<=0). Default is -0.5
#' @param dims Dimensions of reduction to use. Default is 1:50
#' @param kmax The largest k-value. Default is 20
#' @param prune.SNN The cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Default is 1/15
#'
#' @return A Seurat object with a shared-adaptive-KNN graph (named as "aKNN") stored in the respective slot.
#'
#' library(aKNNO)
#' obj<-readRDS(url("https://www.dropbox.com/s/f5khi1zperqkybg/Mouse_Brain_Serial_Section1_SagittalPosterior_rawimage_object.rds?dl=1"))
#' obj <- FindNeighbors_aKNN(obj)

#' @export

FindNeighbors_aKNN<-function(obj,reduction="pca",delta=-0.5,dims=1:50,kmax=20,prune.SNN=1/15){

  if (!reduction %in% names(obj@reductions)) {
    stop("the reduction not present in the Seurat object")
  }

  pca.use<-Embeddings(obj,reduction=reduction)

  if (max(dims) > ncol(x = pca.use)) {
    warning(paste0("More dimensions specified in dims than have been computed; use the computed dimensions:,",ncol(pca.use)))
    dims=1:ncol(pca.use)
  }

  pca.use <- pca.use[,dims]
  nnmtx<-RANN::nn2(pca.use,k = kmax+1)
  ind<-t(apply(nnmtx$nn.dists,1,function(x){x<(sum(sqrt(x))/(kmax-1-delta))^2}))
  ind<-ind[,-ncol(ind)]
  nnmtx$nn.idx<-nnmtx$nn.idx[,-(kmax+1)]

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
  snnmatrix@x<-snnmatrix@x/(kmax*2-snnmatrix@x)
  snnmatrix<-as(snnmatrix,"TsparseMatrix")
  keep<-snnmatrix@x>prune.SNN
  snnmatrix@i<-snnmatrix@i[keep]
  snnmatrix@j<-snnmatrix@j[keep]
  snnmatrix@x<-snnmatrix@x[keep]
  snn.graph<-as.Graph(snnmatrix)

  DefaultAssay(snn.graph) <- DefaultAssay(obj)
  obj@graphs$aKNN<-snn.graph

  return(obj)
}


#'FindOptimalDelta

#' @description Find the optimal delta to construct the adaptive nearest-neighbor graph
#'
#' @param obj A Seurat object
#' @param reduction Reduction to use as input to build an aKNN graph. Default is 'pca'
#' @param delta The range of delta for the grid search (<=0). Default is c(-1,0)
#' @param cutoff The cutoff for the increasing number of communities. Default is 5
#' @param dims Dimensions of reduction to use. Default is 1:50
#' @param kmax The largest k-value. Default is 20
#' @param prune.SNN The cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Default is 1/15
#' @param verbose Print progress message (TRUE or FALSE). Default is TRUE
#'
#' @return A list containing the optimal delta (opt_delta) and a ggplot showing the increased number of communities and singletons along with the decreased delta
#' @examples
#' library(aKNNO)
#' obj<-readRDS(url("https://www.dropbox.com/s/f5khi1zperqkybg/Mouse_Brain_Serial_Section1_SagittalPosterior_rawimage_object.rds?dl=1"))
#' res <- FindOptimalDelta(obj)
#' #plot the increasing number of communities/singletons with decreasing delta values
#' res$plot_delta
#' # the optimal delta is -0.8
#' res$opt_delta
#' @export


FindOptimalDelta<-function(obj,reduction="pca",delta=c(-1,0),cutoff=5,kmax=20,prune.SNN=1/15,dims=1:50,verbose=T){

  if (!reduction %in% names(obj@reductions)) {
    stop("the reduction doesn't exist in the object")
  }
  if (min(delta)>0) {stop("delta should be less than or equal to zero")}

  delta_val<-seq(delta[1],delta[2],by=0.1)
  nint<-length(delta_val)
  result<-data.frame(cbind(delta=delta_val,communities=0,singleton=0))

  for (ind in 1:nint){

    obj<-FindNeighbors_aKNN(obj,kmax=kmax,prune.SNN =prune.SNN,delta= delta_val[ind],dims=dims)
    obj<-FindClusters(obj,graph.name="aKNN",group.singletons = F,verbose=verbose)

    result[ind,3]<-sum(obj@active.ident=="singleton")
    result[ind,2]<-length(levels(droplevels(obj@active.ident,"singleton")))+result[ind,3]


  }

  diffval<-abs(diff(result[,2]))
  if (sum(diffval>cutoff)>0){
    opt_delta<-delta_val[max(which(diffval>5)+1)] } else {opt_delta=delta_val[min(which(result[,3]==0))] }
   if (is.na(opt_delta)){stop("cannot find the optimal delta, considering decrease prune.SNN")}

  plot_delta<-result %>% select(delta,communities,singleton) %>% gather(key="variable",value="Number",-delta) %>% ggplot(aes(x=delta,y=Number,col=variable))+geom_line()+geom_point()+geom_vline(xintercept=opt_delta,linetype=2)+scale_x_reverse(breaks=seq(delta[1],delta[2],0.2))+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())

 return(list(opt_delta=opt_delta,plot_delta=plot_delta))
}


#'aKNNO

#' @description An Adaptive K-Nearest-Neighbors Graph with Optimization (aKNNO) for Clustering of Single Cell and Spatial transcriptomics Data
#'
#' @param obj A Seurat object
#' @param reduction Reduction to use as input to build an aKNN graph. Default is 'pca'
#' @param delta The range of delta for the grid search (<=0). Default is c(-1,0)
#' @param cluster whether to perform clustering on the aKNNO graph (TRUE or FALSE); If TRUE, Louvain clustering at a resolution of 0.8 in the Seurat will be performed; Default is False.
#' @param dims Dimensions of reduction to use. Default is 1:50
#' @param kmax The largest k-value. Default is 20
#' @param prune.SNN The cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Default is 1/15
#' @param cutoff The cutoff for the increasing number of communities. Default is 5
#' @param verbose Print progress message (TRUE or FALSE). Default is TRUE
#'
#' @return A Seurat object with an optimized shared-adaptive-KNN graph (named as aKNN_O) stored in the respective slot.
#'
#' @examples
#' library(aKNNO)
#' obj<-readRDS(url("https://www.dropbox.com/s/f5khi1zperqkybg/Mouse_Brain_Serial_Section1_SagittalPosterior_rawimage_object.rds?dl=1"))
#' obj <- FindNeighbors_aKNNO(obj)
#' @export

FindNeighbors_aKNNO<-function(obj,cluster=FALSE, reduction="pca",delta=c(-1,0), dims=1:50,kmax=20,prune.SNN=1/15, cutoff=5,verbose=T){

  if (!reduction %in% names(obj@reductions)) {
    stop("the reduction doesn't exist in the object")
  }
  if (verbose) {message("find the optimal delta")}
  delta_res<-FindOptimalDelta(obj,reduction=reduction,delta=delta,dims=dims,kmax=kmax,prune.SNN=prune.SNN,cutoff=cutoff,verbose=verbose)
  if (verbose) {message("build the adaptive nearest-neighbor graph at the optimal delta")}
   obj<-FindNeighbors_aKNN(obj,reduction=reduction,delta=delta_res$opt_delta,dims=dims,kmax=kmax,prune.SNN=prune.SNN)

   obj[["aKNN_O"]]<-obj@graphs$aKNN
   obj@graphs$aKNN<-NULL

  obj@misc$aKNNO_delta<-delta_res
  if (cluster) { obj<-FindClusters(obj,graph.name="aKNN_O") }

  return(obj)
}
