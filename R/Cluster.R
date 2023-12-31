#' Cluster things
#'
#' @param Data XX   
#'
#' @return
#' @export
#'
#' @examples
library(cluster)
library(dynamicTreeCut)
library(bestNormalize)
library(tidyverse)

Cluster <- function(Data,min_cluster_size){
  
  standardise <- function(Data){
    b <- Data
    b <- b[,1:ncol(b)-1]
    Chl <- Data[,ncol(Data)]
    b[b==0] <- 1e-6
    b <- b/Chl
    v <- lapply(b,bestNormalize::boxcox)
    return(v)
  }
  
  
  v <- standardise(Data)
  
  L <- length(v)
  
  ndf <- list() 
  for (i in 1:L){
    ndf[[length(ndf)+1]] <- data.frame(v[[i]][1])
  }
  
  S <- Data
  ndf <- do.call("cbind", ndf)
  colnames(ndf) <- colnames(S[,ncol(S)-1])
  
  
  mscluster <- dist(ndf, method = "manhattan")
  mv.hclust <- hclust(mscluster, method = "ward.D2")
  
  plot(mv.hclust)
  
  ev.clust <- Data
  #Change the minClusterSize argument below to change the number of clusters. 
  #You might want to play around with it if you want ~ 6 clusters. 
  #You could try setting it at 1/6th of the total sample number :)
  dynamicCut <- dynamicTreeCut::cutreeDynamic(mv.hclust, cutHeight = 70, minClusterSize=min_cluster_size, method="hybrid", distM=as.matrix(dist(ndf, method="manhattan")), deepSplit=4, pamStage = TRUE, pamRespectsDendro = TRUE,
                              useMedoids = FALSE, maxDistToLabel = NULL,
                              maxPamDist = 50,
                              respectSmallClusters = TRUE,)
  
  ev.clust$Clust <- dynamicCut
  
  L2 <- length(unique(ev.clust$Clust))
  L <- list()
  for (i in 1:L2){
    L[[length(L)+1]] <- dplyr::filter(ev.clust,Clust==i)
  }
  
  L
  
  e <- numeric()
  for (i in 1:length(L)){
    e[[length(e)+1]] <- length(L[[i]][[1]])
  }
  
  
  return(list(L,plot(mv.hclust)))
}
