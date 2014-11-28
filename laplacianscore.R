########
# Laplacian Score
# Reference:He, Xiaofei, Deng Cai, and Partha Niyogi.
# "Laplacian score for feature selection." Advances in neural information processing systems. 2005.
#
# Usage:
# datamatrix <- iris[1:4]
# score <- laplacianScore(datamatrix, k, rho)
# 
# Argument:
# datamatrix: input matrix
# k: Simirality score of more than k-Nearest Neighbor points is 0
# rho: parameter of Gaussian Kernel
#
# Notes:
# When making Simirality Matrix undirected graph, 
# it just use a sum of the original matrix and the tranposed matrix.
# 
# Kazumasa Kaneko kanekokazumasa2@gmail.com 2014.Nov
########

laplacianScore <- function(datamatrix, k = 3,rho=0.1){
  data.dist <- as.matrix(dist(datamatrix, method="euclidean"))
  
  # Calculate Simirality matrix
  smatrix <- t(apply(data.dist,1,similarityMatrix,k,rho))
  # 双方向のリンクがあるとき二倍になっちゃう。後で修正
  smatrix <- smatrix + t(smatrix)
  data.var <- apply(datamatrix,2,variance)
 
  dmatrix <- diag(apply(smatrix,1,sum))
  lmatrix <- dmatrix - smatrix
  
  lscore <- rep(0,ncol(datamatrix))  
  one <- rep(1,nrow(datamatrix))
  for(i in 1:ncol(datamatrix)){
    fr <- datamatrix[,i]
    num <-  fr %*% dmatrix %*% as.matrix(one)
    den <- one %*% dmatrix %*% as.matrix(one)
    frt <- fr - one * (num /den)
    num <- frt %*% lmatrix %*% as.matrix(frt)
    den <- frt %*% dmatrix %*% as.matrix(frt)
    lscore[i] <- num/den
  }
  
  names(lscore) <- colnames(datamatrix)
  return(lscore)
}

similarityMatrix <- function(col,k=3,rho=0.1){
  sorted <- sort(col)
  # get k+1 points because distance from itself is 0
  k.value <- sorted[k+1]
  checked <- sapply(col,function(v,k.value){
      # distance 0 means distance from itself
      if(v<=k.value & v != 0){
        return(exp(-1.0 * v / rho ))
      }else{
        return(0)
      }
    },k.value)
  return(checked)
}

variance <- function(x){
  return(var(x)*(length(x)-1)/length(x))
}

