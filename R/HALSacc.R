# Project: hNMF_git
# 
# Author: nsauwen
# Accelerated hierarchical alternating least squares NMF. This function has been 
# written to be compatible with the 'nmf' function from the 'NMF' package. 
# It can be added to the list of NMF algorithms available to the 'nmf' function 
# using the 'SetNMFMethod' function. For a reference to the method, see N. Gillis, 
# Nonnegative matrix factorization: complexity, algorithms and applications 
# [Section 4.2, Algo. 6], PhD thesis, Université catholique de Louvain, February 2011. 
###############################################################################


#' Accelerated hierarchical alternating least squares NMF. For a reference to the method, see N. Gillis, 
#' Nonnegative matrix factorization: complexity, algorithms and applications 
#' [Section 4.2, Algo. 6], PhD thesis, Université catholique de Louvain, February 2011.
#' @param X Input data matrix, each column represents one data point 
#' and the rows correspond to the different features
#' @param nmfMod Valid NMF model, containing initialized factor matrices
#' (in accordance with the NMF package definition) 
#' @param alpha Nonnegative parameter of the accelerated method
#' @param maxiter Maximum number of iterations
#' @param checkDivergence currently not in use, to be implemented
#' @return Resulting NMF model (in accordance with the NMF package definition) 
#' @author nsauwen
#' @export
HALSacc <- function (X,nmfMod, alpha = 1, maxiter = 1000, checkDivergence = FALSE) {

  # Initialization
  nRows <- nrow(X)
  nCols <- ncol(X)
  iter <- 0
  delta <- 0.01
  relTol <- 1e-5 # Relative convergence tolerance
  relError <- Inf 
  W <- NMF::basis(nmfMod)
  H <- NMF::coef(nmfMod)
  err <- norm(X-W%*%H,'f')
  frobX <- norm(X,'f')^2
  
  while(iter < maxiter && relError > relTol) {    
    # Update of W
    A <- X%*%t(H)
    B <- H%*%t(H)
    W <- HALSupdt(t(X),t(W),t(B),t(A),alpha,delta)
    W <- t(W)
    
    # Update of H
    A <- t(W)%*%X
    B <- t(W)%*%W
    H <- HALSupdt(X,H,B,A,alpha,delta)
    
    # Evaluation of convergence
    errNew <- sqrt((frobX-2*sum(H*A))+sum(B*(H%*%t(H))))
    err <- c(err,errNew)
    relError <- abs(err[length(err)]-err[length(err)-1])/err[length(err)-1]
  
    if(iter%%10 == 0) {
      print(paste("Iteration:",as.character(iter),"Relative error:",as.character(relError)))
    }  
    iter <- iter+1
  }  
  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- H
  #nmfMod <- nmfModel(W=W,H=H) # nmfModel function from the 'NMF' package
  return(nmfMod) 
}


#' Updating step for accelerated HALS NMF
#' @param M Input data matrix
#' @param V Factor matrix to be updated
#' @param UtU Product of the other transposed factor matrix with itself
#' @param UtM Product of the other transposed factor matrix with the input matrix
#' @param alpha Nonnegative parameter of the accelerated method
#' @param delta Convergence parameter
#' @return Updated factor matrix V
#' @author Nicolas Sauwen
#' @export
HALSupdt <- function (M,V,UtU,UtM,alpha,delta) {
  
  nRows <- nrow(V)
  nCols <- ncol(V)
  cnt <- 1
  eps <- 1
  eps0 <- 1
  
  while(eps > sqrt(delta)*eps0) {
    noDelta <- 0
    for(iRow in 1:nRows) {
      deltaV <- pmax((UtM[iRow,]-UtU[iRow,]%*%V)/UtU[iRow,iRow],-V[iRow,])
      V[iRow,] <- V[iRow,] + deltaV
      noDelta <- noDelta + deltaV%*%t(deltaV)
      zeroInds <- which(V[iRow,]==0)
      V[iRow,zeroInds] <- 1e-16*max(V)
    }
    if(cnt == 1) {
      eps0 <- noDelta
    }
    eps <- noDelta
    cnt <- 0
  }  
  return(V)
}















