# Project: hNMF_git
# 
# Author: nsauwen
# NMF by alternative non-negative least squares using projected gradients. This function has been 
# written to be compatible with the 'nmf' function from the 'NMF' package. 
# It can be added to the list of NMF algorithms available to the 'nmf' function 
# using the 'SetNMFMethod' function. For a reference to the method, see C.-J. Lin, 
# "Projected Gradient Methods for Non-negative Matrix Factorization", 
# Neural computation 19.10 (2007): 2756-2779.
###############################################################################


#' NMF by alternating non-negative least squares using projected gradients. 
#' For a reference to the method, see C.-J. Lin, 
#' "Projected Gradient Methods for Non-negative Matrix Factorization", 
#' Neural computation 19.10 (2007): 2756-2779.
#' @param X Input data matrix, each column represents one data point 
#' and the rows correspond to the different features
#' @param nmfMod Valid NMF model, containing initialized factor matrices
#' (in accordance with the NMF package definition) 
#' @param tol Tolerance for a relative stopping condition
#' @param maxIter Maximum number of iterations
#' @param timeLimit Limit of time duration NMF analysis
#' @param checkDivergence Boolean indicating whether divergence checking should be performed
#' Default is TRUE, but it should be set to FALSE when using random initialization
#' @return Resulting NMF model (in accordance with the NMF package definition) 
#' @importFrom NMF basis coef
#' @author nsauwen
#' @export
PGNMF <- function (X, nmfMod, tol = 1e-5, maxIter = 500, timeLimit = 300, checkDivergence = TRUE) {
  
  # Initialization
  startTime <- proc.time()[3]
  W <- NMF::basis(nmfMod)
  H <- NMF::coef(nmfMod)
  err <- norm(X-W%*%H,'f')
  err_diff <- Inf
  pNorm <- c()
  
  # Keep initial W0 for divergence criterion
  W0 <- W
  W_old <- W
  H_old <- H
  hasDiverged <- FALSE
  
  gradW <- W%*%(H%*%t(H)) - X%*%t(H)
  gradH <- (t(W)%*%W)%*%H - t(W)%*%X
  initGrad <- norm(rbind(gradW, t(gradH)),'f')
  
  print(paste("Init gradient norm:",as.character(initGrad)))
  
  tolW <- max(0.001, tol) * initGrad
  tolH <- tolW
  
  for(iter in 1:maxIter) {
    projNorm <- norm(c(gradW[gradW<0 | W>0], gradH[gradH<0 | H>0]), type = "2")
    pNorm <- c(pNorm,projNorm)
  
    if(err_diff < tol | proc.time()[3]-startTime > timeLimit | hasDiverged)  break
    
    nlsOutputW <- nlsSubProb(t(X), t(H), t(W), tolW, 500)
    W <- t(nlsOutputW[[1]])
    gradW <- t(nlsOutputW[[2]])
    iterW <- nlsOutputW[[3]]
    
    if(iterW == 1) {
      tolW <- 0.1*tolW
    }
    
    nlsOutputH <- nlsSubProb(X, W, H, tolH, 500)
    H <- nlsOutputH[[1]]
    gradH <- nlsOutputH[[2]]
    iterH <- nlsOutputH[[3]]
    
    if(iterH == 1) {
      tolH <- 0.1*tolH
    }
    
    err <- c(err,norm(X-W%*%H,'f'))
    err_diff <- abs(err[length(err)] - err[length(err)-1]) / err[length(err)-1]
    
    if(iterW == 1 | iterH ==1) {
      err_diff <- tol
    }
    
    if(iter == 2 | iter == 12) { # First condition is to avoid that initial sources could be used for final result
      W_old <- W
      H_old <- H
    }
    
    if(iter%%10 == 0) {
      print(paste("Relative error:", as.character(err_diff), "after", as.character(iter), "iterations"))
      if(checkDivergence) {
        hasDiverged <- divergenceCheck(W, W0)
        if(hasDiverged) { # then take NMF result from 10 iterations ago:
          W <- W_old
          H <- H_old
          print("Divergence criterion reached")
        }
#        W_old <- W
#        H_old <- H
      }
    }  
  }
  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- H
  return(nmfMod) 
}


#' Algorithm for solving convex non-negative least squares subproblem using projected gradients 
#' 
#' @param X Input data matrix
#' @param W NMF basis matrix
#' @param Hinit Initial NMF coef matrix
#' @param tol Tolerance for a relative stopping condition
#' @param maxIter Maximum number of iterations
#' @return List containing updated H matrix, its gradient and number of iterations
#' @author Nicolas Sauwen
#' @export
nlsSubProb <- function(X, W, Hinit, tol, maxIter) {
  H <- Hinit
  WtX <- t(W)%*%X
  WtW <- t(W)%*%W
  
  alpha <- 1
  beta <- 0.1
  nlsOutput <- list()
  
  for(iter in 1:maxIter) {
    grad <- WtW%*%H - WtX
    projgrad <- norm(grad[grad<0 | H>0], type = "2")
    if(projgrad < tol)  break
    
    # search step size
    for(innerIter in 1:20) {
      Hn <- pmax(H - alpha*grad, 0)
      d <- Hn - H
      gradd <- sum(grad*d)
      dQd <- sum((WtW%*%d)*d)
      suff_decr <- (0.99*gradd + 0.5*dQd) < 0
      if(innerIter == 1) {
        decr_alpha <- !suff_decr
        Hp <- H
      }
      if(decr_alpha) {
        if(suff_decr) {
          H <- Hn
          break
        }
        else {
          alpha <- alpha * beta
        }
      }
      else {
        if(!suff_decr | sum(Hp-Hn) == 0) {
          H <- Hp
          break
        }
        else{
          alpha <- alpha / beta
        }
      }
    }
  }
  if(iter == maxIter) print("Max iter in nlsSubProb")
  nlsOutput[[1]] <- H
  nlsOutput[[2]] <- grad
  nlsOutput[[3]] <- iter
  return(nlsOutput)
}

#' This function performs a divergence check, by comparing the current 
#' NMF sources with the initial ones. 3 divergence criteria are implemented.
#' 
#' @param W Current NMF source matrix
#' @param W0 Initial NMF source matrix
#' @return Boolean value, indicating whether or not one of the 
#' divergence criteria has been reached
#' @author Nicolas Sauwen
#' @export
divergenceCheck <- function(W, W0) {
  
  checkDivergence <- FALSE
  
  # Normalize columns of W and W0:
  N0 <- diag(sqrt(t(W0)%*%W0))
  N_W0 <- matrix(N0,nrow = nrow(W0),ncol = ncol(W0),byrow = T)
  W0 <- W0/N_W0
  N <- diag(sqrt(t(W)%*%W))
  N_W <- matrix(N,nrow = nrow(W),ncol = ncol(W),byrow = T)
  W <- W/N_W
  
  S_W0W0 <- sqrt(t(W0)%*%W0)
  S_WW <- sqrt(t(W)%*%W)
  S_WW0 <- sqrt(t(W)%*%W0)
  
  # Divergence check1: non-diagonal correlation values should not exceed a certain threshold 
  thr_check1 <- 0.97
  nonDiags0 <- S_W0W0[upper.tri(S_W0W0)]
  nonDiags <- S_WW[upper.tri(S_WW)]
  check1_vect <- rep(thr_check1,length(nonDiags))
  check1_vect[nonDiags0 > thr_check1] <- nonDiags0[nonDiags0 > thr_check1] + (1 - nonDiags0[nonDiags0 > thr_check1])/2
  check1 <- which(nonDiags > check1_vect)
  if(length(check1) > 0) {
    checkDivergence <- TRUE
    return(checkDivergence)
  }
  
  # Divergence check2: Correlation values should row-wise be highest on diagonal:
  maxPerRow <- apply(S_WW0,1,which.max)
  check2_vect <- c(1:nrow(S_WW0))
  if(!identical(maxPerRow, check2_vect)){
    checkDivergence <- TRUE
    return(checkDivergence)
  }
  
  # Divergence check3: Sources should not diverge too much from initial sources
  thr_check3 <- 0.90
  check3 <- which(diag(S_WW0) < thr_check3)
  if(length(check3) > 0) {
    checkDivergence <- TRUE
    return(checkDivergence)
  }
  
  return(checkDivergence)
  
}
  
  