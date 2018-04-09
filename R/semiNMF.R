# Project: hNMF_git
# 
# Author: nsauwen
# Semi-NMF based on multiplicative update rules: X=AY, s.t. Y>0;
# Reference:  C. Ding, T. Li, and M.I. Jordan,
# "Convex and semi-nonnegative matrix factorizations",
# IEEE Transations on Pattern Analysis and Machine Intelligence,
# vol. 32, no. 1, pp. 45-55, 2010.
# 
###############################################################################

#' Semi-NMF based on multiplicative update rules. Reference:  C. Ding, T. Li, and M.I. Jordan,
#' "Convex and semi-nonnegative matrix factorizations",
#' IEEE Transations on Pattern Analysis and Machine Intelligence,
#' vol. 32, no. 1, pp. 45-55, 2010.
#' @param X Input data matrix, each column represents one data point 
#' and the rows correspond to the different features
#' @param nmfMod Valid NMF model, containing initialized factor matrices
#' (in accordance with the NMF package definition) 
#' @param maxiter Maximum number of iterations
#' @param checkDivergence currently not in use, to be implemented
#' @return Resulting NMF model (in accordance with the NMF package definition) 
#' @importFrom NMF basis coef
#' @importFrom MASS ginv
#' @author nsauwen
#' @export
semiNMF <- function (X,nmfMod, maxiter = 2000, checkDivergence = FALSE) {
	
	W <- NMF::basis(nmfMod)
	H <- NMF::coef(nmfMod)
	W0 <- W
	W_old <- W
	H_old <- H
	iter <- 1
	hasDiverged <- FALSE
	epsilon <- 1e-9
    
	err <- rep(0, times = maxiter+1)
	err[1] <- norm(X-W%*%H,'f')
	relTol <- 1e-5 # Relative convergence tolerance
	relError <- Inf 
	
	while(iter < maxiter & relError > relTol & !hasDiverged) {  
	  H[rowSums(H)==0,] <- epsilon
	  W_new <- X%*%MASS::ginv(H)
	  # Artificial intervention to slow down potential W divergence
	  for(iCol in 1:ncol(W_new)){
		  W_new[, iCol] <- W_new[, iCol]*max(abs(W[, iCol]))/max(abs(W_new[, iCol]))
	  }
	  W_diff <- (W - W_new)/(W + epsilon)
	  W_diff[W==0] <- 0
	  W_diff_log <- log(4*abs(W_diff)+1)/4*sign(W_diff)
	  W <- W*(1 - W_diff_log)
	  A <- t(X)%*%W
	  Ap <- (abs(A)+A)/2
	  An <- (abs(A)-A)/2
	  B <- t(W)%*%W
	  Bp <- (abs(B)+B)/2
	  Bn <- (abs(B)-B)/2
	  H <- H*sqrt((t(Ap) + t(Bn)%*%H)/(t(An) + t(Bp)%*%H))
	  
	  err[iter+1] <- norm(X-W%*%H,'f')
	  relError <- abs(err[iter+1]-err[iter])/err[iter]
	  
	  if(iter == 1 | iter == 11) { # First condition is to avoid that initial sources could be used for final result
		  W_old <- W
		  H_old <- H
	  }
	  
	  if(iter == 1 || iter%%10 == 0) {
		  print(paste("Iteration:",as.character(iter),"Relative error:",as.character(relError)))
		  if(checkDivergence & iter > 1) {
			  hasDiverged <- divergenceCheck(W, W0)
			  if(hasDiverged) { # then take NMF result from 10 iterations ago:
				  W <- W_old
				  H <- H_old
				  print("Divergence criterion reached")
			  }
		  }
	  }  
	  iter <- iter+1
	}
	NMF::basis(nmfMod) <- W
	NMF::coef(nmfMod) <- H
	return(nmfMod) 
}
