# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################


#' Computation of relative NMF residual per observation
#' @param X Input data matrix, each column represents one observation 
#' @param nmfFit NMF model fitted to the input data in X
#' @return Relative residual per observation, returned as a vector
#' @author nsauwen
#' @importFrom NMF basis coef
#' @export
residualNMF <- function(X, nmfFit) {
	
	H <- NMF::coef(nmfFit)
	W <- NMF::basis(nmfFit)
	residuMat <- X-W%*%H
	absRes <- apply(residuMat^2, 2, sum)
	absX <- apply(X^2, 2, sum)
	relRes <- absRes/absX
	
	return(relRes)
	
}
