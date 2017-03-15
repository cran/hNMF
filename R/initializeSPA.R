# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################


#' The successive projection algorithm, a useful method for initializing the NMF source matrix
#' @param data Input data matrix. The columns correspond to the data points, each row 
#' represents one feature
#' @param nSources Number of sources to be obtained
#' @return Matrix with initialized sources as its columns
#' @author Nicolas Sauwen
#' @export
initializeSPA <- function(data,nSources) {

  inds <- c()
  nCol <- ncol(data)
  dataTemp <- data
  
  if(nCol < nSources) {
    stop("nSources has to be smaller than the number of data points")
  }
  
  S <- apply(data,2,sum)
  noData <- which(S==0)
  
  if(length(noData) != 0) {
    data <- data[,-noData]
  }
  
  while(length(inds) < nSources) {
    nInds <- min(nrow(dataTemp), nSources-length(inds))
    inds <- c(inds,spaSelect(dataTemp,nInds))
    dataTemp[,inds] <- 0
  }
  
  W0 <- data[,inds]
  W0 <- matrix(W0, ncol = nSources)
  return(W0)
  
}