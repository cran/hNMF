# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################


spaSelect <- function(data,nInds) {
  
  nRow <- nrow(data)
  nCol <- ncol(data)
  
  normData <- colSums(data^2)
  maxNorm <- max(normData)
  
  iter <- 1
  inds <- c()
  U <- matrix(0,nRow,nInds)
  
  while(iter <= nInds && max(normData)/maxNorm > 1e-9) {
    # Select the column of data with the largest l2-norm
    newMaxNorm <- max(normData)
    #newMaxInd <- which.max(normData)
    
    if(iter == 1) {
      origNormData <- normData
    }
    
    # Check for ties up to a precision of 1e-6
    select <- which((newMaxNorm-normData)/newMaxNorm < 1e-6)
    # In case of a tie, select the column with the largest norm of the input matrix
    if(length(select) > 1) {
      origMax <- which.max(origNormData(select))
      select <- select(origMax)
    }
    # Add index to set
    inds[iter] <- select
    U[,iter] <- data[,select]
    
    if(iter > 1) {
      for(j in 1:(iter-1)) {
        U[,iter] <- U[,iter] - U[,j]*(t(U[,j])%*%U[,iter])
      }       
    }
    # Normalization
    U[,iter] <- U[,iter]/norm(matrix(U[,iter]),'2')
    
    v <- U[,iter]
    if(iter > 1) {
      for(j in (iter-1):1) {
        v <- v - (t(v)%*%U[,j])*U[,j] 
      }
    }
    # Update the norm of the columns after orthogonal projection
    normData <- normData - (t(v)%*%data)^2
    iter <- iter+1
  } 
  return(inds)
}


















