# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################

## #' Interactive method to select Matlab (.mat) input dataset to perform 
## #' (h)NMF analyses
## #' @return List of input data attributes
## #' @author Nicolas Sauwen
importMatlabMriData <- function() {
  
#  library(R.matlab)
#  library(spatialfil)
  
  matlabFile <- file.choose()
  matlabInputData <- R.matlab::readMat(matlabFile)
  inputDataMRI <- NULL
  inputDataMRI$data <- matlabInputData$X.red
  inputDataMRI$numRows <- matlabInputData$num.rows
  inputDataMRI$numCols <- matlabInputData$num.cols
  inputDataMRI$numSlices <- matlabInputData$num.slices
  inputDataMRI$selectVect <- matlabInputData$select.vect
  inputDataMRI$bgImageTensor <- matlabInputData$T2.data.tensor
  inputDataMRI$W <- matlabInputData$W
  inputDataMRI$H <- matlabInputData$H
  selectVect <- inputDataMRI$selectVect[,,1]
  edgeROI <- spatialfil::applyFilter(selectVect,kernel = spatialfil::convKernel(k = "sharpen"))
  edgeROI[edgeROI < max(edgeROI)-1] <- 0
  edgeROI[edgeROI > 0] <- 1
  inputDataMRI$edgeROI <- edgeROI
  
  return(inputDataMRI)
}
