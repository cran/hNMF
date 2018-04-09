# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################


## #' Visualization of (h)NMF results
## #' 
## #' (h)NMF abundance maps are plotted overlaid on top of a background (input) image
## #' 
## #' @param nmfInput List with NMF input attributes
## #' @param nmfModel NMF model (in accordance with NMF package definition)
## #' @return 
## #' @author Nicolas Sauwen
viewNMFResult <- function(nmfInput,nmfModel) {
  
#  library(NMF)
  
  nRows <- c(nmfInput$numRows)
  nCols <- c(nmfInput$numCols)
  nSlices <- c(nmfInput$numSlices)
  selectVect <- nmfInput$selectVect
  bgImageTensor <- nmfInput$bgImageTensor
  edgeROI <- nmfInput$edgeROI
  H <- NMF::coef(nmfModel)
  nSources <- nrow(H)
  mid <- round(nSlices/2)
  res <- nSources %% 3 # we will show the results on rows with multiples of 3
  overlayTensor <- array(0,dim=c(nRows,nCols,nSlices))
  
  if(res == 0) { # number of sources fits exactly on subfigure rows, we will not show background image
    nFigRows <- nSources %/% 3
    makeFigure(width = 10,height = 3*nFigRows)
    graphics::layout(t(matrix(c(1:nSources),nrow = 3)))
    for(iSource in 1:nSources) {
      overlayTensor[] <- H[iSource,]
      imoverlay(bgImageTensor[,,mid],overlayTensor[,,mid],selectVect[,,mid])
    }
  }
  else {
    nFigRows <- (nSources %/% 3) + 1
    makeFigure(width = 10,height = 3*nFigRows)
    graphics::layout(t(matrix(c(1:(nSources+1)),nrow = 3)))
    imoverlay(bgImageTensor[,,mid],edgeROI)
    for(iSource in 1:nSources) {
      overlayTensor[] <- H[iSource,]
      imoverlay(bgImageTensor[,,mid],overlayTensor[,,mid],selectVect[,,mid])
    }
  }
}
