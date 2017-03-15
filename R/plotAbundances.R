# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################


plotAbundances <- function(nmfMod,source, slice=1) {
  
  # Given an nmfModel object, this function plots the abundances 
  # (coefficients) corresponding to one source (basis vector) on 
  # a particular image slice
  
#  library(rasterImage)
  
  dims <- NMF::misc(nmfMod)
  nRows <- dims$nRows
  nCols <- dims$nCols
  strComp <- match(names(dims),"nSlices")
  
  if(sum(is.na(strComp)) < length(strComp)) {
    nSlices <- dims$nSlices
  }
  else {
    nSlices <- 1
  }
  
  H <- NMF::coef(nmfMod)[source,]
  H <- array(H,c(nRows,nCols,nSlices))
  
  if (length(dim(H))==3) {
  H <- H[,,slice]
  }
  
  grDevices::dev.new()
  graphics::plot(1:nRows, 1:nCols, type = "n")
  
  rasterImage::rasterImage2(z = t(H[nRows:1,]))
  
}
