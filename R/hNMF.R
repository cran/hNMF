# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################


#' Hierarchical non-negative matrix factorization. 
#' @param nmfInput List with NMF input attributes
#' @param nmfMethod String referring to the NMF algorithm to be used.
#' @return Resulting NMF model (in accordance with NMF package definition)
#' @author Nicolas Sauwen
#' @export
hNMF <- function(nmfInput,nmfMethod='HALSacc') {
  
#  library(NMF)
#  library(nnls)
  
  nRows <- c(nmfInput$numRows)
  nCols <- c(nmfInput$numCols)
  nSlices <- c(nmfInput$numSlices)
  NMFdata <- nmfInput$data
  selectVect <- nmfInput$selectVect #This tensor codes which (non-zero) image voxels are included for analysis
  edgeROI <- nmfInput$edgeROI
  bgTensor <- nmfInput$bgImageTensor
  
  # loop below is to verify if the input data have the right input format
  if(ncol(NMFdata) == nRows*nCols*nSlices) {
    NMFdata <- NMFdata[, selectVect==1]
  }
  
  # Make sure there are no negative values in the input matrix:
  
  NMFdata[NMFdata<0] <- 0
  
  nSignals <- ncol(NMFdata)
  nFeatures <- nrow(NMFdata)
  
  # SPA initialization for NMF factor matrices
  W0 <-  initializeSPA(NMFdata,2)
  
  for(iCol in 1:ncol(W0)) {
    W0[,iCol] <- W0[,iCol]/norm(matrix(W0[,iCol]),'2')
  }
  
  H0 <- matrix(0,2,nSignals)
  
  for(iSignal in 1:nSignals) {
    nlsFit <- nnls::nnls(W0,NMFdata[,iSignal]) 
    H0[,iSignal] <- stats::coefficients(nlsFit)
  }
  
#  # Check if NMF method is already available in the NMF package
#  strComp <- match(NMF::nmfAlgorithm(),nmfMethod)
#  if(sum(is.na(strComp)) == length(strComp)) {
#    NMF::setNMFMethod(nmfMethod,get(nmfMethod))
#  }
  
  nmfInit <- NMF::nmfModel(W=W0,H=H0)
  nmfOutputLevel1 <- NMF::nmf(x=NMFdata, rank=2, method=get(nmfMethod), seed=nmfInit)
  
  W2 <- NMF::basis(nmfOutputLevel1)
  H2 <- NMF::coef(nmfOutputLevel1)
  
  for(iCol in 1:ncol(W2)) {
    W2[,iCol] <- W2[,iCol]/norm(matrix(W2[,iCol]),'2')
  }
  
  # NMF result is not unique regarding permutation. Code below takes care of consistent ordering:
  if(W2[1,1] > W2[1,2]) {
    order <-  c(1,2)
  }  else {
      order <- c(2,1)
  }
  
  # Visualizing abundance maps from hNMF level 1, to be able to 
  # determine the ranks for level 2
  H2_1 <- array(0,dim = c(nRows,nCols,nSlices))
  H2_2 <- array(0,dim = c(nRows,nCols,nSlices))
  H2_1[selectVect==1] <- H2[order[1],]
  H2_2[selectVect==1] <- H2[order[2],]
  
  if(nSlices == 1) {
    makeFigure(width=10,height=5)
    graphics::layout(t(c(1,2,3)))
    graphics::plot(1:nRows, 1:nCols, type = "n")
    imoverlay(bgTensor,edgeROI)
    graphics::plot(1:nRows, 1:nCols, type = "n")
    imoverlay(bgTensor,H2_1)
    graphics::plot(1:nRows, 1:nCols, type = "n")
    imoverlay(bgTensor,H2_2)
  }
  else {
    ind1 <- 2
    ind2 <- round(nSlices/2)
    ind3 <- nSlices-1
    makeFigure(width=10,height=10)
    graphics::layout(t(matrix(c(1,2,3,4,5,6,7,8,9),3,3)))
    imoverlay(bgTensor[,,ind1],edgeROI)
    imoverlay(bgTensor[,,ind1],H2_1[,,ind1])
    imoverlay(bgTensor[,,ind1],H2_2[,,ind1])
    imoverlay(bgTensor[,,ind2],edgeROI)
    imoverlay(bgTensor[,,ind2],H2_1[,,ind2])
    imoverlay(bgTensor[,,ind2],H2_2[,,ind2])
    imoverlay(bgTensor[,,ind3],edgeROI)
    imoverlay(bgTensor[,,ind3],H2_1[,,ind3])
    imoverlay(bgTensor[,,ind3],H2_2[,,ind3])
  }
  
  rank1 <- readline("Specify number of tissue types represented by source 1:")
  rank2 <- readline("Specify number of tissue types represented by source 2:")
  rank1 <- as.numeric(rank1)
  rank2 <- as.numeric(rank2)
  
  W <- hNMFloop(NMFdata,nRows,nCols,nSlices,W2[,order],H2[order,],selectVect,rank1,rank2,nmfMethod)
  
  Htemp <- matrix(0,rank1 + rank2,nSignals)
  
  for(iSignal in 1:nSignals) {
    nlsFit <- nnls::nnls(W,NMFdata[,iSignal]) 
    Htemp[,iSignal] <- stats::coefficients(nlsFit)
  }
  
  H <- matrix(0,rank1+rank2,nRows*nCols*nSlices)
  H[,selectVect==1] <- Htemp
  
  nmfMod <- NMF::nmfModel(W = W, H = H)
  return(nmfMod)
}


hNMFloop <- function(data,nRows,nCols,nSlices,W_old,H_old,selectVect,rank1,rank2,nmfMethod) {
  
  H_max <- apply(H_old,1,max)
  H_old <- H_old/kronecker(matrix(1,1,ncol(H_old)),H_max)
  
  # k-means clustering is used to subdivide all the data points 
  # over the 2 sources of the first NMF level 
  kmeansObj <- stats::kmeans(t(H_old),diag(nrow(H_old)),iter.max = 100, algorithm = "Lloyd")
  clusterIdx <- kmeansObj$cluster
  data1 <- data[,clusterIdx == 1]
  data2 <- data[,clusterIdx == 2]
  
  # NMF is now again applied on the 2 created datasets
  
  # SPA initialization for NMF factor matrices
  W01 <-  initializeSPA(data1,rank1)
  
  for(iCol in 1:ncol(W01)) {
    W01[,iCol] <- W01[,iCol]/norm(matrix(W01[,iCol]),'2')
  }
  
  H01 <- matrix(0,rank1,ncol(data1))
  
  for(iSignal in 1:ncol(H01)) {
    nlsFit1 <- nnls::nnls(W01,data1[,iSignal]) 
    H01[,iSignal] <- stats::coefficients(nlsFit1)
  }
  
  nmfInit1 <- NMF::nmfModel(W=W01,H=H01)
  nmfOutput1 <- NMF::nmf(x=data1,rank=rank1,method=get(nmfMethod),seed=nmfInit1)
  
  Wtemp1 <- NMF::basis(nmfOutput1)
  
  
  # SPA initialization for NMF factor matrices
  W02 <-  initializeSPA(data2,rank2)
  
  for(iCol in 1:ncol(W02)) {
    W02[,iCol] <- W02[,iCol]/norm(matrix(W02[,iCol]),'2')
  }
  
  H02 <- matrix(0,rank2,ncol(data2))
  
  for(iSignal in 1:ncol(H02)) {
    nlsFit2 <- nnls::nnls(W02,data2[,iSignal]) 
    H02[,iSignal] <- stats::coefficients(nlsFit2)
  }
  
  nmfInit2 <- NMF::nmfModel(W=W02,H=H02)
  nmfOutput2 <- NMF::nmf(x=data2,rank=rank2,method=get(nmfMethod),seed=nmfInit2)
  
  Wtemp2 <- NMF::basis(nmfOutput2)
  
  W_new <- cbind(Wtemp1,Wtemp2)
  
  for(iCol in 1:ncol(W_new)) {
    W_new[,iCol] <- W_new[,iCol]/norm(matrix(W_new[,iCol]),'2')
  }
    
  return(W_new)
}






























