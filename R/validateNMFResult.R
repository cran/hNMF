# Project: hNMF_git
# 
# Author: nsauwen <nicolas.sauwen@openanalytics.eu>
#
# Description: This function converts NMF abundance maps to actual segmentation results,
# and validates the result with respect to manual segmentation. The manual segmentation results 
# are provided in nifti (.nii) format. Segmentation overlap is quantified by means of the Dice-score.
# Input to the method is the NMF input data and a valid nmfModel object obtained from the NMF package.
# sliceRange must be a 2-element vector specifying the analyzed slice range.
###############################################################################

## #' Validation tool for (h)NMF based segmentation.
## #' 
## #' This function converts NMF abundance maps to actual segmentation results,
## #' and validates the result with respect to manual segmentation. The manual segmentation results 
## #' are provided in nifti (.nii) format. Segmentation overlap is quantified by means of the Dice-score.
## #' Input to the method is the NMF input data and a valid nmfModel object obtained from the NMF package.
## #' sliceRange must be a 2-element vector specifying the analyzed slice range.
## #' 
## #' @param nmfInput List with NMF input attributes
## #' @param nmfModel NMF model (in accordance with NMF package definition)
## #' @param sliceRange 2-element vector containing first and last image slice number
## #' that were analyzed
## #' @param indTumor,indNecrosis,indEdema Indices of the NMF sources representing 
## #' tumor, necrosis and edema, respectively
## #' @param tumorNiftiFile,necrosisNiftiFile,edemaNiftiFile Character strings containing 
## #' the path to the manual segmentation files (in nifti format) for tumor, necrosis and edema
## #' @return Dice scores quantifying the segmentation quality are printed in the R console
## #' @author Nicolas Sauwen
validateNMFResult <- function(nmfInput, nmfModel, sliceRange, indTumor = 0, indNecrosis = 0, indEdema = 0, tumorNiftiFile = NULL, necrosisNiftiFile = NULL, edemaNiftiFile = NULL) {
 
#  library(oro.nifti)
#  library(tcltk)
  
  H <- NMF::coef(nmfModel)
  selectVect <- nmfInput$selectVect
  nRows <- c(nmfInput$numRows)
  nCols <- c(nmfInput$numCols)
  nSlices <- c(nmfInput$numSlices)
  bgImageTensor <- nmfInput$bgImageTensor
  first <- 2
  mid <- round(nSlices/2)
  last <- nSlices-2
  
  # Boolean variables to determine which pathological tissue types where found
  tumorFound <- 0
  necrosisFound <- 0
  edemaFound <- 0
  
  viewNMFResult(nmfInput, nmfModel)
  if(indTumor == 0) {
  	indTumor <- readline("Specify source number for active tumor:")
  }
  if(indNecrosis == 0) {
  	indNecrosis <- readline("Specify source number for necrosis:")
  }
  if(indEdema == 0) {
  	indEdema <- readline("Specify source number for edema:")
  }
  
  indTumor <- as.numeric(indTumor)
  indNecrosis <- as.numeric(indNecrosis)
  indEdema <- as.numeric(indEdema)
  
  if(indTumor != 0) {
    tumorFound <- 1
    if(is.null(tumorNiftiFile)) {
    	tumorNiftiFile <- tcltk::tk_choose.files(default = dirname(file.path(getwd(),".")), caption = "Select tumor segmentation .nii file")
    }
    tumorNiftiFile <- tumorNiftiFile[length(tumorNiftiFile)]
    tumorNifti <- oro.nifti::readNIfTI(tumorNiftiFile,reorient = FALSE)
    tumor <- oro.nifti::img_data(tumorNifti)[,,sliceRange[1]:sliceRange[2]]
    tumor <- aperm(tumor,c(2,1,3))
    tumor[selectVect==0] <- 0
  }
  
  if(indNecrosis != 0) {
    necrosisFound <- 1
    if(is.null(necrosisNiftiFile)) {
    	necrosisNiftiFile <- tcltk::tk_choose.files(default = dirname(file.path(getwd(),".")), caption = "Select necrosis segmentation .nii file")
    }
    necrosisNiftiFile <- necrosisNiftiFile[length(necrosisNiftiFile)]
    necrosisNifti <- oro.nifti::readNIfTI(necrosisNiftiFile,reorient = FALSE)
    necrosis <- oro.nifti::img_data(necrosisNifti)[,,sliceRange[1]:sliceRange[2]]
    necrosis <- aperm(necrosis,c(2,1,3))
    necrosis[selectVect==0] <- 0
  }
  
  if(indEdema != 0) {
    edemaFound <- 1
    if(is.null(edemaNiftiFile)) {
   	 edemaNiftiFile <- tcltk::tk_choose.files(default = dirname(file.path(getwd(),".")), caption = "Select edema segmentation .nii file")
  	}
    edemaNiftiFile <- edemaNiftiFile[length(edemaNiftiFile)]
    edemaNifti <- oro.nifti::readNIfTI(edemaNiftiFile,reorient = FALSE)
    edema <- oro.nifti::img_data(edemaNifti)[,,sliceRange[1]:sliceRange[2]]
    edema <- aperm(edema,c(2,1,3))
    edema[selectVect==0] <- 0
  }
  
  HMax <- apply(H,1,max)
  HScaled <- H/kronecker(matrix(1,1,ncol(H)),HMax)
  HScaledRed <- HScaled[,selectVect==1]
  
  # Conversion from NMF abundance maps to segmentation result is based 
  # on kmeans clustering:
  kmeansObj <- stats::kmeans(t(HScaledRed),diag(nrow(HScaledRed)),iter.max = 100, algorithm = "Lloyd")
  clusterIdx <- kmeansObj$cluster
  
  NMFTumor <- array(0, dim = c(nRows, nCols, nSlices))
  NMFNecrosis <- array(0, dim = c(nRows, nCols, nSlices))
  NMFEdema <- array(0, dim = c(nRows, nCols, nSlices))
  
  if(tumorFound) {
    NMFTumorRed <- rep(0,ncol(HScaledRed))
    NMFTumorRed[clusterIdx == indTumor] <- 1
    NMFTumor[selectVect == 1] <- NMFTumorRed
    if(nSlices>1) {
      makeFigure(width = 10,height = 6)
      graphics::layout(t(matrix(c(1:3))))
      imoverlay_2masks(bgImageTensor[,,first],NMFTumor[,,first],color1 = c(0,1,0),tumor[,,first],color2 = c(0,0,1))
      imoverlay_2masks(bgImageTensor[,,mid],NMFTumor[,,mid],color1 = c(0,1,0),tumor[,,mid],color2 = c(0,0,1))
      graphics::title("Tumor segmentation, green = manual, blue = NMF")
      imoverlay_2masks(bgImageTensor[,,last],NMFTumor[,,last],color1 = c(0,1,0),tumor[,,last],color2 = c(0,0,1))
    }
  }
  
  if(necrosisFound) {
    NMFNecrosisRed <- rep(0,ncol(HScaledRed))
    NMFNecrosisRed[clusterIdx == indNecrosis] <- 1
    NMFNecrosis[selectVect == 1] <- NMFNecrosisRed
    if(nSlices>1) {
      makeFigure(width = 10,height = 6)
      graphics::layout(t(matrix(c(1:3))))
      imoverlay_2masks(bgImageTensor[,,first],NMFNecrosis[,,first],color1 = c(0,1,0),necrosis[,,first],color2 = c(0,0,1))
      imoverlay_2masks(bgImageTensor[,,mid],NMFNecrosis[,,mid],color1 = c(0,1,0),necrosis[,,mid],color2 = c(0,0,1))
      graphics::title("Necrosis segmentation, green = manual, blue = NMF")
      imoverlay_2masks(bgImageTensor[,,last],NMFNecrosis[,,last],color1 = c(0,1,0),necrosis[,,last],color2 = c(0,0,1))
    }
  }
  
  if(edemaFound) {
    NMFEdemaRed <- rep(0,ncol(HScaledRed))
    NMFEdemaRed[clusterIdx == indEdema] <- 1
    NMFEdema[selectVect == 1] <- NMFEdemaRed
    if(nSlices>1) {
      makeFigure(width = 10,height = 6)
      graphics::layout(t(matrix(c(1:3))))
      imoverlay_2masks(bgImageTensor[,,first],NMFEdema[,,first],color1 = c(0,1,0),edema[,,first],color2 = c(0,0,1))
      imoverlay_2masks(bgImageTensor[,,mid],NMFEdema[,,mid],color1 = c(0,1,0),edema[,,mid],color2 = c(0,0,1))
      graphics::title("Edema segmentation, green = manual, blue = NMF")
      imoverlay_2masks(bgImageTensor[,,last],NMFEdema[,,last],color1 = c(0,1,0),edema[,,last],color2 = c(0,0,1))
    }
  }
  
  # Computation of the Dice scores
dice_tumor <- 100*2*length(which(tumor + NMFTumor == 2))/(length(which(tumor != 0)) + length(which(NMFTumor != 0)))
dice_necrosis <- 100*2*length(which(necrosis + NMFNecrosis == 2))/(length(which(necrosis != 0)) + length(which(NMFNecrosis != 0)))
dice_edema <- 100*2*length(which(edema + NMFEdema == 2))/(length(which(edema != 0)) + length(which(NMFEdema != 0)))
dice_core <- 100*2*(length(which(tumor + NMFTumor == 2)) + length(which(necrosis + NMFNecrosis == 2)) +
      length(which(tumor + NMFNecrosis == 2)) + length(which(necrosis + NMFTumor == 2)) - 
      length(which(tumor + necrosis +  NMFNecrosis == 3)) - length(which(tumor + necrosis +  NMFTumor == 3))) / 
      (length(which(tumor + necrosis != 0)) + length(which(NMFTumor + NMFNecrosis != 0)))
dice_total <- 100*2*(length(which(tumor + NMFTumor == 2)) + length(which(tumor + NMFNecrosis == 2)) +
      length(which(tumor + NMFEdema == 2)) + length(which(necrosis + NMFTumor == 2)) +
      length(which(necrosis + NMFNecrosis == 2)) + length(which(necrosis + NMFEdema == 2)) +
      length(which(edema + NMFTumor == 2)) + length(which(edema + NMFNecrosis == 2)) +
      length(which(edema + NMFEdema == 2)) - length(which(tumor + necrosis +  NMFTumor == 3)) -
      length(which(tumor + necrosis +  NMFNecrosis == 3)) - length(which(tumor + necrosis +  NMFEdema == 3)) -
      length(which(tumor + edema +  NMFTumor == 3)) - length(which(tumor + edema +  NMFNecrosis == 3)) -
      length(which(tumor + edema +  NMFEdema == 3)) - length(which(necrosis + edema +  NMFTumor == 3)) -
      length(which(necrosis + edema +  NMFNecrosis == 3)) - length(which(necrosis + edema +  NMFEdema == 3))) /
      (length(which(tumor + necrosis + edema != 0)) + length(which(NMFTumor + NMFNecrosis + NMFEdema != 0)))
print(paste("dice tumor:",as.character(round(dice_tumor, digits = 2)),"%"))
print(paste("dice necrosis:",as.character(round(dice_necrosis, digits = 2)),"%"))
print(paste("dice edema:",as.character(round(dice_edema, digits = 2)),"%"))
print(paste("dice core:",as.character(round(dice_core, digits = 2)),"%"))
print(paste("dice total:",as.character(round(dice_total, digits = 2)),"%"))  

dice_scores <- c(dice_tumor, dice_necrosis, dice_edema, dice_core, dice_total)
return(dice_scores)
  
}
