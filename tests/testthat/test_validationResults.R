# Project: hNMF_git
# 
# Author: nsauwen <nicolas.sauwen@openanalytics.eu>
#
# Description:
#
###############################################################################


library(hNMF)
context("Computation of Dice-scores from segmentation result")

setwd(system.file("extdata",package="hNMF"))

#tumorNiftiFile <- system.file("extdata","tumor.nii",package="hNMF")
tumorNifti <- oro.nifti::readNIfTI(file.path(system.file('extdata', package = 'hNMF'), "tumor.nii"),reorient = FALSE)
tumor <- oro.nifti::img_data(tumorNifti)
tumor <- aperm(tumor,c(2,1,3))


#necrosisNiftiFile <- system.file("extdata","necrosis.nii",package="hNMF")
necrosisNifti <- oro.nifti::readNIfTI(file.path(system.file('extdata', package = 'hNMF'), "necrosis.nii"),reorient = FALSE)
necrosis <- oro.nifti::img_data(necrosisNifti)
necrosis <- aperm(necrosis,c(2,1,3))

#edemaNiftiFile <- system.file("extdata","edema.nii",package="hNMF")
edemaNifti <- oro.nifti::readNIfTI(file.path(system.file('extdata', package = 'hNMF'), "edema.nii"),reorient = FALSE)
edema <- oro.nifti::img_data(edemaNifti)
edema <- aperm(edema,c(2,1,3))

nRows <- dim(tumor)[1]
nCols <- dim(tumor)[2]
nSlices <- dim(tumor)[3]
bgImageTensor <- array(0,dim=dim(tumor))
selectVect <- tumor + necrosis + edema
sliceRange <- c(1,nSlices)

tumor[selectVect==0] <- 0
necrosis[selectVect==0] <- 0
edema[selectVect==0] <- 0

H <- matrix(0,nrow=3,ncol=nCols*nRows*nSlices)
H[1,] <- tumor
H[2,] <- necrosis
H[3,] <- edema

nmfInput <- NULL
nmfInput$numRows <- nRows
nmfInput$numCols <- nCols
nmfInput$numSlices <- nSlices
nmfInput$bgImageTensor <- bgImageTensor
nmfInput$selectVect <- selectVect

nmfMod <- NMF::nmfModel(H=H)

test_that("Dice-scores are 100% when segmentation masks are identical", {
      dice_scores <- validateNMFResult(nmfInput, nmfMod, sliceRange, indTumor = 1, indNecrosis = 2, indEdema = 3, tumorNiftiFile = "tumor.nii", necrosisNiftiFile = "necrosis.nii", edemaNiftiFile = "edema.nii")
      expect_equal(dice_scores, c(100,100,100,100,100))
    })
