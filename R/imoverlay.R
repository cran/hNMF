# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################


#' Overlay a mask or a color scaled image on top of a background image
#' @param image A matrix, background image
#' @param overlay A matrix, serving as the overlay mask or figure
#' @param selectVect A matrix (binary values), specifying which matrix elements are to be overlaid
#' @param color 3-element vector, defining the RGB color to be used in case the overlay is a mask
#' @return 
#' @author Nicolas Sauwen
#' @export
imoverlay <- function(image, overlay, selectVect = NULL, color = c(0,1,0)) {
  
  # This function takes a background image and an overlay 
  # image, and visualizes the assembled image in a new graph. If 
  # the overlay i mage is a mask, then it will be overlaid in the 
  # color specified by color.
  
  nCols <- ncol(image)
  nRows <- nrow(image)
  
  colorScale1 <- grDevices::gray(0:255/255) # Color scale for background image
  image <- round(image/max(image)*255)+1
  imageRGB <- matrix(data = colorScale1[image],nrow = nRows, ncol = nCols)
  
  # if-loop below verifies whether overlay image is a mask or an actual image
  if (length(unique(c(overlay)))>2) {
    colorScale2 <- rasterImage::colorPalette(256)
    if(is.null(selectVect)) {
      overlayNonzero <- overlay[overlay!=0]
      overlayNonzero <- round(overlayNonzero/max(overlayNonzero)*255+1)
      overlayNonzeroRGB <- colorScale2[overlayNonzero]
      imageRGB[overlay != 0] <- overlayNonzeroRGB
    }
    else {
      overlayNonzero <- overlay[selectVect!=0]
      overlayNonzero <- round(overlayNonzero/max(overlayNonzero)*255+1)
      overlayNonzeroRGB <- colorScale2[overlayNonzero]
      imageRGB[selectVect != 0] <- overlayNonzeroRGB
    }
  }
  else {
    overlayNonzero <- overlay[overlay!=0]
    overlayNonzeroRGB <- grDevices::rgb(color[1],color[2],color[3])
    imageRGB[overlay != 0] <- overlayNonzeroRGB
  }
  
#  imageRGB[overlay != 0] <- overlayNonzeroRGB
  
#  dev.new()
  graphics::plot(1:nRows, 1:nCols, type = "n", xlab = "", ylab = "") #, axes = FALSE)
  
  graphics::rasterImage(imageRGB,xleft = 1,ybottom = 1,xright = nCols,ytop = nRows)
  
}
