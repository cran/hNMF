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
#' @importFrom grDevices rgb
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
  alpha <- 0.5 # Opacity of overlay
  if (length(unique(c(overlay)))>2) {
    colorScale2 <- rasterImage::colorPalette(256)
    if(is.null(selectVect)) {
      overlayNonzero <- overlay[overlay!=0]
      overlayNonzero <- round(overlayNonzero/max(overlayNonzero)*255+1)
      overlayNonzeroRGB <- colorScale2[overlayNonzero]
      imageRGB_red <- strtoi(x = substr(imageRGB,2,3), base = 16)
      imageRGB_green <- strtoi(x = substr(imageRGB,4,5), base = 16)
      imageRGB_blue <- strtoi(x = substr(imageRGB,6,7), base = 16)
      imageRGB_red <- imageRGB_red[selectVect != 0]
      imageRGB_green <- imageRGB_green[selectVect != 0]
      imageRGB_blue <- imageRGB_blue[selectVect != 0]
      overlayNonzero_red <- strtoi(x = substr(overlayNonzeroRGB,2,3), base = 16)
      overlayNonzero_green <- strtoi(x = substr(overlayNonzeroRGB,4,5), base = 16)
      overlayNonzero_blue <- strtoi(x = substr(overlayNonzeroRGB,6,7), base = 16)
      imageRGB_combined <- rgb(alpha * imageRGB_red + (1 - alpha) * overlayNonzero_red,
          alpha * imageRGB_green + (1 - alpha) * overlayNonzero_green,
          alpha * imageRGB_blue + (1 - alpha) * overlayNonzero_blue, alpha = 255, maxColorValue=255)
      imageRGB[selectVect != 0] <-  imageRGB_combined
    }
    else {
      overlayNonzero <- overlay[selectVect!=0]
      overlayNonzero <- round(overlayNonzero/max(overlayNonzero)*255+1)
      overlayNonzeroRGB <- colorScale2[overlayNonzero]
      imageRGB_red <- strtoi(x = substr(imageRGB,2,3), base = 16)
      imageRGB_green <- strtoi(x = substr(imageRGB,4,5), base = 16)
      imageRGB_blue <- strtoi(x = substr(imageRGB,6,7), base = 16)
      imageRGB_red <- imageRGB_red[selectVect != 0]
      imageRGB_green <- imageRGB_green[selectVect != 0]
      imageRGB_blue <- imageRGB_blue[selectVect != 0]
      overlayNonzero_red <- strtoi(x = substr(overlayNonzeroRGB,2,3), base = 16)
      overlayNonzero_green <- strtoi(x = substr(overlayNonzeroRGB,4,5), base = 16)
      overlayNonzero_blue <- strtoi(x = substr(overlayNonzeroRGB,6,7), base = 16)
      imageRGB_combined <- rgb(alpha * imageRGB_red + (1 - alpha) * overlayNonzero_red,
        alpha * imageRGB_green + (1 - alpha) * overlayNonzero_green,
        alpha * imageRGB_blue + (1 - alpha) * overlayNonzero_blue, alpha = 255, maxColorValue=255)
      imageRGB[selectVect != 0] <-  imageRGB_combined
    }
  }
  else {
    overlayNonzero <- overlay[overlay!=0]
    overlayNonzeroRGB <- grDevices::rgb(color[1],color[2],color[3])
    imageRGB[overlay != 0] <- overlayNonzeroRGB
  }
  
  graphics::plot(1:nRows, 1:nCols, type = "n", xlab = "", ylab = "", legend = T, axes = FALSE)  
  graphics::rasterImage(imageRGB,xleft = 1,ybottom = 1,xright = nCols,ytop = nRows)
  
#  tiff("Tissue abundance.tiff", type ="cairo")
#  graphics::plot(1:nRows+50, 1:nCols, type = "n", xlab = "", ylab = "", legend = T, axes = FALSE)  
#  graphics::rasterImage(imageRGB,xleft = 18,ybottom = 1,xright = nCols+19,ytop = nRows)
#  color.legend(210,20,215,195,legend=c("0","0.4", "0.8", "1.2"), cex = 0.9, col = "white", rect.col=colorScale2, align = "rb", gradient = "y")
#  text(x = 233, y = 108, labels = "Tissue abundance", col = "white", cex=1, srt = 90)
#  dev.off()
  
}
