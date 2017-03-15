# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################


imoverlay_2masks <- function(image, overlay1, color1 = c(0,1,0), overlay2, color2 = c(1,0,0)) {
  
  # This function takes a background image and an overlay 
  # image, and visualizes the assembled image in a new graph. If 
  # the overlay i mage is a mask, then it will be overlaid in the 
  # color specified by color.
  
  nCols <- ncol(image)
  nRows <- nrow(image)
  
  colorScale1 <- grDevices::gray(0:255/255) # Color scale for background image
  image <- round(image/max(image)*255)+1
  imageRGB <- matrix(data = colorScale1[image],nrow = nRows, ncol = nCols)
  
  overlay2[overlay2 != 0] <- 2
  overlay <- overlay1 + overlay2
  
  overlayNonzeroRGB1 <- grDevices::rgb(color1[1],color1[2],color1[3])
  imageRGB[overlay == 1] <- overlayNonzeroRGB1
  
  overlayNonzeroRGB2 <- grDevices::rgb(color2[1],color2[2],color2[3])
  imageRGB[overlay == 2] <- overlayNonzeroRGB2
  
  overlayNonzeroRGB3 <- grDevices::rgb(color1[1]+color2[1],color1[2]+color2[2],color1[3]+color2[3])
  imageRGB[overlay == 3] <- overlayNonzeroRGB3
  
#  imageRGB[overlay != 0] <- overlayNonzeroRGB
  
#  dev.new()
  graphics::plot(1:nRows, 1:nCols, type = "n", xlab = "", ylab = "") #, axes = FALSE)
  
  graphics::rasterImage(imageRGB,xleft = 1,ybottom = 1,xright = nCols,ytop = nRows)
  
}
