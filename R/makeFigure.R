# Project: hNMF_git
# 
# Author: nsauwen <nicolas.sauwen@openanalytics.eu>
#
# Description:
#
###############################################################################


makeFigure <- function(width, height) {
  
  errX11 <- try(grDevices::x11(width = width, height = height),TRUE)
  if(!is.null(errX11)){
    grDevices::dev.new(width = width, height = height)
  }   
}
