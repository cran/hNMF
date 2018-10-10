# Project: hNMF_git
# 
# Author: nsauwen
###############################################################################

#' Perform Non-Negative Matrix factorization
#' 
#' @param X input matrix. Each column represents one observation 
#' and the rows correspond to the different features
#' @param rank number of NMF components to be found
#' @param initData either of the NMF factor matrices, with initial values 
#' @param method name of the NMF method to be used. "PGNMF" (default) and "HALSacc" 
#' are available by default. Any method from the NMF package can also be specified
#' @param nruns number of NMF runs. It is recommended to run the NMF analyses multiple
#' times when random seeding is used, to avoid a suboptimal solution
#' @param checkDivergence Boolean indicating whether divergence checking should be performed
#' @return Scaled NMF model (in accordance with the NMF package definition)
#' @importFrom stats runif
#' @author Nicolas Sauwen
#' @export
#' @examples
#' 
#' # random data
#' X <- matrix(runif(10*20), 10,20)
#' 
#' # run NMF with default algorithm, 5 runs with random initialization
#' NMFresult1 <- oneLevelNMF(X, rank=2, nruns=5)
#' 
#' # run NMF with specified algorithm and with initialized sources
#' W0 <- initializeSPA(X,3)
#' NMFresult2 <- oneLevelNMF(X, rank=3, method="HALSacc", initData = W0)

oneLevelNMF <- function(X, rank, initData = NULL, method = "PGNMF", nruns = 10, checkDivergence = TRUE){
	
	X <- preProcesInputData(X)
	seed <- initializeNMF(X, initData)
	
	if( is.null( seed ) ) {
#    NMFResult               <- nmf( X , rank = rank, method = get(method), nrun = nruns, checkDivergence = F, .options="-p")
		
		NMFResult <- NULL
		residu <- Inf
		
		for(iRun in 1:nruns) {	
			W0 <- matrix(runif(rank*nrow(X)), nrow = nrow(X), ncol = rank)	
			H0 <- matrix(runif(rank*ncol(X)), nrow = rank, ncol = ncol(X))	
			seed <- NMF::nmfModel( W = W0 , H = H0 )
			if(method == "PGNMF"){
				NMFResult_temp <- PGNMF(X, nmfMod = seed, checkDivergence = F)
			}
			else if(method == "HALSacc"){
				NMFResult_temp <- HALSacc(X, nmfMod = seed, checkDivergence = F)
			}
			W_temp <- NMF::basis(NMFResult_temp)
			H_temp <- NMF::coef(NMFResult_temp)
			residu_temp <- norm(X - W_temp%*%H_temp,'f')
			if(residu_temp < residu) {
				NMFResult <- NMFResult_temp
			}
		}
	}
	else {
		# Check if NMF initialization is consistent with specified rank
		W0                      <- NMF::basis(seed)
		if(ncol(W0) != rank) {
			stop("Number of provided pure component X does not correspond to specified NMF rank.")    
		}
		if(is.null(checkDivergence)) {
#      NMFResult             <-  nmf( X , rank = rank , method = get(method) , nrun = 1 , seed = seed , checkDivergence = T, .options = "-p" )
			if(method == "PGNMF"){
				NMFResult <- PGNMF(X, nmfMod = seed, checkDivergence = T)
			}
			else if(method == "HALSacc"){
				NMFResult <- HALSacc(X, nmfMod = seed, checkDivergence = T)
			}
		}
		else{
#      NMFResult             <-  nmf( X , rank = rank , method = get(method) , nrun = 1 , seed = seed , checkDivergence = checkDivergence, .options = "-p" )
			if(method == "PGNMF"){
				NMFResult <- PGNMF(X, nmfMod = seed, checkDivergence = checkDivergence)
			}
			else if(method == "HALSacc"){
				NMFResult <- HALSacc(X, nmfMod = seed, checkDivergence = checkDivergence)
			}
		}
	}
	
	
	NMFResult <- scaleNMFResult( NMFResult )
	return(NMFResult)
	
}


#' Condition input data matrix properly for NMF
#' 
#' @param X input matrix
#' @return matrix with non-zero elements
#' @export
preProcesInputData <- function(X){
	
	zeroInds <- which(X < 0)
	if(length(zeroInds) > 0){
		X[ X < 0 ]    <-  0
		warning("Input data contain negative values. These will be reset to zero for NMF analysis")
	}
	
	zeroRow                   <-  which( ( apply( X , 1 , sum ) ) == 0 )
	X[ zeroRow , ]     <- 1e-16
	X[ X == 0 ] <- 1e-16
	return(X)
	
}


#' Initialize NMF model with initial spectral data
#' 
#' @param X input matrix
#' @param initData source or abundance matrix with initial values
#' @importFrom nnls nnls
#' @importFrom stats coefficients approx
#' @importFrom NMF nmfModel
#' @export
initializeNMF     <- function(X, initData = NULL) {
	
	if(is.null(initData)) {
		NMFInit            <- NULL
		return(NMFInit)
	}
	else if(is.matrix( initData ) & nrow(initData) == nrow(X)) {
		W0                 <-  initData
		rank               <-  ncol( W0 )
	}
	else if(is.matrix( initData ) & ncol(initData) == ncol(X)) {
		H0                 <-  initData
		rank               <-  nrow( H0 )
	}
	else {
		warning("NMF initialization is not provided in the right format. NMF will be ran with random initialization instead")
		NMFInit            <- NULL
		return(NMFInit)
	}
	
	if(!exists("H0")){
		nSignals           <-  ncol(X)
		H0                 <-  matrix( 0 , rank , nSignals)
		for( iSignal in 1:nSignals ) {
			nlsFit           <-  nnls::nnls( W0 , X[,iSignal] )   
			H0[,iSignal]     <-  coefficients(nlsFit)
		}
	}
	
	if(!exists("W0")){
		nFeatures           <-  nrow(X)
		W0                  <-  matrix( 0 , rank, nFeatures)
		for( iFeature in 1:nFeatures ) {
			nlsFit           <-  nnls::nnls( t(H0) , X[iFeature, ] )   
			W0[,iFeature]     <-  coefficients(nlsFit)
		}
		W0 <- t(W0)
	}

	NMFInit            <-  NMF::nmfModel( W = W0 , H = H0 )
	NMFInit            <-  scaleNMFResult(NMFInit)
	
	return(NMFInit)
}

#' Apply fixed scaling to NMF model matrices by normalizing the basis vectors
#' 
#' @param NMFResult Fitted NMF model
#' @return NMFResult Rescaled NMF model
#' @importFrom NMF basis coef
#' @author Nicolas Sauwen
#' @export
scaleNMFResult      <- function(NMFResult) {
	
	W                 <- NMF::basis(NMFResult)
	N                 <- sqrt(diag(t(W)%*%W))
	N_W               <- matrix(N,nrow = nrow(W),ncol = ncol(W),byrow = T)
	W                 <- W/N_W
	H                 <- NMF::coef(NMFResult)
	N_H               <- matrix(N,nrow = nrow(H),ncol = ncol(H))
	H                 <- H*N_H
	NMF::basis(NMFResult)  <- W
	NMF::coef(NMFResult)   <- H
	return(NMFResult)
	
}