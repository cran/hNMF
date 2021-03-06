% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PGNMF.R
\name{PGNMF}
\alias{PGNMF}
\title{NMF by alternating non-negative least squares using projected gradients. 
For a reference to the method, see C.-J. Lin, 
"Projected Gradient Methods for Non-negative Matrix Factorization", 
Neural computation 19.10 (2007): 2756-2779.}
\usage{
PGNMF(X, nmfMod, tol = 1e-05, maxIter = 500, timeLimit = 300,
  checkDivergence = TRUE)
}
\arguments{
\item{X}{Input data matrix, each column represents one data point 
and the rows correspond to the different features}

\item{nmfMod}{Valid NMF model, containing initialized factor matrices
(in accordance with the NMF package definition)}

\item{tol}{Tolerance for a relative stopping condition}

\item{maxIter}{Maximum number of iterations}

\item{timeLimit}{Limit of time duration NMF analysis}

\item{checkDivergence}{Boolean indicating whether divergence checking should be performed
Default is TRUE, but it should be set to FALSE when using random initialization}
}
\value{
Resulting NMF model (in accordance with the NMF package definition)
}
\description{
NMF by alternating non-negative least squares using projected gradients. 
For a reference to the method, see C.-J. Lin, 
"Projected Gradient Methods for Non-negative Matrix Factorization", 
Neural computation 19.10 (2007): 2756-2779.
}
\author{
nsauwen
}
