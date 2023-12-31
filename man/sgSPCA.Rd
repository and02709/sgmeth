% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgSPCA.R
\name{sgSPCA}
\alias{sgSPCA}
\title{This function is the parent function for performing sparse groupwise supervised PCA}
\usage{
sgSPCA(
  X,
  Y,
  npc,
  kernel = c("linear", "delta"),
  groups = NULL,
  nonzero.groups = NULL,
  ind.names = NULL,
  alpha = 0,
  niter = 50,
  trace = F
)
}
\arguments{
\item{X}{nxp matrix of predictors}

\item{Y}{nx1 response vector}

\item{npc}{number of desired principal components}

\item{kernel}{designation for either linear or delta kernel}

\item{groups}{vector of groups to which each predictor belongs}

\item{nonzero.groups}{number of desired nonzero groups}

\item{ind.names}{vector of each observation label}

\item{alpha}{l1 penalty term}

\item{niter}{number of iterations for the algorithm}

\item{trace}{displays information about progress of algorithm}
}
\description{
This function is the parent function for performing sparse groupwise supervised PCA
}
\examples{
sgSPCA(X=X,Y=Y,npc=2, kernel="linear", groups=group.list, nonzero.groups=5, ind.names=NULL, alpha=0.5, niter=200, trace=F)
}
\keyword{Analysis}
\keyword{Component}
\keyword{Princpal}
\keyword{Supervised}
