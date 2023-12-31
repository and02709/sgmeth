% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.gSIMPLS.R
\name{cv.gSIMPLS}
\alias{cv.gSIMPLS}
\title{cross validation function for groupwise SIMPLS}
\usage{
cv.gSIMPLS(
  X,
  Y,
  npc,
  n.folds,
  groups = NULL,
  nonzero.groups = NULL,
  kernel = c("linear", "delta"),
  parallel = F,
  n.cores = NULL,
  niter = 50,
  trace = F,
  part.balance = T,
  mc.method = T,
  ind.names = NULL
)
}
\arguments{
\item{X}{nxp predictor matrix}

\item{Y}{nx1 response vector}

\item{npc}{number of desired principal components}

\item{n.folds}{number of folds to perform cross validation}

\item{groups}{vector of groups to which each predictor belongs}

\item{nonzero.groups}{number of desired nonzero groups}

\item{kernel}{specification of the response kernel}

\item{parallel}{flag for parallel process to parApply}

\item{n.cores}{number of cores to be used in parallel process}

\item{niter}{number of iterations in the SMD algorithm}

\item{trace}{display algorithm progress}

\item{part.balance}{flag for whether folds process should balance factors}

\item{mc.method}{flag for whether the parallel method should use mclapply}

\item{ind.names}{vector of each observation label}
}
\description{
cross validation function for groupwise SIMPLS
}
\examples{
cv.gSIMPLS(X, Y ,npc, n.folds, groups=NULL, nonzero.groups=NULL, kernel=c("linear", "delta"), parallel=F, n.cores=NULL, niter=50, trace=F, part.balance=T,ind.names=NULL)
}
\keyword{Analysis}
\keyword{Component}
\keyword{Princpal}
\keyword{Supervised}
