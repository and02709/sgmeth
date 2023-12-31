% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.partition.group.simpls.R
\name{cv.partition.group.simpls}
\alias{cv.partition.group.simpls}
\title{cross validation function for groupwise SIMPLS}
\usage{
cv.partition.group.simpls(
  arg.group,
  df.partition,
  npc,
  n.folds,
  groups,
  kernel = c("linear", "delta"),
  ind.names,
  niter,
  trace
)
}
\arguments{
\item{arg.group}{vector for folds and sparse parameter}

\item{df.partition}{list of folds}

\item{npc}{number of desired principal components}

\item{n.folds}{number of folds to perform cross validation}

\item{groups}{vector of groups to which each predictor belongs}

\item{kernel}{specification of the response kernel}

\item{ind.names}{vector of each observation label}

\item{niter}{number of iterations in the SMD algorithm}

\item{trace}{display algorithm progress}
}
\description{
cross validation function for groupwise SIMPLS
}
\examples{
cv.partition.group.simpls(arg.group, df.partition, npc, n.folds, groups, kernel=c("linear", "delta"), ind.names, niter,trace)
}
\keyword{Analysis}
\keyword{Component}
\keyword{Princpal}
\keyword{Supervised}
