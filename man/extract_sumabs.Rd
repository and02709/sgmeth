% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_sumabs.R
\name{extract_sumabs}
\alias{extract_sumabs}
\title{This extracts and returns the left and right singular vectors}
\usage{
extract_sumabs(X, n, p, sumabsv, niter, trace, npc)
}
\arguments{
\item{X}{predictor matrix}

\item{n}{number of observations}

\item{p}{number of predictors}

\item{niter}{number of iterations}

\item{trace}{whether to output progress}

\item{npc}{which pc is being extracted}

\item{sumabs}{value of sum of absolute values for loadings desired}
}
\description{
This extracts and returns the left and right singular vectors
}
\examples{
extract_sumabs(X=X,n=200,p=1000,sumabsv=sqrt(10),niter=200,trace=F,npc=2)
}
\keyword{Sparse}
\keyword{elast}
\keyword{partial}
\keyword{squares}
