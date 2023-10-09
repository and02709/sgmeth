#' This function deflates the matrix X based on components V
#' @param X predictor matrix
#' @param V component vectors
#' @param n number of observations
#' @param p number of predictors
#' @keywords Sparse partial least squares
#' @export
#' @examples deflate(X=X,V=V,n=200,p=1000)

deflate <- function(X,V,n,p){
  # Encode vectors to data matrix
  Z <- Rfast::mat.mult(X,as.matrix(V))
  
  # Generate spanning vector
  #P <- t(X)%*%Z%*%solve(t(Z)%*%Z)
  P <- Rfast::mat.mult(Rfast::Crossprod(X,Z),Rfast::spdinv(Rfast::Crossprod(Z,Z)))
  
  # Generate orthogonal projection matrix
  #Q <- diag(x=1,nrow = p) - P%*%solve(t(P)%*%P)%*%t(P)
  Q <- diag(x=1,nrow = p) - Rfast::mat.mult(P,Rfast::Tcrossprod(Rfast::spdinv(Rfast::Crossprod(P,P)),P))
  
  return(Rfast::mat.mult(X,Q))
}