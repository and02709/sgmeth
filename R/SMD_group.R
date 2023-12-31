#' This function peforms penalized matrix decomposition for one pair of
#' u and v vectors
#' Uses number of nonzero groups for the shrinkage mechanism
#' @param X matrix to undergo PMD
#' @param d supplied singular value
#' @param u supplied left vector
#' @param v supplied right vector
#' @param n number of observations
#' @param p number of predictors
#' @param groups vector of groups to which each predictor belongs
#' @param nonzero.groups number of desired nonzero groups
#' @param unique.group vector of each unique group label
#' @param num.group number of unique groups
#' @param niter number of iterations for the algorithm
#' @param trace displays information about progress of algorithm
#' @param npc number of desired principal components
#' @keywords Supervised Princpal Component Analysis
#' @export
#' @examples SMD.group(X=Xuse, d=temp$d, u=temp$u, v=temp$v, n=n, p=p, groups=groups, nonzero.groups=nonzero.groups, unique.group=unique.group, num.group=num.group,niter=niter, trace=trace, npc=i)

SMD_group <- function(X, d, u, v, n, p, groups, nonzero.groups, unique.group, 
                      num.group, niter, trace, npc){
  oldv <- rnorm(p,0,1)
  oldu <- rnorm(n,0,1)
  if(trace) cat("Vector ", npc, ": ")
  for(iter in 1:niter){
    if((sum(abs(oldv-v)) < 1e-7) && (sum(abs(oldu-u)) < 1e-7)) break
    oldv <- v
    oldu <- u
    if(trace) cat(iter," ", fill=F)
    # calculate vector u
    u <- X%*%v
    # normalize vector u
    u <- matrix(matrix(u/sgmeth::l2n(u)),ncol=1)
    # Get measurements for each group of predictors
    temp <- sgmeth::get.Xk(X=X,u=u,groups=groups,
                   unique.group=unique.group,num.group=num.group)
    # Get lambda based on the relative strength of each group
    #   we calculate the lambda that gives us the desired number of
    #   nonzero groups
    lambda <- sgmeth::calculate.lambda(Xtu.vec=temp$Xtu.vec, 
                               Xk.col.vec=temp$Xk.col.vec, 
                               nonzero.groups=nonzero.groups)
    # Shrink each vector by the approprite amount determined by lambda
    shrunk.list <- sgmeth::vector.shrinkage(Xtu.list = temp$Xtu.list, Xtu.vec=temp$Xtu.vec, 
                                    Xk.col.vec=temp$Xk.col.vec, lambda=lambda)
    # Replace shrunk values for each group back into the vector v
    v <- sgmeth::embed.v(shrunk.Xtu.list=shrunk.list,vold=v)
    # normalize vector v
    v <- matrix(matrix(v/sgmeth::l2n(v)),ncol=1)
  }
  if(trace) cat("\n")
  d <- as.numeric(t(u)%*%(X%*%v))
  return(list(d=d, u=u, v=v))
}