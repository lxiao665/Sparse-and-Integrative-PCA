soda <-
function(x){
## ------------------------------------------------------------------------
## Sparse Orthogonal Decomposition(SOD) Function of a matrix
## ------------------------------------------------------------------------
##
## DESCRIPTION:
##          Sparse Orthogonal Decomposition Algorithm (SODA) of the prespecified matrix
##
##
## INPUT ARGUMENTS:
## x            The matrix which SOD will be applied on
##
##
## OUTPUT ARGUMENTS:
## x            The SOD of the original matrix
##
##
##
##
##
## AUTHORS:
##
##
## REFERENCE:
##
##
##
##
##
##
##
## DEVELOPMENT:
##
##


  if(is.vector(x)==TRUE) x <- matrix(x,length(x),1);
  std <- function(x){ 
    if(sum(x^2)>0)
    return(x/sqrt(sum(x^2)));
    return(x);
    }
  
  Dim <- dim(x);
  p <- Dim[1];
  rank <- Dim[2];
  d <- rep(NA,rank);
  ## identify non-zero rows
  index <- apply(x^2,1,sum);
  index <- (1:p)[index!=0];
  
  x <- x[index,]
  x <- matrix(x,length(index),length(x)/length(index))
  d[1] <- sqrt(sum(x[,1]^2))
  y <- x ## temporary holder
  x[,1] = std(x[,1])
  if(rank>1){
  for(r in 2:rank){
    vec <- x[,r]
    ind <- (1:length(vec))[vec!=0]
    vec.ind <- vec[ind]
    for(t in 1:(r-1))
    { vec.iter <- x[ind,t]
      if(sum(vec.iter^2)>0)
      vec.ind <- vec.ind - sum(vec.ind*vec.iter)/sum(vec.iter^2)*vec.iter
    }
    
    if(sum(vec.ind^2)==0){
      ind <- 1:length(vec)
      vec.ind <- vec[ind]
      for(t in 1:(r-1))
      { vec.iter <- x[ind,t]
        if(sum(vec.iter^2)>0)
          vec.ind <- vec.ind - sum(vec.ind*vec.iter)/sum(vec.iter^2)*vec.iter
      }
      
    }
    d[r] <- sqrt(sum(vec.ind^2))
    x[ind,r] <- std(vec.ind)
  }
  }
Q <- x
R <- t(Q)%*%y

return(list(Q=Q,R = R,index = index,d=d))
}
