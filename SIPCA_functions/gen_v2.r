gen_v2<-function(pseq,J,j0,bool=TRUE){
  v=rep(0,sum(pseq))
  
  for (k in J){
    if (k==1){
      indice=1:(sum(pseq[1:k]))
    }
    else{
      indice=(sum(pseq[1:(k-1)])+1):(sum(pseq[1:k]))
    }
    
    if(k==j0){
      if(isTRUE(bool)){
        indice=indice[which(indice%%2==0)]
      }
      else{
        indice=indice[which(indice%%2==1)]
      }
    }
    
    nnz_index=sample(indice,0.5*length(indice))
    
    v[nnz_index]=runif(length(nnz_index),min=0.5,max=1)
    v[nnz_index]=v[nnz_index]/(sqrt(sum(v[nnz_index]^2))*sqrt(length(J))) #normalize
  }
  
  return(v)
}


