gen_v3<-function(pseq,J,j0,bool=TRUE){
    v=rep(0,sum(pseq))
    
    if (j0==1){
      indice=1:(sum(pseq[1:j0]))
    }
    else{
      indice=(sum(pseq[1:(j0-1)])+1):(sum(pseq[1:j0]))
    }
    
    indice=indice[which(indice%%2==0)]
    indice=indice[1:(length(indice)-1)]

    if(isTRUE(bool)){
      v[indice]=rep(1/2,length(indice))
    } else{
      v[indice]=rep(c(1/2,-1/2),length(indice)/2)
    }
    
    v[indice]=v[indice]/(sqrt(length(J))*sqrt(sum(v[indice]^2)))
    

    for (k in J){
        if (k!=j0){
          if (k==1){
            indice=1:(sum(pseq[1:k]))
          }else{
            indice=(sum(pseq[1:(k-1)])+1):(sum(pseq[1:k]))
          }
          
          
          nnz_index=sample(indice,0.5*length(indice))
          
          v[nnz_index]=runif(length(nnz_index),min=0.5,max=1)
          v[nnz_index]=v[nnz_index]/(sqrt(sum(v[nnz_index]^2))*sqrt(length(J)))     #normalize
        }
        }
     
    
    return(v)
  }
  
  
  