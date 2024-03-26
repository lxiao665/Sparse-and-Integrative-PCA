gen_weight<-function(pseq){
  I=length(pseq)
  w=matrix(0,I,I)
  for(i in 1:I){
    w[i,i]=pseq[i]
    if(i<I){
      for(j in (i+1):I){
        w[i,j]=sqrt(pseq[i]*pseq[j])
        w[j,i]=w[i,j]
      }
    }
  }
  return(w)
}


BlockProj<-function(W,H1,w,rind,cind,rho,rho2){
  
  #w: w_kl
  M_kl=(1/rho)*W[rind,cind]+H1[rind,cind]
  cons=max(0,1-rho2*w/(rho*(norm(M_kl,"F"))))
  
  return(cons*M_kl)
  
}


BlockProx<-function(W,H1,w,rind,cind,rho,rho1,rho2){
  
  M_kl=(1/rho)*W[rind,cind]+H1[rind,cind]
  
  
  #first step:elementwise soft thresholding
  A_kl=matrix(unlist(sapply(c(M_kl),FUN=function(u){
    SoftThreshold(u,rho1/rho)
  })),length(rind),length(cind))
  
  #groupwise thresholding
  if(norm(A_kl,'F')==0){
    cons=0
  }
  else{
    cons=max(0,1-rho2*w/(rho*norm(A_kl,'F')))
  }
  
  return(cons*A_kl)
  
}


simplex_sum<-function(r,theta){
  vec=sapply(r,FUN=function(u){
    min(max(u-theta,0),1)
  })
  return (sum(vec))
}


#find the right endpoint of the interval which contains theta
find_first<-function(r,knots,f,ndim){
  nk=length(knots)
  for (i in 1:nk){
    if(f(r,knots[i])<ndim){
      return(i)
      break
    }
  }
}


#using interpolation to solve the piecewise linear equation
GetTheta<-function(r,ndim){
  #concatenation vector r-1 and r
  knots=as.numeric(as.character((as.data.frame(table(c(r-1,r)))[,1])))
  index=find_first(r,knots,simplex_sum,ndim) #return the index of the right endpoint 
  ia=knots[index]   #right endpoint
  ib=knots[index-1]
  fa=simplex_sum(r,ia)
  fb=simplex_sum(r,ib)
  theta=ia+(ib-ia)*(ndim-fa)/(fb-fa)  #interpolation
  return(theta)
  
}


starnorm<-function(H,pseq,w){
  #w: weight matrix
  I=length(pseq)
  Norm=0
  for(k in 1:I){
    for(l in 1:I){
      if(k==1){
        rind=1:pseq[1]
      }
      else{
        rind=(sum(pseq[1:(k-1)])+1):(sum(pseq[1:k]))   #rows corresponding to block(k,l)
      }
      if(l==1){
        cind=1:pseq[1]
      }
      else{
        cind=(sum(pseq[1:(l-1)])+1):(sum(pseq[1:l]))   #columns corresponding to block (k,l)
      }
      Norm=Norm+w[k,l]*norm(as.matrix(H[rind,cind]),'F')
      
    }
  }
  
  return(Norm)
  
}



SoftThreshold<-function(x,lambda){
  value=sign(x)*max(abs(x)-lambda,0)
  return(value)
}


Kindice<-function(n,fold){
  ### N:number of observations
  set.seed(1)
  ###folds:number of desired splits
  ##generate indices of holdout observations
  indice=split(sample(1:n),1:fold)
  return (indice)
  
  
}


GenerateRho0<-function(Sj,n_rho0,q=0.95,minv=0,maxv=NULL,skew=NULL){
  #nsol: number of grids
  p=dim(Sj)[1]
  
  diag(Sj)=rep(NA,p)
  if(is.null(minv)){
    minv=0
  }
  
  if(is.null(maxv)){
    maxv= as.numeric(quantile(as.numeric(na.omit(abs(c(Sj)))),q))
  }
  
  #generate mixing parameter sequence
  step=1/(n_rho0-1)
  if(is.null(skew)){
    skew=1
  }
  
  mixseq=unlist(sapply(seq(from=0,to=1,by=step),FUN=function(u){
    (1-exp(skew*u))/(1-exp(skew))
  }))
  
  rho0seq=(1-mixseq)*minv+mixseq*maxv
  
  return(rho0seq)
}

gen_v_unif<-function(pseq,J,alpha,i){
  v=rep(0,sum(pseq))
  
  for (k in J){
    if (k==1){
      indice=1:(sum(pseq[1:k]))
    }
    else{
      indice=(sum(pseq[1:(k-1)])+1):(sum(pseq[1:k]))
    }
    set.seed(i)
    nnz_index=sample(indice,alpha*length(indice))
    v[nnz_index]=runif(length(nnz_index),min=0,max=1)
    #v[nnz_index]=rnorm(length(nnz_index),0,1)
    v[nnz_index]=v[nnz_index]/(sqrt(sum(v[nnz_index]^2))*sqrt(length(J))) #normalize
  }
  
  return(v)
}


gen_v_unif2<-function(pseq,J,alpha,i){
  v=rep(0,sum(pseq))
  
  for (k in J){
    if (k==1){
      indice=1:(sum(pseq[1:k]))
    }
    else{
      indice=(sum(pseq[1:(k-1)])+1):(sum(pseq[1:k]))
    }
    set.seed(i)
    nnz_index=sample(indice,alpha*length(indice))
    v[nnz_index]=runif(length(nnz_index),min=0.5,max=1)
    #v[nnz_index]=rnorm(length(nnz_index),0,1)
    v[nnz_index]=v[nnz_index]/(sqrt(sum(v[nnz_index]^2))*sqrt(length(J))) #normalize
  }
  
  return(v)
}


gen_v_normal<-function(pseq,J,alpha,i){
  v=rep(0,sum(pseq))
  
  for (k in J){
    if (k==1){
      indice=1:(sum(pseq[1:k]))
    }
    else{
      indice=(sum(pseq[1:(k-1)])+1):(sum(pseq[1:k]))
    }
    set.seed((i+k))
    nnz_index=sample(indice,alpha*length(indice))
    #v[nnz_index]=runif(length(nnz_index),min=0.5,max=1)
    v[nnz_index]=rnorm(length(nnz_index),0,1)
    # v[nnz_index]=v[nnz_index]/(sqrt(sum(v[nnz_index]^2))*sqrt(length(J))) #normalize
  }
  
  v = v/(sqrt(sum(v^2)))
  return(v)
}


#only group-wise sparsity
gen_vg<-function(pseq,J){
  v=rep(0,sum(pseq))
  
  for (k in J){
    if (k==1){
      indice=1:(sum(pseq[1:k]))
    }
    else{
      indice=(sum(pseq[1:(k-1)])+1):(sum(pseq[1:k]))
    }
    nnz_index=sample(indice,length(indice))
    v[nnz_index]=runif(length(nnz_index),min=0.5,max=1)
    v[nnz_index]=v[nnz_index]/(sqrt(sum(v[nnz_index]^2))*sqrt(length(J))) #normalize
  }
  
  return(v)
}




block_norm<-function(M,pseq){
  I=length(pseq)
  Norm=matrix(0,I,I)
  
  for(k in 1:I){
    for(l in 1:I){
      if(k==1){
        rind=1:pseq[1]
      }else{
        rind=(sum(pseq[1:(k-1)])+1):(sum(pseq[1:k]))   #rows corresponding to block(k,l)
      }
      if(l==1){
        cind=1:pseq[1]
      }else{
        cind=(sum(pseq[1:(l-1)])+1):(sum(pseq[1:l]))   #columns corresponding to block (k,l)
      }
      
      Norm[k,l]=norm(as.matrix(M[rind,cind]),'F')
      
    }
  }
  return(Norm)
}



v_fac<-function(pseq){
  v=rep(1,pseq[1])
  I=length(pseq)
  for (i in 2:I){
    v=c(v,rep(i,pseq[i]))  
  }
  
  return(v)
  
}


eig_bnorm<-function(eig,pseq){
  I=length(pseq)
  fac=v_fac(pseq) 
  eig_M=cbind(eig,fac)
  
  eig_blnorm<-rep(0,I)
  
  
  for(i in 1:I){
    eig_sub=eig_M[which(eig_M[,2]==i),1]
    
    eig_blnorm[i]=norm(as.matrix(eig_sub),'2')
    
  }
  return(eig_blnorm)
}

l2norm<-function(eig,v){
  
  
  if(t(Re(eig))%*%v>0)
  {  Norm=norm(as.matrix(Re(eig)-v),'2')}
  
  else{
    Norm=norm(as.matrix(Re(eig)+v),'2')
  }
  
  return(Norm)
  
}


#false positive
##arrange eig into matrix form, then calculate the F norm of each block
FP<-function(eig,v,tol){
  
  fp=0
  index=which(v==0)
  for(i in index){
    if(abs(Re(eig[i]))>tol){
      fp=fp+1
    }
    
  }
  return(fp)
  
}

#false negative
FN<-function(eig,v,tol){
  
  fn=0
  index=which(v!=0)
  
  for(i in index){
    if(abs(Re(eig[i]))<tol){
      fn=fn+1
    }
  }
  
  return(fn)
  
}


v_fac<-function(pseq){
  v=rep(1,pseq[1])
  I=length(pseq)
  for (i in 2:I){
    v=c(v,rep(i,pseq[i]))  
  }
  
  return(v)
  
}

FP_block<-function(eig,v,pseq,tol){
  I=length(pseq)
  fac=v_fac(pseq) 
  eig_M=cbind(eig,fac)
  v_M=cbind(v,fac)
  eig_blnorm<-rep(0,I)
  v_blnorm<-rep(0,I)
  
  for(i in 1:I){
    eig_sub=eig_M[which(eig_M[,2]==i),1]
    v_sub=v_M[which(v_M[,2]==i),1]
    eig_blnorm[i]=norm(as.matrix(eig_sub),'2')
    v_blnorm[i]=norm(as.matrix(v_sub),'2')  
  }
  
  
  fp=0
  
  index=which(v_blnorm==0)
  
  for(i in index){
    if(abs(eig_blnorm[i])>tol){
      fp=fp+1
    }
    
  }
  return(fp)
  
}

#create the new eigen estimate after thresholding
##calculate block norm first, decide the index of the blocks which needed to be included
post_est<-function(eig,pseq,tol){
  
  I=length(pseq)
  fac=v_fac(pseq) 
  eig_M=cbind(eig,fac)
  eig_blnorm<-rep(0,I)
 
  for(i in 1:I){
   
    eig_sub=eig_M[which(eig_M[,2]==i),1]
    eig_blnorm[i]=norm(as.matrix(eig_sub),'2')
    
  }
  
  bindex=which(eig_blnorm<tol)
  
  for (i in bindex){
    
    if(i==1){
      eindex=1:pseq[1]
    }else{
      eindex=(sum(pseq[1:(i-1)])+1):(sum(pseq[1:i]))   #rows corresponding to block(k,l)
    }
    
    eig[eindex]<-rep(0,length(eindex))
  }
  return (eig)
}





FN_block<-function(eig,v,pseq,tol){
  
  I=length(pseq)
  fac=v_fac(pseq) 
  eig_M=cbind(eig,fac)
  v_M=cbind(v,fac)
  eig_blnorm<-rep(0,I)
  v_blnorm<-rep(0,I)
  
  for(i in 1:I){
    eig_sub=eig_M[which(eig_M[,2]==i),1]
    v_sub=v_M[which(v_M[,2]==i),1]
    eig_blnorm[i]=norm(as.matrix(eig_sub),'2')
    v_blnorm[i]=norm(as.matrix(v_sub),'2')  
  }
  
  fn=0
  index=which(v_blnorm!=0)
  
  for(i in index){
    if(abs(eig_blnorm[i])<tol){
      fn=fn+1
    }
  }
  
  return(fn)
  
}


adp_weight_elem<-function(eig,tol){
  
  M<-eig%*%t(eig)
  
  vecM=abs(c(M))
  vecM[which(vecM<tol)]=tol
  M.new=matrix(vecM,nrow(M),ncol(M))
  w_e=mean(c(M.new))/abs(M.new)
  
  return(w_e)
  
}

adp_weight_block<-function(eig,pseq,tol){
 
   M<-eig%*%t(eig)
  
  adp.bl<-block_norm(M,pseq)
  #vectorize
  vec.adp.bl=c(adp.bl)
  vec.adp.bl[which(vec.adp.bl<tol)]=tol
  adp.bl=matrix(vec.adp.bl,nrow(adp.bl),ncol=ncol(adp.bl))
  adp.bl.weg=mean(c(adp.bl))/adp.bl
  w_b=w*adp.bl.weg
  
  return(w_b)
  
}

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


#bool: true is the top half, ow,the lower half

gen_v1<-function(pseq,J,j0,bool=TRUE){
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
        indice=indice[1]:indice[length(indice)/2]
      }
      else{
        indice=indice[length(indice)/2+1]:indice[length(indice)]
      }
    }
    
    
    nnz_index=sample(indice,0.5*length(indice))
    
    v[nnz_index]=runif(length(nnz_index),min=0.5,max=1)
    v[nnz_index]=v[nnz_index]/(sqrt(sum(v[nnz_index]^2))*sqrt(length(J))) #normalize
  }
  
  return(v)
}


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

standardizeX <- function(X, pvec, center = F){
  d <- length(pvec)
  norms <- rep(0,d)
  svec <- rep(0,d)
  pcum <- cumsum(pvec)
  # Center each column
  if (center){
    Xmean <- colMeans(X)
    X <- X - matrix(Xmean, nrow(X), ncol(X), byrow = T) 
  }else{
    Xmean <- rep(0, ncol(X))
  }
  for (i in 1:d){
    if (i == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[i-1]+1):pcum[i])
    }
    norms[i] <- sum(X[,index]^2)
    # Scale the dataset
    X[, index] <- X[,index]/sqrt(norms[i])
    # Calculate largest singular value
    svec[i] <- max(svd(X[, index])$d)
  }
  return(list(X = X, svec = svec, norms = norms, Xmean = Xmean))
}




