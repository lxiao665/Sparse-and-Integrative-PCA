FantopeProj<-function(mat,PrevPi_Eig=NULL,PrevPi_d=NULL,ndim){
  
  #PrevPi==NULL, then do fantope projection, otherwise, do deflated fantope projection
  if(is.null(PrevPi_Eig)==FALSE){
    p=dim(PrevPi)[1]
    # Eig=eigen(diag(1,p)-PrevPi,symmetric=TRUE)
    U=PrevPi_Eig[[2]]
    U=U[,(1:(p-PrevPi_d))]  #U : the orthogonal complement basis of PrevPi
    #mat=t(U)%*%mat%*%U    #the matrix that needs to be projected to the fantope
    mat=eigenMapMatMult3(t(U),mat,U)
    mat=(t(mat)+mat)/2    #force to be symmetric
    
  }
  
  #spectral decomposition of mat
  # U_Eig=eigen(mat,symmetric=TRUE)
  
  # U_Eig=getEigenDecomp(mat)
  
  #U_value=U_Eig[[1]]  #extract eigenvalue
  
  #if using  CPP, note that we should reverse the order of eigenvalues
  # U_value=rev(U_value)
  
  # U_vec=U_Eig[[2]]    #extract eigenvector
  
  #if  using cpp, reverse the order of eigenvectors
  #U_vec=U_vec[,ncol(U_vec):1]
  
  #U_value=eigen(mat,symmetric=TRUE,only.values = TRUE)$values
  
  
  
  #nnz_index=which(new.values!=0)
  
  # U_vec=eigs_sym(A=mat,k=length(nnz_index))$vectors
  
  U_eig=eigs_sym(A=mat,k=10)
  U_value=U_eig$values
  U_vec=U_eig$vectors
  
  theta=GetTheta(U_value,ndim)       #get theta by solving corresponding piecewise linear equation
  
  #compute the "eigenvalues" of the projection
  new.values=unlist(sapply(U_value,FUN=function(u){
    min(max(u-theta,0),1)
  }))
  #form the projection
  
  #newmat=U_vec%*%diag(new.values)%*%t(U_vec)
  # newmat=U_vec%*%Diagonal(x=new.values)%*%t(U_vec)
  #newmat=as.matrix(newmat)
  
  newmat=eigenMapMatMult3(U_vec,diag(new.values),t(U_vec))
  
  #newmat=eigenMapMatMult3(U_vec,diag(new.values),t(U_vec))
  
  
  if(is.null(PrevPi_Eig)==FALSE){
    # newmat=U%*%newmat%*%t(U)
    newmat=eigenMapMatMult3(U,newmat,t(U))
  }
  
  return(newmat)
  
}

