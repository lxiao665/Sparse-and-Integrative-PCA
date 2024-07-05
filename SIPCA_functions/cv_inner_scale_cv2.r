cv_inner<-function(Y,PrevPi_Eig=NULL,PrevPi_d=NULL,pseq,w=NULL,w_e=NULL,rho,rho1=NULL,rho2=NULL,indice,ndim,t,K,rhostep,eps){
  fold=length(indice)
  #a vector used to store the cross validated inner product for each fold.
  V=rep(0,fold)
  
  #findice:fold indice
  for(i in 1:length(indice)){
    Y_train=Y[-indice[[i]],]
    Y_test=Y[indice[[i]],]
    Strain=cov(Y_train)  #sample covariance matrix for training dataset
    Stest=cov(Y_test)   #sample covariance matrix for test dataset
    
    tr_Strain=sum(diag(Strain))
    Strain=Strain*(p/tr_Strain)
    
    tr_Stest=sum(diag(Stest))
    Stest=Stest*(p/tr_Stest)
    
    n_indice=length(indice[[i]])
    ProjH=LA_seqadmm(S=Strain,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,pseq=pseq,w=w,w_e=w_e,rho=rho,rho1=rho1,rho2=rho2,ndim=ndim,t=t,K=K,rhostep=rhostep,eps=eps)
    
    ProjH.eig=getEigenDecomp(ProjH)[[2]][,ncol(ProjH)]
    
    V[i]=eigenMapMatMult3(t(ProjH.eig),Stest,ProjH.eig)
    
  }
  
  return (V)
  
}