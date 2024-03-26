getRho<-function(Y,PrevPi_Eig=NULL,PrevPi_d=NULL,pseq,w=NULL,w_e=NULL,rho,alphaseq,rho0seq,fold,ndim,eps,maxiter=10000){
  #Y:n*p observation matrix
  N=dim(Y)[1]
  #generate indice for fold
  indice=Kindice(N,fold)
  
  
  n_alpha=length(alphaseq)
  n_rho0=length(rho0seq)
  rho1.rho2.candi<-rep(0,2)
  
  for(i in 1:n_alpha){
    
    unit=cbind(alphaseq[i]*rho0seq,(1-alphaseq[i])*rho0seq)
    rho1.rho2.candi=rbind(rho1.rho2.candi,unit)
    
  }
  rho1.rho2.candi=rho1.rho2.candi[-1,]
  
  #construct the candidate set for rho1 and rho2, combination of alpha and rho
  
  Vsplit=apply(rho1.rho2.candi,1,FUN=function(u){
    cv_inner(Y=Y,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,pseq=pseq,w=w,w_e=w_e,rho=rho,rho1=u[1],rho2=u[2],indice=indice,ndim=ndim,eps=eps,maxiter=maxiter)
  }) 
  
  
  
  return(Vsplit)
}