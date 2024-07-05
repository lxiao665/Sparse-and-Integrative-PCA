spca<-function(Y,S,PrevPi_Eig=NULL,PrevPi_d=NULL,alphaseq,rho0seq,pseq,w=NULL,w_e=NULL,rho,fold,ndim,t_cv,t_est,K,rhostep,eps_cv,eps_est,parallel=F,no_cores=8){
  
  #cross validation
  if(parallel==F){
    param=CV(Y=Y,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,pseq=pseq,w=w,w_e=w_e,rho=rho,alphaseq=alphaseq,rho0seq=rho0seq,fold=fold,ndim=ndim,t=t_cv,K=K,rhostep=rhostep,eps=eps_cv)
  }else{
    param=CV_parallel(Y=Y,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,pseq=pseq,w=w,w_e=w_e,rho=rho,alphaseq=alphaseq,rho0seq=rho0seq,fold=fold,ndim=ndim,t=t_cv,K=K,rhostep=rhostep,eps=eps_cv,no_cores=no_cores)
  }
  
  
  
  alpha_cv=param$alpha_cv
  rho0_cv=param$rho0_cv
  rho1_cv=alpha_cv*rho0_cv
  rho2_cv=(1-alpha_cv)*rho0_cv
  
  cv_sum=param$cv_sum
  sd_cv=param$sd_cv
  
  
  H=LA_seqadmm(PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,S=S,pseq=pseq,w=w,w_e=w_e,rho=rho,rho1=rho1_cv,rho2=rho2_cv,ndim=ndim,t=t_est,K=K,rhostep=rhostep,eps=eps_est)
  eig=getEigenDecomp(H)[[2]][,ncol(H)]
  
  
  Vsplit = param$Vsplit
  
  return(list(Vsplit = Vsplit, alpha_cv= alpha_cv, rho0_cv= rho0_cv, eig=eig, cv_sum=cv_sum, sd_cv=sd_cv))
  
}


