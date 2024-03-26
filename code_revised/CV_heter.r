CV<-function(Y,Sigma0_est=NULL,PrevPi_Eig=NULL,PrevPi_d=NULL,pseq,w=NULL,w_e=NULL,rho,alphaseq,rho0seq,fold,ndim,t,K,rhostep,eps){
 
  
  
 n_rho0=length(rho0seq)
 n_alpha=length(alphaseq)
  
  
  Vsplit=getRho(Y=Y,Sigma0_est=Sigma0_est,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,pseq=pseq,w=w,w_e=w_e,rho=rho,alphaseq=alphaseq,rho0seq=rho0seq,fold=fold,ndim=ndim,t=t,K=K,rhostep=rhostep,eps=eps)
  
  
  #one standard error rule
  cv_sum = Re(apply(Vsplit, 2, sum))
  cv_mean = Re(apply(Vsplit, 2, mean))
  max_index=which(cv_sum==max(cv_sum))[1]
  
  sd_cv=sd(Vsplit[,max_index])/sqrt(fold)
  
  
  
  if(max_index%%n_rho0!=0){
    
    alphaindex=max_index%/%n_rho0+1
    rho0index=max_index%%n_rho0
    
    alpha_cv=alphaseq[alphaindex]
    
    rho0_cv=rho0seq[rho0index]
  }else{
    
    alphaindex=max_index%/%n_rho0
    rho0index=n_rho0
    alpha_cv=alphaseq[alphaindex]
    
    rho0_cv=rho0seq[rho0index]
    
  }
  

  
  m_cv_mean=t(matrix(cv_mean,n_rho0,n_alpha))

  max_cv_mean=m_cv_mean[alphaindex,]

  rho0_1sd_index=max(which(max_cv_mean>max(cv_mean)-1*sd_cv))
  
  alpha_cv_1sd=alpha_cv
  rho0_cv_1sd=rho0seq[rho0_1sd_index]
  
  
  
 return(list(Vsplit=Vsplit,alpha_cv=alpha_cv,rho0_cv=rho0_cv,alpha_cv_1sd=alpha_cv_1sd,rho0_cv_1sd=rho0_cv_1sd,cv_sum=cv_sum,sd_cv=sd_cv))
  
  
}