CV_parallel<-function(Y,PrevPi_Eig=NULL,PrevPi_d=NULL,pseq,w=NULL,w_e=NULL,rho,alphaseq,rho0seq,fold,ndim,t,K,rhostep,eps,no_cores){
 
  cv_inner<-function(j){
    rho1 = rho1.rho2.candi[j,1]
    rho2 = rho1.rho2.candi[j,2]
    
    fold=length(indice)
    #a vector used to store the cross validated inner product for each fold.
    V=rep(0,fold)
    
    #findice:fold indice
    for(i in 1:length(indice)){
      Y_train=Y[-indice[[i]],]
      Y_test=Y[indice[[i]],]
      Strain=cov(Y_train)  #sample covariance matrix for training dataset
      Stest=cov(Y_test)   #sample covariance matrix for test dataset
      
      
      n_indice=length(indice[[i]])
      ProjH=LA_seqadmm(S=Strain,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,pseq=pseq,w=w,w_e=w_e,rho=rho,rho1=rho1,rho2=rho2,ndim=ndim,t=t,K=K,rhostep=rhostep,eps=eps)
      
      ProjH.eig=getEigenDecomp(ProjH)[[2]][,ncol(ProjH)]
      
      V[i]=eigenMapMatMult3(t(ProjH.eig),Stest,ProjH.eig)
      
    }
    
    return (V)
    
  }
  
  
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
  
 n_rho0 = length(rho0seq)
 n_alpha = length(alphaseq)
 n_param = n_rho0 * n_alpha
 
 
 
 cl=makeCluster(no_cores,type="FORK")
 Vsplit=parSapply(cl, 1:n_param, cv_inner)
 stopCluster(cl)
  
  
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
 
  # 
  
 # return(list(Vsplit=Vsplit,alpha_cv=alpha_cv,rho0_cv=rho0_cv,alpha_cv_1sd=alpha_cv_1sd,rho0_cv_1sd=rho0_cv_1sd,cv_sum=cv_sum,sd_cv=sd_cv))
 return(list(Vsplit=Vsplit,alpha_cv=alpha_cv,rho0_cv=rho0_cv,alpha_cv_1sd=alpha_cv_1sd,rho0_cv_1sd=rho0_cv_1sd,cv_sum=cv_sum,sd_cv=sd_cv))
  
}