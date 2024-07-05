


LA_seqadmm<-function(H20=NULL,S,PrevPi_Eig=NULL,PrevPi_d=NULL,pseq,w=NULL,w_e=NULL,rho=NULL,rho1=NULL,rho2=NULL,ndim,t,K,rhostep,eps=1e-2){
  
  norm1=rep(0,K-2)
  norm2=rep(0,K-2)
  
  
  
  
  for (i in 1:K){
    
    H=seqADMM_perstage(H20=H20,S=S,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,pseq=pseq,w=w,w_e=w_e,rho=rho,rho1=rho1,rho2=rho2,ndim=ndim,t=t)    
    
    
    H1=H[[1]]
    H2=H[[2]]
    W=H[[3]]
    
    rho=rho*rhostep
    
    if(i>2){
      
      #stopping criterion
      norm1[i-2]=(norm(H1-H2,'F'))^2    #primal residual
      norm2[i-2]=(norm(H2-H20,'F'))^2   #dual residual
      maxnorm=max(norm1[i-2],norm2[i-2])
      
      if(maxnorm<eps^2){
        break
      }
      
    }
    
    #reset parameter
    
    H20=H2     #store primal variable
    W20=W    #store dual variable
    
  }
  
  if(i==K){
    return(-1)    #not converge
  }
  
  else{
    return(H2)
  }
  
}
