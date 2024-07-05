seqADMM_perstage<-function(H20=NULL,S,PrevPi_Eig=NULL,PrevPi_d=NULL,pseq,w=NULL,w_e=NULL,rho=NULL,rho1=NULL,rho2=NULL,ndim,t){
  
  
  p=dim(S)[1]               #dimension of sample covariance matrix
  I=length(pseq)            #number of datasets
  
  
  if(is.null(w_e)==TRUE){
    w_e=matrix(1,p,p)
  }
  
  if(is.null(w)==TRUE){
    w=matrix(1,I,I)
  }
  
  #starting value
  niter=0
  
  
  #if no initial value of H20 is given, then the default value is zero.
  if(is.null(H20)==TRUE){
    H20=matrix(0,p,p)
  }
  
  
  H2=matrix(0,p,p)   ## used to store  value for H2 at the current step
  
  W0=matrix(0,p,p)   ##intial value for dual variable , used to store value of dual variable at the previous step
  W=matrix(0,p,p)    ##used to store value of dual variable at the current step
  
  H1_sum=matrix(0,p,p)  ## used to store the sum of all H1 in each iteration of this stage
  H2_sum=matrix(0,p,p) ## used to store the sum of all H2 in each iteration of this stage
  W_sum=matrix(0,p,p) ## used to store the sum of all dual variable W in each iteration of this stage
  
  #iteration
  while(niter<t){
    #update primal variable
    ##update H_1
    mat=H20-(1/rho)*(W0-S)
    mat=(mat+t(mat))/2    #force symmetric
    H1=FantopeProj(mat=mat,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,ndim=ndim)
    
    ##update H_2
    #op==1,with elementwise 1-1 norm
 
    
    
    #op=2 or op=3,update each block of H2 respectively
    
   
      ##identify the rows and columns index corresponding to block(k,l)
      
      for (k in 1:I){
        for (l in k:I){
          
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
          
          #update the block (k,l) of H2
          
          ## with block-wise 1-1 penalty 
     
          
          
          ##with element-wise and block-wise 1-1 penalty
          # M_kl=(1/rho)*W0[rind,cind]+H1[rind,cind]
          # we_kl=w_e[rind,cind]
          # 
          # H2[rind,cind]=BlockProx(M_kl,w[k,l],we_kl,rho,rho1,rho2)    
            
          H2[rind,cind]=BlockProx(W0,H1,w[k,l],w_e=w_e,rind,cind,rho,rho1,rho2)                                     
          
          
          if(l!=k){
            H2[cind,rind]=t(H2[rind,cind])   #H2 should be symmetric
          }
        }
      
      
    }
    
    #update dual variable W
    W=W0+rho*(H1-H2)
    
    #reset H20,W0
    H20=H2
    W0=W
    
    
    #update sum
    H1_sum=H1_sum+H1
    H2_sum=H2_sum+H2
    W_sum=W_sum+W
    
    niter=niter+1
    
  }
  
  return(list(H1_sum/t,H2_sum/t,W_sum/t))
  
  
}

