# BlockProx<-function(M_kl,w_kl,we_kl,rho,rho1,rho2){
#   
#   
#   A_kl=matrix(0,nrow(M_kl),ncol(M_kl))
#   
#   #first step:elementwise soft thresholding
#   for(i in 1:nrow(M_kl)){
#     for(j in 1:ncol(M_kl)){
#       A_kl[i,j]=SoftThreshold(M_kl[i,j],we_kl[i,j]*rho1/rho)
#     }
#   }
#   
#   #groupwise thresholding
#   if(norm(A_kl,'F')==0){
#     cons=0
#   }
#   else{
#     cons=max(0,1-rho2*w/(rho*norm(A_kl,'F')))
#   }
#   
#   return(cons*A_kl)
#   
# }

BlockProx<-function(W,H1,w,w_e,rind,cind,rho,rho1,rho2){
  
  M_kl=(1/rho)*W[rind,cind]+H1[rind,cind]
  we_kl=w_e[rind,cind]
  A_kl=matrix(0,length(rind),length(cind))
  
  #first step:elementwise soft thresholding
  for(i in 1:length(rind)){
    for(j in 1:length(cind)){
      A_kl[i,j]=SoftThreshold(M_kl[i,j],we_kl[i,j]*rho1/rho)
    }
  }
  
  #groupwise thresholding
  if(norm(A_kl,'F')==0){
    cons=0
  }
  else{
    cons=max(0,1-rho2*w/(rho*norm(A_kl,'F')))
  }
  
  return(cons*A_kl)
  
}

