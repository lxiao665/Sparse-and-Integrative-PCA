setwd('/SIPCA_functions')


library(MASS)
library(Rcpp)
library(RSpectra)
library(RMTstat)


source('auxilliaryfunction.r')
source('BlockProx.r')
source('cv_inner_noscale_cv2.r')
source('CV.r')
source('FantopeProj_k10.r')
source('getRho.r')
source('LA_seqadmm.r')
source('seqADMM_perstage.r')
source('spca.r')
sourceCpp('eigendecomp.cpp')
sourceCpp('MatrixMtp.cpp')


###################### Generate signal/noise eigenvectors and covariance matrix with heteroskedastic noise ###########
# I: number of blocks
I=10
# p0: dimension of each block
p0=50
pseq=rep(p0,I)       #[P1,P2,...PI]
pcum <- cumsum(pseq)
pcum_ex <- c(0, pcum)
p=sum(pseq)
ndim=1
# eigenvalues
lambda1=40
lambda2=20

# n: sample size
n=200

#index of nonzero data blocks for the first two eigenvectors
Jv1<-1:3
Jv2<-4:6

# generate sparse eigenvectors of the signal matrix
i = 12345
set.seed((i))

v1=gen_v_normal(pseq,Jv1,alpha=0.5,i)    #leading eigenvector
v2=gen_v_normal(pseq,Jv2,alpha=0.5,i)    #second eigenvector
V=cbind(v1,v2)

lambda=c(lambda1,lambda2)
Lambda=diag(lambda)
Sigma_signal = V%*%Lambda%*%t(V)
# mean vector
mu = rep(0,p)
Lambda=diag(lambda)
# signal covariancee
Sigma_signal = V%*%Lambda%*%t(V)
Y_signal = mvrnorm(n,mu,Sigma_signal) 


# generate noise matrix with homoskedastic noise level across different data blocks

sigma2 <- 0.5

Sigma_noise = diag(rep(sigma2, p))
Y_noise = mvrnorm(n,mu,Sigma_noise) 
Sigma = Sigma_signal+Sigma_noise

# observational matrix
Y = Y_signal + Y_noise


########################## Specifying hyperparameters ######################
# weight in block-wise penalty
w=gen_weight(pseq)
# parameters in admm algorithm
rho=100 
t_cv=5
t_est=5
rhostep=2
K=50
eps_cv=1e-2
eps_est=1e-3

# cross validation
fold=5
n_alpha=5
n_rho0=10


##################### Estimate sparse eigenvectors #########################
Y = scale(Y,scale=FALSE)  #demean
S = cov(Y)
#### leading eigenvector###
# candidates for cross validation tuning parameters
alphaseq1=seq(from=0,to=1,length.out=n_alpha)
rho0seq1=GenerateRho0(Sj=S,n_rho0=n_rho0,minv=0.01)
out1<-spca(Y=Y,S=S,PrevPi_Eig=NULL,PrevPi_d=NULL,alphaseq=alphaseq1,rho0seq=rho0seq1,pseq=pseq,w=w,
           rho=rho,fold=fold,ndim=ndim,t_cv=t_cv,t_est=t_est,K=K,rhostep=rhostep,eps_cv=eps_cv,eps_est=eps_est)

eig1=out1$eig



#### second eigenvector ####

PrevPi=eig1%*%t(eig1)
PrevPi_Eig=getEigenDecomp(diag(1,p)-PrevPi)

PrevPi_Eig[[1]]<-rev(PrevPi_Eig[[1]])
PrevPi_Eig[[2]]<-PrevPi_Eig[[2]][,ncol(PrevPi_Eig[[2]]):1]

PrevPi_d=1

Sj=eigenMapMatMult3((diag(1,p)-PrevPi),S,(diag(1,p)-PrevPi))
alphaseq2=seq(from=0,to=1,length.out=n_alpha)
rho0seq2=GenerateRho0(Sj=Sj,n_rho0=n_rho0,minv=0.01)

out2<-spca(Y=Y,S=S,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,alphaseq=alphaseq2,rho0seq=rho0seq2,pseq=pseq,w=w,
           rho=rho,fold=fold,ndim=ndim,t_cv=t_cv,t_est=t_est,K=K,rhostep=rhostep,eps_cv=eps_cv,eps_est=eps_est)

eig2=out2$eig


# ####################### Calculate the eigenspace estimation error #################
V <- cbind(v1, v2)
V_est <- cbind(eig1, eig2)

l2_err = sqrt(sum((crossprod(t(V),t(V))-crossprod(t(V_est),t(V_est)))^2)/sum((crossprod(t(V),t(V)))^2))


  

  
 

