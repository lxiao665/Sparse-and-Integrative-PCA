setwd('/SIPCA_functions')

library(MASS)
library(Rcpp)
library(RSpectra)
library(doParallel)
library(doRNG)
library(RMTstat)


source("AuxillaryFunctions.R")
source('auxilliaryfunction.r')
source('BlockProx.r')
source('cv_inner_noscale_cv2_heter.r')
source('CV_heter.r')
source('CV_parallel_heter.r')
source('FantopeProj_k10.r')
source('getRho_heter.r')
source('LA_seqadmm.r')
source('seqADMM_perstage.r')
source('spca_heter.r')

sourceCpp('eigendecomp.cpp')
sourceCpp('MatrixMtp.cpp')

# Load pre-processed TCGA BRCA data
geneExp <- read.table("/TCGA data/Preprocessed_GeneExp.txt", sep = "\t")
#dim(geneExp) # 645 by 348
meth <- read.table("/TCGA data/Preprocessed_Meth.txt", sep = "\t")
#dim(meth) # 574 by 348
miRNA <- read.table("/TCGA data/Preprocessed_miRNA.txt", sep = "\t")
#dim(miRNA) # 423 by 348
protein <- read.table("/TCGA data/Preprocessed_Protein.txt", sep = "\t")
#dim(protein) # 171 by 348



#################### Concatenate the data and form pvec #####################
pseq <- c(nrow(geneExp), nrow(meth), nrow(miRNA), nrow(protein))
Xall <- cbind(t(geneExp), t(meth), t(miRNA), t(protein))
rem_idx <- c(50, 137, 144, 181, 308)
Xall <- Xall[-rem_idx, ]
n <- nrow(Xall)


#################### hyperparameters ##########################
rho=100 #initial value
fold=5
K=50
t_cv=5
t_est=5
rhostep=2
I=4

pcum <- cumsum(pseq)
pcum_ex <- c(0, pcum)
p=sum(pseq)
ndim=1

w=gen_weight(pseq)

eps_cv=1e-2
eps_est=1e-3

n_alpha = 1
n_rho0 = 25
n_eig = 53


############################## Normalization ###############################
standardizeX <- function(X, pvec, center = F){
  n <- nrow(X)
  d <- length(pvec)
  norms <- rep(0,d)
  svec <- rep(0,d)
  pcum <- cumsum(pvec)
  # Center each column
  if (center){
    Xmean <- colMeans(X)
    X <- X - matrix(Xmean, nrow(X), ncol(X), byrow = T) 
  }else{
    Xmean <- rep(0, ncol(X))
  }
  for (i in 1:d){
    if (i == 1){
      index <- c(1:pvec[1])
    }else{
      index <- c((pcum[i-1]+1):pcum[i])
    }
    norms[i] <- sum(X[,index]^2)
    # Scale the dataset
    X[, index] <- n*X[,index]/sqrt(norms[i])
    # Calculate largest singular value
    svec[i] <- max(svd(X[, index])$d)
  }
  return(list(X = X, svec = svec, norms = norms, Xmean = Xmean))
}

out_s <- standardizeX(Xall, pseq, center = T)
Y <- out_s$X
S <- cov(Y)


############################# Estimating top n_eig eigenvectors #########################
output <- vector('list', length = n_eig)
eig_list <- vector('list', length = n_eig)


PrevPi_d=0
PrevPi = NULL


for (i in 1:length(output)){
  if(i==1){
    PrevPi_d=0
    PrevPi = NULL
    Sj=S
    PrevPi_Eig = NULL
  }else{
    PrevPi_d = PrevPi_d + 1
    if(is.null(PrevPi)){
      PrevPi =  eig%*%t(eig)
    }else{
      PrevPi =  PrevPi + eig%*%t(eig)
    }
    PrevPi_Eig=getEigenDecomp(diag(1,p)-PrevPi)
    
    PrevPi_Eig[[1]]<-rev(PrevPi_Eig[[1]])
    PrevPi_Eig[[2]]<-PrevPi_Eig[[2]][,ncol(PrevPi_Eig[[2]]):1]
    
    Sj=eigenMapMatMult3((diag(1,p)-PrevPi),S,(diag(1,p)-PrevPi))
  }
  
  
  
  alphaseq=seq(from=0,to=1,length.out=n_alpha)
  rho0seq=GenerateRho3(Sj=Sj,n_rho0=n_rho0,minv=0.01)
  rho0seq=c(0,rho0seq)
  
  
  out<-spca(Y=Y,Sigma0_est=sigma2_est_diag,S=S,PrevPi_Eig=PrevPi_Eig,PrevPi_d=PrevPi_d,alphaseq=alphaseq,rho0seq=rho0seq,pseq=pseq,w=w,
            rho=rho,fold=fold,ndim=ndim,t_cv=t_cv,t_est=t_est,K=K,rhostep=rhostep,eps_cv=eps_cv,eps_est=eps_est,parallel=T)
  output[[i]] <- out
  
  eig<- out$eig_1sd
  
  eig_list[[i]] <- eig
  
  
}







