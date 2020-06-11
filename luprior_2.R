## luprior_2.R
## log prior distribution of the parameters
## ver 2) fix transformation 

## arguments:
## x : parameter on the uniform and log transformed range
## cv : scalar or vector of prior cv for normal or truncated normal distribution
## hyper: defining the range and whether log-transformation is needed
## value: value of truncated and independent normal densities based on the range
require(truncnorm)
dtlnorm <- function(x,cv,l,r,log=TRUE)
{
  if(log){
    bound <- log(c(l,r))
    mu <- mean(bound)
    sigma <-  sqrt(log(1+(cv/100)^2))
    output <- dtruncnorm(log(x),a=bound[1],b=bound[2],mean=mu,sd=sigma)
  }
  else{
    mu <- mean(c(l,r))
    output <- dtruncnorm(x,a=l,b=r,mean=mu,sd=cv*mu/100)
  }
  log(output)
}
test1 <- function()
{
  rm(list=ls())
  source("luprior_1.R")
  xseq <- seq(0.01,0.21,len=99)
  dseq <- dtlnorm(xseq,50,0.01,0.21,log=TRUE)
  plot(xseq,dseq,type="l",ylim=c(0,1))
  dseq2 <- dtlnorm(xseq,50,0.01,0.21,log=FALSE)
  plot(xseq,dseq2,type="l")
}
luprior <- function(    ###CORRECT DEPTH and TOTAL MICROSITE HERE
  x,cv=500,hyper=list(
    KmO2_Rh=c(2.45,24500,1),
    alphaVmaxSx=c(2300000000,23000000000000,1),
    EaVmaxSx=c(7.2,7200,1),
    KMSx=c(10,100000,1),
    Depth=c(5,30,1),
    total.microsite=c(1000,100000,1),
    alphaVmaxCH4prod=c(30000,300000000,1),
    EaVmaxCH4prod=c(20,100,1),
    KMSx_CH4prod=c(3.5,35000,1),
    alphaVmaxCH4ox=c(0.007,700,1),
    EaVmaxCH4ox=c(10,50,0),
    KMCH4ox=c(1.00E-04,1.00,1),
    KmO2_CH4ox=c(0.043,4300,1),
    kl_CH4prod=c(0.03,300,1)
  )
)
{
  stopifnot(length(x)==length(hyper))
  
  ## Parameter Transformation
  # stopifnot(!any(is.na(x)))
  z <- damm_tran(x=x,hyper=hyper) ## back-transform to restricted range
  
  output <- rep(NA,length(x))
  if(length(cv)==1){
    cv <- rep(cv,length(x))
  }
  for(i in 1:length(x)){
    output[i] <- dtlnorm(z[i],cv[i],hyper[[i]][1],hyper[[i]][2],log=as.logical(hyper[[i]][3]))
  }
  sum(output)
}

test2 <- function()
{
  rm(list=ls())
  source("damm_utils_3.R")
  source("luprior_2.R")
  hyper=list(
    KmO2_Rh=c(0.01,0.21,0),
    alphaVmaxSx=c(1666666.66666667,166666666.666667,1),
    EaVmaxSx=c(30,100,0),
    KMSx=c(8.00E-08,8.00E-06,1),
    alphaVmaxCH4prod=c(1722222222222.22,172222222222222,1),
    EaVmaxCH4prod=c(15,75,0),
    KMSx_CH4prod=c(1666666.66666667,166666666.666667,1),
    alphaVmaxCH4ox=c(8.32660531276481,832.660531276481,1),
    EaVmaxCH4ox=c(10,50,0),KMCH4ox=c(5,50,0),
    KmO2_CH4ox=c(0.01,0.21,1),
    kl_CH4prod=c(0.001,0.01,1)
  )
  x0 <-   c(
    lKmO2_Rh=logit_tran(0.1,0.01,0.21),
    lalphaVmaxSx=logit_tran(16666666.6666667,hyper[["alphaVmaxSx"]][1],hyper[["alphaVmaxSx"]][2],log=TRUE), 
    lEaVmaxSx=logit_tran(68,hyper[["EaVmaxSx"]][1],hyper[["EaVmaxSx"]][2]),
    lKMSx=logit_tran(0.0000008,hyper[["KMSx"]][1],hyper[["KMSx"]][2],log=TRUE),
    lalphaVmaxCH4prod=logit_tran((6.2*(10^16)*(1/3600)),hyper[["alphaVmaxCH4prod"]][1],hyper[["alphaVmaxCH4prod"]][2],log=TRUE), 
    lEaVmaxCH4prod=logit_tran(50,hyper[["EaVmaxCH4prod"]][1],hyper[["EaVmaxCH4prod"]][2]),
    lKMSx_CH4prod=logit_tran(2*16666666.6666667,1666666.66666667,166666666.666667,log=TRUE),
    lalphaVmaxCH4ox=logit_tran(299757.791259533*(1/3600),hyper[["alphaVmaxCH4ox"]][1],hyper[["alphaVmaxCH4ox"]][2],log=TRUE), 
    lEaVmaxCH4ox=logit_tran(30,hyper[["EaVmaxCH4ox"]][1],hyper[["EaVmaxCH4ox"]][2]),
    lKMCH4ox=logit_tran(7,hyper[["KMCH4ox"]][1],hyper[["KMCH4ox"]][2]),
    lKmO2_CH4ox=logit_tran(0.10,0.01,0.21),
    lkl_CH4prod=logit_tran(0.005,0.001,0.01,log=TRUE)
  )
  luprior(x0,cv=1)
  #[1] -Inf
  luprior(x0,cv=100)
  #[1]-19.47363
}

mcmctest <- function()
{
  rm(list=ls())
  rm(list=ls())
  source("damm_utils_3.R")
  source("luprior_2.R")
  hyper=list(
    KmO2_Rh=c(0.01,0.21,0),
    alphaVmaxSx=c(1666666.66666667,166666666.666667,1),
    EaVmaxSx=c(30,100,0),
    KMSx=c(8.00E-08,8.00E-06,1),
    alphaVmaxCH4prod=c(1722222222222.22,172222222222222,1),
    EaVmaxCH4prod=c(15,75,0),
    KMSx_CH4prod=c(1666666.66666667,166666666.666667,1),
    alphaVmaxCH4ox=c(8.32660531276481,832.660531276481,1),
    EaVmaxCH4ox=c(10,50,0),KMCH4ox=c(5,50,0),
    KmO2_CH4ox=c(0.01,0.21,1),
    kl_CH4prod=c(0.001,0.01,1)
  )
  x0 <-   c(
    lKmO2_Rh=logit_tran(0.1,0.01,0.21),
    lalphaVmaxSx=logit_tran(16666666.6666667,hyper[["alphaVmaxSx"]][1],hyper[["alphaVmaxSx"]][2],log=TRUE), 
    lEaVmaxSx=logit_tran(68,hyper[["EaVmaxSx"]][1],hyper[["EaVmaxSx"]][2]),
    lKMSx=logit_tran(0.0000008,hyper[["KMSx"]][1],hyper[["KMSx"]][2],log=TRUE),
    lalphaVmaxCH4prod=logit_tran((6.2*(10^16)*(1/3600)),hyper[["alphaVmaxCH4prod"]][1],hyper[["alphaVmaxCH4prod"]][2],log=TRUE), 
    lEaVmaxCH4prod=logit_tran(50,hyper[["EaVmaxCH4prod"]][1],hyper[["EaVmaxCH4prod"]][2]),
    lKMSx_CH4prod=logit_tran(2*16666666.6666667,1666666.66666667,166666666.666667,log=TRUE),
    lalphaVmaxCH4ox=logit_tran(299757.791259533*(1/3600),hyper[["alphaVmaxCH4ox"]][1],hyper[["alphaVmaxCH4ox"]][2],log=TRUE), 
    lEaVmaxCH4ox=logit_tran(30,hyper[["EaVmaxCH4ox"]][1],hyper[["EaVmaxCH4ox"]][2]),
    lKMCH4ox=logit_tran(7,hyper[["KMCH4ox"]][1],hyper[["KMCH4ox"]][2]),
    lKmO2_CH4ox=logit_tran(0.10,0.01,0.21),
    lkl_CH4prod=logit_tran(0.005,0.001,0.01,log=TRUE)
  )
  require(mcmc)
  (z0 <- damm_tran(x0,hyper=hyper))
  luprior(x0,cv=50)
  #[1] -16.6281
  out <- metrop(luprior,initial=x0,nbatch=1000,cv=50)
  out$accept # 0.04
  out <- metrop(out,scale=0.5,cv=50)
  out$accept # 0.224
  burnin <- metrop(out,nbatch=5000,nspac=2,cv=50,scale=0.5)
  burnin$accept # 0.2364
  final <- metrop(burnin,nbatch=5000,nspac=2,cv=50,scale=0.5)
  final$accept # 2396
  plot(ts(final$batch[,1:6]),plot.type="multiple")
  plot(ts(final$batch[,6+1:6]),plot.type="multiple")
  x1 <- final$final
  (z1 <- damm_tran(x1))
}
