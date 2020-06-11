## lupost_2.R
##  log un-normalized posterior density
## ver 2) return scalar values
## Arguments:
## x     : transformed parameters on the real line
##
## Prior Arguments:
## cv : scalar or vector of prior cv for normal or truncated normal distribution
## hyper: defining the range and whether log-transformation is needed
##
## Additional Likelihood Arguments:
## formula: a one sided formula defining the observed soil moisture
##         soil temperature, CO2 and methane in that order
## data : input data frame
## nlevels: number of levels for bivariate distribution (check consistency)
## R, O2airfrac,  BD, PD,  Porosity,
## p, RQ,  Depth,
## microsite,total.microsite,time.step,
## flagO2=TRUE,  flagCH4=TRUE, ## whether run O2 and CH4 to steady states
## iO2.lst=NULL,iCH4.lst=NULL,
## reltol=sqrt(.Machine$double.eps),verbose=0,
## sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
## soilm.cv=20, soilm.range=c(70,200), soilm.nlevel=10, lognorm=TRUE,## iO2.lst: a list of initial O2 values (default 0.209?)
## iCH4.lst: a list of initial CH4 values (default 1.8?)
## nu   : prior degree of freedom for 2x2 precision (default=2)
## R0   : prior concentration matrix for 2x2 precision (default identity)
## value
##    scalar un-normalized posterior density

lupost <- function(
  x,cv,formula,data,nlevels,
  R, O2airfrac,  BD, PD,  Porosity,
  p, RQ,  Depth,    ###CORRECT DEPTH and TOTAL MICROSITE HERE
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=TRUE, ## whether run O2 and CH4 to steady states
  iO2.lst=NULL,iCH4.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,200), soilm.nlevel=10, lognorm=TRUE,
  nu=2,R0=diag(2),
  hyper=list(   ###CORRECT DEPTH and TOTAL MICROSITE HERE
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
){
  
  ## prior density
  llik1 <- luprior(x=x,cv=cv,hyper=hyper)
  ## log likelihood
  llik2 <- lulik(     ###CORRECT DEPTH and TOTAL MICROSITE HERE
    x=x,formula=formula,data=data,nlevels=nlevels,
    R=R,O2airfrac = O2airfrac, BD=BD,PD=PD, Porosity = Porosity,
    p=p, RQ=RQ, Depth= Depth, microsite = microsite,
    total.microsite = total.microsite, time.step = time.step,
    flagO2 = flagO2, flagCH4 = flagCH4,
    iO2.lst = iO2.lst, iCH4.lst = iCH4.lst,
    reltol = reltol, verbose = verbose,
    sxtot.range = sxtot.range, nquantile = nquantile,
    mult = mult, soilm.cv = soilm.cv, soilm.range = soilm.range,
    soilm.nlevel = soilm.nlevel, lognorm = lognorm,
    nu = nu, R0 = R0, hyper=hyper)
  
  ## un-normalized log posterior
  llik2$value <- llik1 + llik2$value
  llik2$value
}

mcmctest <- function()
{
  rm(list=ls())
  source("damm_utils_3.R")
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
  trench <- read.csv("2015data.csv",head=T)
  library(doParallel)
  source("damm_10c.R")
  source("luprior_2.R")
  source("lulik_3.R")
  source("lupost_2.R")
  
  # r0 <- lupost(
  #   x=x0,cv=50,formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
  #   R=0.008314472, O2airfrac=0.209, BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
  #   p=0.00046,RQ=1, Depth=10, microsite = 0,total.microsite = 10000, time.step = 300
  # )
  # 
  # r0$elapsed
  
  library(doParallel)
  cl <- makeCluster(31)
  registerDoParallel(cl)
  sink <- clusterEvalQ(cl,source("damm_10c.R"))
  sink <- clusterEvalQ(cl,source("damm_utils_3.R"))

  
  r1 <- lupost(
    x=x0,cv=50,formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth=10, microsite = 0,total.microsite = 10000, time.step = 300
  )

  ## first try optimization
  o1 <- optim(
    par=x0,fn=lupost, gr=NULL, method = "SANN",control=list(trace=9,fnscale=-1),
    cv=50,formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth=10, microsite = 0,total.microsite = 10000, time.step = 300)
  
  require(mcmc)  
  ## adaptation staget
  ## save iteration and keep burn-in to run on cluster, increase tolerance.
  nsweep <- 500 ## each calibration batch size
  tol <- 0.1    ## relative tolerance to target
  target <- 0.2   ## target acceptance rate
  delta <- 9999 + tol ## check
  scale__ <- 1    ## default proposal scale
  count__ <- 1    ## counter
  max_iter__ <- 9  ## maximum number of iterations to adapt
  while(delta > tol){
    if(count__>max_iter__){
      break
    }
    out <- metrop(
      obj=lupost,initial=x0,nbatch=nsweep,scale=scale__,  
      cv=50,formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
      R=0.008314472, O2airfrac=0.209, BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
      p=0.00046,RQ=1, Depth=10, microsite = 0,total.microsite = 10000, time.step = 300
    )
    save(out,file=paste("Scratch/r_adapt",count__,".rda",sep=""))
    ## update proposal scale based on acceptance rate
    cat("scale = ",scale__,"accept =",out$accept,"\n")
    if(out$accept>target){
      scale__ <- scale__ * 2
    }
    else{
      scale__ <- scale__ / 2
    }
    ## how close relative to target in magnitude
    delta <- abs(out$accept - target)/target 
  }
  # scale =  1 accept = 0.028
  # scale =  0.5 accept = 0.02
  # scale =  0.25 accept = 0.04
  # scale =  0.125 accept = 0.05
  # scale =  0.0625 accept = 0.032
  # scale =  0.03125 accept = 0.096
  # scale =  0.015625 accept = 0.154
  # scale =  0.0078125 accept = 0.346
  # scale =  0.015625 accept = 0.206
  
  ## burnin stage
  out <- metrop(
    obj=lupost,initial=x0,nbatch=burnin,nspac=thin,scale=scale__,  
    cv=50,formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth=10, microsite = 0,total.microsite = 10000, time.step = 300
  )
  save(out,file=paste("Scratch/rburnin",count__,".rda",sep=""))
  
  ## production stage
  out <- metrop(
    obj=out,nbatch=iters-burnin,nspac=thin,scale=scale__,
    cv=50,formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth=10, microsite = 0,total.microsite = 10000, time.step = 300
  )
  save(out,file=paste("Scratch/rprod",count__,".rda",sep=""))
  
  stopCluster(cl)
}