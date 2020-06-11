## lulik_utils_2.R
## sub ver b)
## minor change 9/28/17: be consistent with ver 4b) of lulik
## utility function to run non-linear least square using minpack.lm
## and the conventional Nelder-Mead algorithm using optim

## most arugments are defined in damm_xx.R and lulik_xx.R
## the first two arugments have been modified to allow piece-wise optimization
## par: the interested paramemter
## z  : the ancillary parameter
## for rh, CH4 updating was not performed
## value: the residual as defined by minpack.lm
lulik.rh <- function(         ###CORRECT DEPTH and TOTAL MICROSITE HERE
  par,z,formula,data,nlevels,
  R, O2airfrac, CH4airfrac, BD, PD,  Porosity,
  p, RQ,  Depth_O2, Depth_CH4, Depth,   
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=FALSE, ## whether run O2 and CH4 to steady states
  iO2.lst=NULL,iCH4.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,130), soilm.nlevel=10, lognorm=TRUE,
  nu=2,R0=diag(2),
  hyper=list(
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
  ))
{
  r <- lulik(    ###CORRECT DEPTH and TOTAL MICROSITE HERE
    x=c(par,z),formula=formula,data=data,nlevels=nlevels,
    R=R, O2airfrac=O2airfrac, CH4airfrac = CH4airfrac,
    BD=BD, PD=PD,  Porosity=Porosity,
    p=p, RQ=RQ,  Depth_O2 = Depth_O2, Depth_CH4 = Depth_CH4,Depth = Depth,   
    microsite=microsite,total.microsite=total.microsite,time.step=time.step,
    flagO2=flagO2,  flagCH4=flagCH4, 
    iO2.lst=iO2.lst,iCH4.lst=iCH4.lst,
    reltol=reltol,verbose=verbose,
    sxtot.range=sxtot.range, nquantile=nquantile, mult=mult,
    soilm.cv=soilm.cv, soilm.range=soilm.range, soilm.nlevel=soilm.nlevel, 
    lognorm=lognorm,
    nu=nu,R0=R0,
    hyper=hyper)
  
  resid <- r$data[,3] - r$fit_Rh #This indicates this is 3rd in l.mat
  resid[!is.na(resid)]
}

lulik.rh.optim <- function(   ###CORRECT DEPTH and TOTAL MICROSITE HERE
  par,z,formula,data,nlevels,
  R, O2airfrac, CH4airfrac, BD, PD,  Porosity,
  p, RQ,  Depth_O2, Depth_CH4,Depth,
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=FALSE, ## whether run O2 and CH4 to steady states
  iO2.lst=NULL,iCH4.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,130), soilm.nlevel=10, lognorm=TRUE,
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
  ))
{
  r <- lulik(    ###CORRECT DEPTH and TOTAL MICROSITE HERE
    x=c(par,z),formula=formula,data=data,nlevels=nlevels,
    R=R, O2airfrac=O2airfrac, CH4airfrac = CH4airfrac,
    BD=BD, PD=PD,  Porosity=Porosity,
    p=p, RQ=RQ,  Depth_O2=Depth_O2, Depth_CH4 = Depth_CH4,Depth = Depth,
    microsite=microsite,total.microsite=total.microsite,time.step=time.step,
    flagO2=flagO2,  flagCH4=flagCH4, 
    iO2.lst=iO2.lst,iCH4.lst=iCH4.lst,
    reltol=reltol,verbose=verbose,
    sxtot.range=sxtot.range, nquantile=nquantile, mult=mult,
    soilm.cv=soilm.cv, soilm.range=soilm.range, soilm.nlevel=soilm.nlevel, 
    lognorm=lognorm,
    nu=nu,R0=R0,
    hyper=hyper)
  
  resid <- r$data[,3] - r$fit_Rh
  resid <- resid[!is.na(resid)]
  as.vector(crossprod(resid))
}

lulik.ch4 <- function(    ###CORRECT DEPTH and TOTAL MICROSITE HERE
  par,z,formula,data,nlevels,
  R, O2airfrac, CH4airfrac, BD, PD,  Porosity,
  p, RQ,  Depth_O2, Depth_CH4,Depth,
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=TRUE, ## whether run O2 and CH4 to steady states
  iO2.lst=NULL,iCH4.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,130), soilm.nlevel=10, lognorm=TRUE,
  nu=2,R0=diag(2),
  hyper=list(    ###CORRECT DEPTH and TOTAL MICROSITE HERE
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
  ))
{
  r <- lulik(    ###CORRECT DEPTH and TOTAL MICROSITE HERE
    x=c(z,par),formula=formula,data=data,nlevels=nlevels,
    R=R, O2airfrac=O2airfrac, CH4airfrac = CH4airfrac,
    BD=BD, PD=PD,  Porosity=Porosity,
    p=p, RQ=RQ,  Depth_O2=Depth_O2, Depth_CH4 = Depth_CH4, Depth = Depth,
    microsite=microsite,total.microsite=total.microsite,time.step=time.step,
    flagO2=flagO2,  flagCH4=flagCH4, 
    iO2.lst=iO2.lst,iCH4.lst=iCH4.lst,
    reltol=reltol,verbose=verbose,
    sxtot.range=sxtot.range, nquantile=nquantile, mult=mult,
    soilm.cv=soilm.cv, soilm.range=soilm.range, soilm.nlevel=soilm.nlevel, 
    lognorm=lognorm,
    nu=nu,R0=R0,
    hyper=hyper)
  
  resid <- r$data[,4] - r$fit_CH4 # this indicated it is 4th in l.mat
  resid[!is.na(resid)]
}

lulik.ch4.optim <- function(   ###CORRECT DEPTH and TOTAL MICROSITE HERE
  par,z,formula,data,nlevels,
  R, O2airfrac, CH4airfrac,  BD, PD,  Porosity,
  p, RQ,  Depth_O2, Depth_CH4,Depth,
  microsite,total.microsite,time.step,
  flagO2=TRUE,  flagCH4=TRUE, ## whether run O2 and CH4 to steady states
  iO2.lst=NULL,iCH4.lst=NULL,
  reltol=sqrt(.Machine$double.eps),verbose=0,
  sxtot.range=c(0.01,0.15), nquantile=10, mult=1/4,
  soilm.cv=20, soilm.range=c(70,130), soilm.nlevel=10, lognorm=TRUE,
  nu=2,R0=diag(2),
  hyper=list( ###CORRECT DEPTH and TOTAL MICROSITE HERE
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
  ))
{
  r <- lulik(   ###CORRECT DEPTH and TOTAL MICROSITE HERE
    x=c(z,par),formula=formula,data=data,nlevels=nlevels,
    R=R, O2airfrac=O2airfrac, CH4airfrac = CH4airfrac,
    BD=BD, PD=PD,  Porosity=Porosity,
    p=p, RQ=RQ,  Depth_O2=Depth_O2, Depth_CH4 = Depth_CH4,Depth = Depth,
    microsite=microsite,total.microsite=total.microsite,time.step=time.step,
    flagO2=flagO2,  flagCH4=flagCH4, 
    iO2.lst=iO2.lst,iCH4.lst=iCH4.lst,
    reltol=reltol,verbose=verbose,
    sxtot.range=sxtot.range, nquantile=nquantile, mult=mult,
    soilm.cv=soilm.cv, soilm.range=soilm.range, soilm.nlevel=soilm.nlevel, 
    lognorm=lognorm,
    nu=nu,R0=R0,
    hyper=hyper)
  
  resid <- r$data[,4] - r$fit_CH4 
  resid <- resid[!is.na(resid)]
  as.vector(crossprod(resid))
}

test <- function() #Error in if (log) { : argument is of length zero
{
  rm(list=ls())
  ## Session > Work Dir > Source File Loc
  trench <- read.csv("2015-16data.csv",head=T)
  ## Initial values to start estimation
  z0 <-   c(
    KmO2_Rh=245,alphaVmaxSx=230000000000, EaVmaxSx=72,
    KMSx=1000,
    alphaVmaxCH4prod=3000000,  EaVmaxCH4prod=100,
    KMSx_CH4prod=350,
    alphaVmaxCH4ox=7, EaVmaxCH4ox=30,
    KMCH4ox=1.00E-02,  KmO2_CH4ox=43,kl_CH4prod=3
  )
  
  hyper=list(
    KmO2_Rh=c(2.45,24500,1),
    alphaVmaxSx=c(2300000000,23000000000000,1),
    EaVmaxSx=c(7.2,7200,1),
    KMSx=c(10,100000,1),
    alphaVmaxCH4prod=c(30000,300000000,1),
    EaVmaxCH4prod=c(0.3,3000,1),
    KMSx_CH4prod=c(3.5,35000,1),
    alphaVmaxCH4ox=c(0.07,700,1),
    EaVmaxCH4ox=c(10,50,0),KMCH4ox=c(1.00E-04,1.00,1),
    KmO2_CH4ox=c(0.043,4300,1),
    kl_CH4prod=c(0.03,300,1)
  )
  
  source("damm_11a-non-steady-state.R")
  source("damm_utils_3.R")
  x0 <-   damm_btran(z0,hyper)

  source("lulik_4b.R")
  source("lulik_utils_2b.R")
  r0 <- lulik.rh(
    par=x0[1:4],z=x0[5:12],formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,
    microsite = 0,total.microsite = 10000, time.step = 300
  )
  hist(r0)
  
  r1 <- lulik.ch4(
    par=x0[5:12],z=x0[1:4],formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,
    microsite = 0,total.microsite = 10000, time.step = 300
  )
  hist(r1)
}