## model_framework_14.R
## ver 14) parameter estimation
## sub ver a) Liquid steady state model
rm(list=ls())
#### Note the following packages are required to run this bayesian analysis
## truncnorm, doParallel
source("steady_state_4b.R")

## initial range and whether on log-scale for each parameter
z0 <-   c(
  KmO2_Rh=245,alphaVmaxSx=230000000000, EaVmaxSx=72,
  KMSx=1000, Depth=10, total.microsite=10000,
  alphaVmaxCH4prod=3000000,  EaVmaxCH4prod=75,
  KMSx_CH4prod=350,
  alphaVmaxCH4ox=0.07, EaVmaxCH4ox=30,
  KMCH4ox=1.00E-02,  KmO2_CH4ox=43,kl_CH4prod=3
)

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
)

## parameter estimation starts
est0 <- steady.state(
  z0=z0,hyper=hyper,trench.file="flux.csv",
  outfile="out",tol=5e-4,opt=TRUE)
