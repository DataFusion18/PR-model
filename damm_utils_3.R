## damm_utils_3.R
## ver 3) add-in new parameters to optimize over
## utility function for damm

## original scale -> unrestricted
## x : variate in uniform scale
## z : transformed scale on real line
## l : lower limit of the uniform interval
## r : upper limit of the uniform interval
## log: whether to apply log to x 
logit_tran <- function(x,l,r,log=FALSE)
{
  if(log){
    x <- log(x)
    l <- log(l)
    r <- log(r)
  }
  log((x-l)/(r-x))
}
## unrestricted scale -> original
expit_tran <- function(z,l,r,log=FALSE)
{
  if(log){
    l <- log(l)
    r <- log(r)
  }
  l.z <- exp(z)
  output <- ifelse(l.z>0 & is.infinite(l.z),r,(l+r*l.z)/(1+l.z))
  if(log){
    output <- exp(output)
  }
  output
}

## helper function to transform real parameters back to range
damm_tran <- function( ###ADD DEPTH and TOTAL MICROSITE HERE
  x,hyper=list(
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
  ## x: unrestricted parameter on real line
  ## hyper: the range and whether or not to apply log transformation 
  ## value: parameter on restricted range
  names(x) <- NULL         ###ADD DEPTH and TOTAL MICROSITE HERE
  lKmO2_Rh <- x[1]
  lalphaVmaxSx <- x[2]
  lEaVmaxSx <- x[3]
  lKMSx <- x[4]
  lDepth <- x[5]
  ltotal.microsite <- x[6]
  lalphaVmaxCH4prod <- x[7]
  lEaVmaxCH4prod <- x[8]
  lKMSx_CH4prod <- x[9]
  lalphaVmaxCH4ox <- x[10]
  lEaVmaxCH4ox <- x[11]
  lKMCH4ox <- x[12]
  lKmO2_CH4ox <- x[13]
  lkl_CH4prod <- x[14]
  
  ## KmO2_Rh ###ADD DEPTH and TOTAL MICROSITE HERE
  tmp <- hyper[["KmO2_Rh"]]
  KmO2_Rh <- expit_tran(lKmO2_Rh,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxSx
  tmp <- hyper[["alphaVmaxSx"]]
  alphaVmaxSx <- expit_tran(lalphaVmaxSx,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## EaVmaxSx
  tmp <- hyper[["EaVmaxSx"]]
  EaVmaxSx <- expit_tran(lEaVmaxSx,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMSx
  tmp <- hyper[["KMSx"]]
  KMSx <- expit_tran(lKMSx,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## Depth
  tmp <- hyper[["Depth"]]
  Depth <- expit_tran(lDepth,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## total.microsite
  tmp <- hyper[["total.microsite"]]
  total.microsite <- expit_tran(ltotal.microsite,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxCH4prod
  tmp <- hyper[["alphaVmaxCH4prod"]]
  alphaVmaxCH4prod <- expit_tran(lalphaVmaxCH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxCH4prod
  tmp <- hyper[["EaVmaxCH4prod"]]
  EaVmaxCH4prod <- expit_tran(lEaVmaxCH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMSx_CH4prod
  tmp <- hyper[["KMSx_CH4prod"]]
  KMSx_CH4prod <- expit_tran(lKMSx_CH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxCH4ox
  tmp <- hyper[["alphaVmaxCH4ox"]]
  alphaVmaxCH4ox <- expit_tran(lalphaVmaxCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxCH4ox
  tmp <- hyper[["EaVmaxCH4ox"]]
  EaVmaxCH4ox <- expit_tran(lEaVmaxCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## KMCH4ox
  tmp <- hyper[["KMCH4ox"]]
  KMCH4ox <- expit_tran(lKMCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## KmO2_CH4ox
  tmp <- hyper[["KmO2_CH4ox"]]
  KmO2_CH4ox <- expit_tran(lKmO2_CH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## kl_CH4prod
  tmp <- hyper[["kl_CH4prod"]]
  kl_CH4prod <- expit_tran(lkl_CH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
 
   r <- c(KmO2_Rh=KmO2_Rh,alphaVmaxSx=alphaVmaxSx,EaVmaxSx=EaVmaxSx,KMSx=KMSx,
         Depth=Depth, total.microsite=total.microsite,
         alphaVmaxCH4prod=alphaVmaxCH4prod,EaVmaxCH4prod=EaVmaxCH4prod,
         KMSx_CH4prod=KMSx_CH4prod,alphaVmaxCH4ox=alphaVmaxCH4ox,
         EaVmaxCH4ox=EaVmaxCH4ox,KMCH4ox=KMCH4ox,
         KmO2_CH4ox=KmO2_CH4ox,kl_CH4prod=kl_CH4prod) ###ADD DEPTH and TOTAL MICROSITE HERE
  #names(r) <- names(x)
  r
}

## helper function to back transform restricted parameters to real
damm_btran <- function( ###ADD DEPTH and TOTAL MICROSITE HERE
  z,hyper=list(
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
  ## z: parameter on uniform range
  ## value: parameter on real line
  names(z) <- NULL
  
  KmO2_Rh <- z[1]   ###ADD DEPTH and TOTAL MICROSITE HERE
  alphaVmaxSx <- z[2]
  EaVmaxSx <- z[3]
  KMSx <- z[4]
  Depth <- z[5]
  total.microsite <- z[6]
  alphaVmaxCH4prod <- z[7]
  EaVmaxCH4prod <- z[8]
  KMSx_CH4prod <- z[9]
  alphaVmaxCH4ox <- z[10]
  EaVmaxCH4ox <- z[11]
  KMCH4ox <- z[12]
  KmO2_CH4ox <- z[13]
  kl_CH4prod <- z[14]
  
  ## KmO2_Rh       ###ADD DEPTH and TOTAL MICROSITE HERE
  tmp <- hyper[["KmO2_Rh"]]
  lKmO2_Rh <- logit_tran(KmO2_Rh,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxSx
  tmp <- hyper[["alphaVmaxSx"]]
  lalphaVmaxSx=logit_tran(alphaVmaxSx,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxSx
  tmp <- hyper[["EaVmaxSx"]]
  lEaVmaxSx=logit_tran(EaVmaxSx,tmp[1],tmp[2],log=as.logical(tmp[3])) 
  ## KMSx
  tmp <- hyper[["KMSx"]]
  lKMSx=logit_tran(KMSx,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## Depth
  tmp <- hyper[["Depth"]]
  lDepth=logit_tran(Depth,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## total.microsite
  tmp <- hyper[["total.microsite"]]
  ltotal.microsite=logit_tran(total.microsite,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxCH4prod
  tmp <- hyper[["alphaVmaxCH4prod"]]
  lalphaVmaxCH4prod <- logit_tran(alphaVmaxCH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxCH4prod
  tmp <- hyper[["EaVmaxCH4prod"]]
  lEaVmaxCH4prod=logit_tran(EaVmaxCH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMSx_CH4prod
  tmp <- hyper[["KMSx_CH4prod"]]
  lKMSx_CH4prod <- logit_tran(KMSx_CH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## alphaVmaxCH4ox
  tmp <- hyper[["alphaVmaxCH4ox"]]
  lalphaVmaxCH4ox <- logit_tran(alphaVmaxCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## EaVmaxCH4ox
  tmp <- hyper[["EaVmaxCH4ox"]]
  lEaVmaxCH4ox=logit_tran(EaVmaxCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KMCH4ox
  tmp <- hyper[["KMCH4ox"]]
  lKMCH4ox=logit_tran(KMCH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## KmO2_CH4ox
  tmp <- hyper[["KmO2_CH4ox"]]
  lKmO2_CH4ox <- logit_tran(KmO2_CH4ox,tmp[1],tmp[2],log=as.logical(tmp[3]))
  ## kl_CH4prod
  tmp <- hyper[["kl_CH4prod"]]
  lkl_CH4prod <- logit_tran(kl_CH4prod,tmp[1],tmp[2],log=as.logical(tmp[3]))
  
  r <- c(lKmO2_Rh=lKmO2_Rh,     ###ADD DEPTH and TOTAL MICROSITE HERE
         lalphaVmaxSx=lalphaVmaxSx,lEaVmaxSx=lEaVmaxSx,lKMSx=lKMSx,
         lDepth=lDepth,ltotal.microsite=ltotal.microsite,
         lalphaVmaxCH4prod=lalphaVmaxCH4prod,lEaVmaxCH4prod=lEaVmaxCH4prod,
         lKMSx_CH4prod=lKMSx_CH4prod,lalphaVmaxCH4ox=lalphaVmaxCH4ox,
         lEaVmaxCH4ox=lEaVmaxCH4ox,lKMCH4ox=lKMCH4ox,
         lKmO2_CH4ox=lKmO2_CH4ox,lkl_CH4prod=lkl_CH4prod)
  #names(r) <- names(z)
  r
}

#####Error in damm_tran(x1) : object 'total.microsite' not found
tran_test_08082017 <- function() 
{
  rm(list=ls())
  source("damm_utils_3.R")
  hyper=list(
    KmO2_Rh=c(0.01,0.21,0),
    alphaVmaxSx=c(1666666.66666667,166666666.666667,1),
    EaVmaxSx=c(30,100,0),
    KMSx=c(8.00E-08,8.00E-06,1),
    Depth=c(5,30,0),
    total.microsite=c(1000,100000,0),
    alphaVmaxCH4prod=c(1722222222222.22,172222222222222,1),
    EaVmaxCH4prod=c(15,75,0),
    KMSx_CH4prod=c(1666666.66666667,166666666.666667,1),
    alphaVmaxCH4ox=c(8.32660531276481,832.660531276481,1),
    EaVmaxCH4ox=c(10,50,1),KMCH4ox=c(5,50,1),
    KmO2_CH4ox=c(0.01,0.21,1),
    kl_CH4prod=c(0.001,0.01,1)
  )
  x1=c(
    KmO2_Rh=1.417843,
    alphaVmaxSx=15.1561754430518,
    EaVmaxSx=0.224663729,
    KMSx=-15.46972856,
    Depth=10,
    total.microsite=10000,
    alphaVmaxCH4prod=15.56470445,
    EaVmaxCH4prod=1.086290931,
    KMSx_CH4prod=-15.46972856,
    alphaVmaxCH4ox=5.130463522,
    EaVmaxCH4ox=-0.025674405,
    KMCH4ox=-2.77212177,
    KmO2_CH4ox=0,
    kl_CH4prod=0
  )
  
  damm_tran(x1)
  x0 <-   c(
    lKmO2_Rh=logit_tran(0.1,0.01,0.21),
    lalphaVmaxSx=logit_tran(16666666.6666667,hyper[["alphaVmaxSx"]][1],hyper[["alphaVmaxSx"]][2],log=TRUE), 
    lEaVmaxSx=logit_tran(68,hyper[["EaVmaxSx"]][1],hyper[["EaVmaxSx"]][2]),
    lKMSx=logit_tran(0.0000008,hyper[["KMSx"]][1],hyper[["KMSx"]][2],log=TRUE),
    lDepth=logit_tran(10,hyper[["Depth"]][1],hyper[["Depth"]][2],log=TRUE),
    ltotal.microsite=logit_tran(10000,hyper[["total.microsite"]][1],hyper[["total.microsite"]][2],log=TRUE),
    lalphaVmaxCH4prod=logit_tran((6.2*(10^16)*(1/3600)),hyper[["alphaVmaxCH4prod"]][1],hyper[["alphaVmaxCH4prod"]][2],log=TRUE), 
    lEaVmaxCH4prod=logit_tran(50,hyper[["EaVmaxCH4prod"]][1],hyper[["EaVmaxCH4prod"]][2]),
    lKMSx_CH4prod=logit_tran(2*16666666.6666667,1666666.66666667,166666666.666667,log=TRUE),
    lalphaVmaxCH4ox=logit_tran(299757.791259533*(1/3600),hyper[["alphaVmaxCH4ox"]][1],hyper[["alphaVmaxCH4ox"]][2],log=TRUE), 
    lEaVmaxCH4ox=logit_tran(30,hyper[["EaVmaxCH4ox"]][1],hyper[["EaVmaxCH4ox"]][2]),
    lKMCH4ox=logit_tran(7,hyper[["KMCH4ox"]][1],hyper[["KMCH4ox"]][2]),
    lKmO2_CH4ox=logit_tran(0.10,0.01,0.21),
    lkl_CH4prod=logit_tran(0.005,0.001,0.01,log=TRUE)
  )
  options(digits=2,warn=1)
  (z <- damm_tran(x0))
  (x1 <- damm_btran(z))
  x0-x1
  z <- c(KmO2_Rh=0.10,
         alphaVmaxSx=2666666.66666667,
         EaVmaxSx=70,
         KMSx=5.00E-06,
         Depth=10,
         total.microsite=10000,
         alphaVmaxCH4prod=(6.2*(10^16)*(1/3600)),
         EaVmaxCH4prod=25,
         KMSx_CH4prod=2*16666666.6666667,
         alphaVmaxCH4ox=83.2660531276481,
         EaVmaxCH4ox=30,
         KMCH4ox=7,
         KmO2_CH4ox=0.10,
         kl_CH4prod=0.005)
  (x <- damm_btran(z))
  (z2 <- damm_tran(x))
  z2-z
}