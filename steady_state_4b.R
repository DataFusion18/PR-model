## steady_state_4.R
## version 4)
## helper script to optimize parameters of the liquid steady state model 
## using observed time series
## of SoilM, SoilT, Rh and CH4 data

## Arguments:
## z0 : initial values to start MCMC
## hyper: initial range and whether on log-scale for each parameter
## trench.file: excel data file name
## microsite.file: excel data file for fixed microsite
##                 0 to generate bivariate distribution automatically
## num.node      : number of cores for parallel computing,
##                default total number of cores minus one
## nsweep        : number of sweeps through the parameter space
## outfile       : root name for output files ending with
## tol           : relative tolerence for each conditional optimization
## ...           : additional argument to lulik
## value
##          : a time series plot of observed and fitted Rh and CH4 values
## param.csv: for final parameter values
## main.csv:  for data point and microsite specific output
## Rh.csv:   for data point specific Rh value
## CH4.csv:  for data point specific CH4 value
#######################################################

## load required packages
require(truncnorm)
require(doParallel)

## Load R functions
source("damm_11a-non-steady-state.R")
source("damm_utils_3.R")
source("lulik_4b.R")
source("lulik_utils_2b.R")
source("M3D-DAMM.R")
source("luprior_2.R")
source("lupost_2.R")

steady.state <- function(z0,hyper,trench.file,outfile,tol,
                         microsite.file=0,num.node=NULL,nsweep=999,
                         opt=TRUE,...)
{
  ## transformed initial parameters on the real line
  x0 <- damm_btran(z=z0,hyper=hyper)
  
  ## Load data
  trench <- read.csv(trench.file,head=T)
  
  ## prepare the high performance cluster
  if(is.null(num.node)){
    num.node <- detectCores() - 1 
  }
  cl <- makeCluster(num.node)
  registerDoParallel(cl)
  sink <- clusterEvalQ(cl,source("damm_11a-non-steady-state.R"))
  sink <- clusterEvalQ(cl,source("damm_utils_3.R"))
  
  ## Storage for the algorithm
  history <- matrix(NA,nsweep,length(x0))
  fnval <- vector("list",nsweep)
  O2.val <- vector("list",nsweep)
  CH4.val <- vector("list",nsweep)
  if(opt){
  ## start the iterations with optimizing O2 while fixing CH4
  nls.O2.0 <- optim( ####UPDATE NUMBERS HERE
    par=x0[1:6],fn=lulik.rh.optim,gr=NULL,method="Nelder-Mead",control=list(reltol=tol),
    z=x0[7:14],formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth=10, microsite = 0,
    total.microsite = 10000, time.step = 300, ####REMOVE total.microsite and Depth
    iO2.lst=NULL,iCH4.lst=NULL,hyper=hyper,...
  )
  ## evaluate loglik at the last values
  nls.O2.r0 <- lulik(  ####UPDATE NUMBERS HERE
    x=c(nls.O2.0$par,x0[7:14]),formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth=10,microsite = 0,total.microsite = 10000, time.step = 300,
    flagO2=TRUE, flagCH4=FALSE,hyper=hyper,...) ####REMOVE total.microsite and Depth
  nls.O2.x0 <- c(nls.O2.0$par,x0[7:14])
  
  ## Note the initial O2 was not updated successively because minpack.lm does not allow
  ## that, only after convergence
  
  for(i in 1:nsweep){
    ## continue iteration with CH4 while fixing O2
    cat("iter ",i," CH4 given O2.\n")
    #browser()
    nls.CH4 <- optim(  ####UPDATE NUMBERS HERE
      par=nls.O2.x0[7:14],fn=lulik.ch4.optim,gr=NULL,method="Nelder-Mead",control=list(reltol=tol),
      z=nls.O2.x0[1:6],formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
      R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
      BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
      p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth=10,
      microsite = 0,total.microsite = 10000, time.step = 300, ####REMOVE total.microsite and Depth
      iO2.lst=nls.O2.r0$O2.lst,iCH4.lst=nls.O2.r0$CH4.lst,
      hyper=hyper,...)
    ## evaluate loglik at the last values
    nls.CH4.r1 <- lulik(  ####UPDATE NUMBERS HERE
      x=c(nls.O2.x0[1:6],nls.CH4$par),formula=~SoilM+SoilT+Rh+CH4, #UPDATE NUMBERS 1:4 HERE
      data=trench,nlevels=100,
      R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
      BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
      p=0.00046,RQ=1, Depth_O2=10, Depth_CH4 = 10,Depth = 10,
      microsite = 0,total.microsite = 10000, time.step = 300, ####REMOVE total.microsite and Depth
      iO2.lst=nls.O2.r0$O2.lst,iCH4.lst=nls.O2.r0$CH4.lst,
      hyper=hyper,...)
    nls.CH4.x1 <- c(nls.O2.x0[1:6],nls.CH4$par) #UPDATE NUMBERS 1:4 HERE
    
    ## store results
    fnval[[i]] <- nls.CH4
    (history[i,] <- nls.CH4.r1$tpar) ## original scale parameters
    O2.val[[i]] <- nls.CH4.r1$O2.lst
    CH4.val[[i]] <- nls.CH4.r1$CH4.lst
    cat(history[i,],", convergence=",nls.CH4$convergence,", fn=",nls.CH4$value,"\n",
        file=paste(outfile,"_obj.txt",sep=""),append=TRUE)
    save(O2.val,CH4.val,history,i,file=paste(outfile,"_dump.rda",sep=""))
    
    ## repeat iteration with O2 while fixing CH4
    cat("iter ",i," O2 given CH4.\n")
    nls.O2.1 <- optim(   ####UPDATE NUMBERS HERE
      par=nls.CH4.x1[1:6],fn=lulik.rh.optim,gr=NULL,method="Nelder-Mead",control=list(reltol=tol),
      z=nls.CH4.x1[7:14], formula=~SoilM+SoilT+Rh+CH4, #UPDATE NUMBERS 1:4 HERE
      data=trench,nlevels=100,
      R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
      BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
      p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth=10,
      microsite = 0,total.microsite = 10000, time.step = 300, #UPDATE NUMBERS 1:4 HERE
      iO2.lst=nls.CH4.r1$O2.lst,iCH4.lst=nls.CH4.r1$CH4.lst,
      hyper=hyper,...
    )
    ## evaluate the latest likelihood function
    nls.O2.r1 <- lulik(
      x=c(nls.O2.1$par,nls.CH4.x1[7:14]), #UPDATE NUMBERS 1:4 HERE
      formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
      R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8,
      BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
      p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth=10,
      microsite = 0,total.microsite = 10000, time.step = 300, ####REMOVE total.microsite and Depth
      flagO2=TRUE, flagCH4=FALSE,
      iO2.lst=nls.CH4.r1$O2.lst,iCH4.lst=nls.CH4.r1$CH4.lst,
      hyper=hyper,...)
    nls.O2.x1 <- c(nls.O2.1$par,nls.CH4.x1[7:14]) #UPDATE NUMBERS 1:4 HERE
    
    ## stop condition
    delta <- abs(nls.O2.x1-nls.O2.x0)
    cat(delta,"\n")
    if(all(delta < tol * abs(nls.O2.x0))){
      cat("Converged a iteration ",i,"\n")
      break
    }
    
    ## update
    nls.O2.x0 <- nls.O2.x1
  }
  }
  else{
    nls.O2.x1 <- x0
  }
  
  ## prepare output
  x1 <- nls.O2.x1
  output <- lulik(
    x=x1, formula=~SoilM+SoilT+Rh+CH4,data=trench,nlevels=100,
    R=0.008314472, O2airfrac=0.209, CH4airfrac = 1.8,
    BD=1.009, PD=2.52,Porosity=(1-0.13)*0.6,
    p=0.00046,RQ=1, Depth_O2=10, Depth_CH4=10,Depth=10,
    microsite = 0,total.microsite = 10000, time.step = 300, ####REMOVE total.microsite and Depth
    hyper=hyper,...)
  save(fnval,history,O2.val,CH4.val,output, file=paste(outfile,".rda",sep=""))  ## save compressed R data
  write.csv(output$tpar, file=paste(outfile,"_parm.csv",sep=""),row.names = FALSE)
  write.csv(output$output,file=paste(outfile,"_main.csv",sep=""),row.names = FALSE)
  write.csv(output$fit_Rh, file=paste(outfile,"_Rh.csv",sep=""),row.names = FALSE)
  write.csv(output$fit_CH4, file=paste(outfile,"_CH4.csv",sep=""),row.names = FALSE)
  
  ## fit versus observed data #FOR PLOTTING
  png(file=paste(outfile,".png",sep=""),width=3000,height=3000,res=400)
  op <- par(mar=c(4.1,4.1,1.1,1.1),mfrow=c(2,1)) 
  with(output,{
    plot(data[,"Rh__"],bg="steelblue1",col=1,pch=21,cex=2,
         xlab="Time",ylab="Rh",main="",
         ylim=range(data[,"Rh__"],fit_Rh,na.rm=T))
    lines(fit_Rh,lty=2)
    plot(data[,"CH4__"],bg="salmon",col=1,pch=21,cex=2,
         xlab="Time",ylab="CH4",main="",
         ylim=range(data[,"CH4__"],fit_CH4,na.rm=T))
    lines(fit_CH4,lty=2)
  })
  par(op)
  dev.off()
  ## stop cluster
  stopCluster(cl)
  
  output
}

test <- function()
{
  ## in model_framework_14_a.R
}
