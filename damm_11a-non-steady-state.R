#Non-steady state model for CH4-Liq

source("soilm_pdf_3.R") ## load script for bi-variate distribution

## Arguments
damm <- function( 
  R, O2airfrac, CH4airfrac, BD, PD,  Porosity,
  p, RQ, 
  SoilM.o,   SoilT ,  Depth_O2, Depth_CH4, Depth,   ###correct for diff depth
  microsite,total.microsite,
  time.step,
  iO2, iCH4, ## initial O2 and CH4 values
  flagO2,  flagCH4, ## whether run O2 and CH4 to steady states
  KmO2_Rh, alphaVmaxSx,   EaVmaxSx,   KMSx,   
  alphaVmaxCH4prod, EaVmaxCH4prod,  KMSx_CH4prod,
  alphaVmaxCH4ox, EaVmaxCH4ox,   KMCH4ox, KmO2_CH4ox, kl_CH4prod,
  reltol,verbose,
  sxtot.range=NULL, nquantile=NULL, mult=NULL,
  soilm.cv=NULL, soilm.range=NULL, soilm.nlevel=NULL, lognorm=NULL
)
{
  ## Build microsite bi-variate distribution
  if(microsite==0){
    ## Check site specific limit
    if(SoilM.o*soilm.range[2]/100 > (Porosity*100)){
      soilm.range[2] <- (Porosity*100)*95/SoilM.o
      if(verbose>3){
        cat("Soil Moisture upper limit updated to ",soilm.range[2],"\n")
      }
    }

    soilm1 <- soilm.microsite(
      soilm=SoilM.o,soilm.cv=soilm.cv,soilm.range = soilm.range,
      soilm.nlevel = soilm.nlevel,lognorm=lognorm)
    sxtot1 <- sxtot.microsite(
      sxtot.range = sxtot.range,nquantile = nquantile,nsite=total.microsite,
      mult=mult,verbose = FALSE)
    microsite <- soilm.pdf(sxtot1,soilm1)
  }
  #
  #
  Sxtot <- microsite$sxtot
  SoilM <- microsite$soilm
  Freq <- microsite$freq

  ## initial values
  if(length(iO2)==1){
    iO2 <- rep(iO2,length(Sxtot))
  }
  if(length(iCH4)==1){
    iCH4 <- rep(iCH4,length(Sxtot))
  }

  #### Need to check for possible changes to be done
  ## per micro-site steady state outputs and interim results
  vRh <- rep(NA,length(Sxtot))
  O2 <- rep(NA,length(Sxtot))
  vCH4 <- rep(NA,length(Sxtot))
  Rh_calc <- rep(NA,length(Sxtot))
  O2_calc <- rep(NA,length(Sxtot))
  CH4_diff <- rep(NA,length(Sxtot))
  CH4_calc <- rep(NA,length(Sxtot))
  a_calc <- rep(NA,length(Sxtot)) #need to remove
  flag <- rep(NA,length(Sxtot)) #need to remove
  output <- matrix(NA,length(Sxtot),18) ###CORRECT HERE to 18 # 16 without D_P and flag

  #Need to change here
  colnames(output) <- c("O2","O2_l_micromolperliter","O2_l_micromolpermicrosite","Rh_l_micromolpermicrosite","Rh",
                        "O2_l_micromolpermicrosite_conc","O2_l_micromolperliter_conc","CH4","CH4_l_micromolperliter",
                        "CH4_l_micromolpermicrosite","CH4prod","CH4_l_conc_micromolpermicrosite",
                        "CH4_l_conc_micromolperliter","CH4ox","CH4Flux_l_micromolpermicrosite",
                     "CH4Flux_l_micromolperliter","D_P","flag")
                       
    for(i in 1:length(Sxtot)){
    damm_return <- damm_work(
      R=R, O2airfrac=O2airfrac,CH4airfrac=CH4airfrac,BD=BD,PD=PD,
      Porosity=Porosity,p=p,RQ=RQ,
      SoilM.o=SoilM.o,SoilT=SoilT,Depth_O2=Depth_O2,Depth_CH4=Depth_CH4,Depth=Depth, ###correct for diff depth
      Sxtot=Sxtot[i], SoilM=SoilM[i],
      total.microsite=total.microsite,
      time.step=time.step,
      iO2=iO2[i], iCH4=iCH4[i], ## initial O2 and CH4 values
      flagO2=flagO2,  flagCH4=flagCH4, ## whether run O2 and CH4 to steady states
      KmO2_Rh=KmO2_Rh, alphaVmaxSx=alphaVmaxSx,   EaVmaxSx=EaVmaxSx,
      KMSx=KMSx,alphaVmaxCH4prod=alphaVmaxCH4prod,
      EaVmaxCH4prod=EaVmaxCH4prod,KMSx_CH4prod=KMSx_CH4prod,
      alphaVmaxCH4ox=alphaVmaxCH4ox, EaVmaxCH4ox=EaVmaxCH4ox,
      KMCH4ox=KMCH4ox, KmO2_CH4ox=KmO2_CH4ox, kl_CH4prod=kl_CH4prod,
      reltol=reltol,verbose=verbose)
    output[i,] <- damm_return
    vRh[i] <- damm_return[5] #5
    O2[i] <- damm_return[1] #1
    vCH4[i] <- damm_return[16] #16 for CH4Flux_l_micromolperliter
    a_calc[i] <- damm_return[17] #17
    flag[i] <- damm_return[18] # 18
    Rh_calc[i] <- vRh[i]*Freq[i]
    CH4_calc[i] <- vCH4[i] * Freq[i]
    if(verbose>0) {
      cat(" micro site ",i," Rh=",vRh[i],"O2=",O2[i],"CH4=",vCH4[i]," flag=",flag[i],"\n")
    }

  }
  microsite.vol <- (3.14*12.5*12.5*Depth)/(total.microsite)
  Henry_constant_CO2 <- 3.4*10^(-2)*exp(2400*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  Henry_constant_CH4 <- 1.4*10^(-3)*exp(1700*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  
  damm_final_return <- rep(NA,2)

  damm_final_return[1] <- (sum(Rh_calc)/sum(Freq))*(1/Henry_constant_CO2)*(10^-6)
  damm_final_return[1] <- damm_final_return[1]*12*(1/1000)*1000*Depth*10000*(3600/time.step)

  damm_final_return[2] <- (sum(CH4_calc)/sum(Freq))*(1/Henry_constant_CH4)*(10^-6)
  damm_final_return[2] <- damm_final_return[2]*12*(1/1000)*1000*Depth*10000*(3600/time.step)
  
  list(value=damm_final_return,flag=flag,output=output,microsite = microsite)
}

test2 <- function()
{
  rm(list=ls())
  source("damm_11a-non-steady-state.R")
  microsite <- read.csv("microsite.csv")
  sink("./damm_liq_debug_1.txt")
  r <- damm( 
    R=0.008314472, O2airfrac=0.209, CH4airfrac=1.8, BD=1.009, PD=2.52,
    Porosity=(1-0.13)*0.6, p=0.00046,RQ=1,
    SoilM.o=30,   SoilT=20,   Depth_O2=10,Depth_CH4=10,Depth=10,  
    microsite=microsite,  total.microsite = 10000, time.step = 3600,
    
    iO2=rep(0.209,nrow(microsite)), 
    iCH4=rep(1.8,nrow(microsite)),
    flagO2=TRUE, flagCH4=TRUE,
    
    KmO2_Rh=245, alphaVmaxSx=230000000000,   EaVmaxSx=72,   KMSx=1000,   
    alphaVmaxCH4prod=3000000, EaVmaxCH4prod=100,  KMSx_CH4prod=350,
    alphaVmaxCH4ox=0.07, EaVmaxCH4ox=30,   KMCH4ox=1.00E-02, KmO2_CH4ox=43, kl_CH4prod=3,
    
    reltol=sqrt(.Machine$double.eps),verbose=99
  )
  sink()
  r$value
}

damm_work <- function( 
  R, O2airfrac,  CH4airfrac, BD, PD,  Porosity,
  p, RQ, 
  SoilM.o, SoilM, SoilT,  Sxtot, Depth_O2,Depth_CH4,Depth,  ###correct for diff depth
  total.microsite,  time.step,
  iO2, iCH4, ## initial O2 and CH4 values
  flagO2,  flagCH4, ## whether run O2 and CH4 to steady states
  KmO2_Rh, alphaVmaxSx,   EaVmaxSx,   KMSx,   
  alphaVmaxCH4prod, EaVmaxCH4prod,  KMSx_CH4prod,
  alphaVmaxCH4ox, EaVmaxCH4ox,   KMCH4ox, KmO2_CH4ox, kl_CH4prod,
  reltol,verbose
)
{
  if(verbose>3){
    cat("Soil M observed=",SoilM.o,", Soil T=",SoilT,
        "Sxtot =",Sxtot, ",Soil M =",SoilM,"\n")
  }
  if(verbose>3) {
    cat("KmO2_Rh=",KmO2_Rh," ",
        "alphaVmaxSx=",alphaVmaxSx, " ", 
        "EaVmaxSx=",EaVmaxSx, " ", 
        "KMSx=",KMSx, " ",
        "Depth=", Depth, " ", ###Include Depth_O2 and Depth_CH4
        "total.microsite=", total.microsite, " ",
        "alphaVmaxCH4prod=",alphaVmaxCH4prod," ",
        "EaVmaxCH4prod=", EaVmaxCH4prod, " ",
        "KMSx_CH4prod", KMSx_CH4prod, " ",
        "alphaVmaxCH4ox=",alphaVmaxCH4ox," ", 
        "EaVmaxCH4ox=",EaVmaxCH4ox," ",
        "KMCH4ox=",KMCH4ox," ",
        "KmO2_CH4ox", KmO2_CH4ox, " ",
        "kl_CH4prod", kl_CH4prod, " ", 
        "\n") 
  }
  
  #browser()
  ## derived constants 
  Dliq <- 1/Porosity^3 
  ## derived parameters
  VmaxCH4prod <- (alphaVmaxCH4prod)*exp(-EaVmaxCH4prod/(R*(SoilT+273.15))) 
  VmaxCH4ox <- (alphaVmaxCH4ox)*exp(-EaVmaxCH4ox/(R*(SoilT+273.15))) 
  
  microsite.vol <- (3.14*12.5*12.5*Depth)/(total.microsite)
  microsite.area <- (3.14*12.5*12.5)/(total.microsite)
  
  Henry_constant_O2 <- 1.3*10^(-3)*exp(1700*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  Henry_constant_CO2 <- 3.4*10^(-2)*exp(2400*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  Henry_constant_CH4 <- 1.4*10^(-3)*exp(1700*((1/(273.15+SoilT))-1/298.15)) ## derived constant and parameters
  
  ## derived constant and parameters
  D_P <-(2*0.261^3+0.04*0.261)*((Porosity-SoilM/100)/0.261)^(2+(3/2.9))
  #D_P <- ((Porosity-SoilM/100)^(4.0/3.0))*(((SoilT+273.15)/293.15)^1.75)##a4_3 #update temp sensitivity in steadystate model also
  a4_3_max <- ((Porosity-0/100)^(4/3))*(((20+273.15)/293.15)^1.75)
  Dgas <- 1/a4_3_max
  DCH4 <- Dgas
  DO2 <- Dgas  

  kl_CH4prod <- kl_CH4prod
  kMO2_Rh <- KmO2_Rh
  KMCH4ox <- KMCH4ox
  kMO2_CH4ox <- KmO2_CH4ox
  
  Sxtot_micromolperliter <- Sxtot*p*(100/SoilM)*1000*(10^6)*(1/12)
  Sx <- Sxtot_micromolperliter*Dliq*(SoilM/100)^(3.0) ## derived constant and parameters
  MMSx <- Sx/(KMSx+Sx)
  VmaxSx <- (alphaVmaxSx)*exp(-EaVmaxSx/(R*(SoilT+273.15))) 
  
  if(verbose>3) {
    cat("Dliq=",Dliq," ",
        "VmaxCH4prod=",VmaxCH4prod," ",
        "VmaxCH4ox=",VmaxCH4ox," ",
        "microsite.vol=",microsite.vol," ",
        "microsite.area=",microsite.area," ",
        "D_P=",D_P," ",
        "DO2=",DO2," ",
        "DCH4=",DCH4," ",
        "kl_CH4prod=",kl_CH4prod," ",
        "kMO2_Rh=",kMO2_Rh," ",
        "KMCH4ox=",KMCH4ox," ", 
        "kMO2_CH4ox=",kMO2_CH4ox," ",
        "Sx=",Sx," ",
        "MMSx=",MMSx," ",
        "VmaxSx=",VmaxSx,"\n") 
  }
  
  ## Initialize
  #Need to change here
  iO2 <- O2airfrac*DO2*D_P #(for initialization)
  iO2_l_micromolperliter <- iO2*Henry_constant_O2*(10^6)
  iO2_l_micromolpermicrosite <- iO2_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
  iRh_l_micromolpermicrosite <- VmaxSx*MMSx*(iO2_l_micromolperliter/(KmO2_Rh+iO2_l_micromolperliter))*(SoilM/100)*(1/1000)*(microsite.vol)*(time.step)
  iRh <- iRh_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol)
  iO2_l_micromolpermicrosite_conc <- iO2_l_micromolpermicrosite-(iRh_l_micromolpermicrosite*(1/RQ))
  iO2_l_micromolperliter_conc <- iO2_l_micromolpermicrosite_conc*1000*(100/SoilM)*(1/microsite.vol)

  #Need to change here
  iCH4 <- CH4airfrac*DCH4*D_P #(for initialization)
  iCH4_l_micromolperliter <- iCH4*Henry_constant_CH4
  iCH4_l_micromolpermicrosite <- iCH4_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
  iCH4prod <- ((VmaxCH4prod*(1/1000)*(SoilM/100)*microsite.vol)/(1+(iO2_l_micromolperliter_conc/kl_CH4prod)))*(Sx/(Sx+KMSx_CH4prod))*(time.step)
  iCH4_l_conc_micromolpermicrosite <- iCH4_l_micromolpermicrosite+iCH4prod 
  iCH4_l_conc_micromolperliter <- iCH4_l_conc_micromolpermicrosite*1000*(100/SoilM)*(1/microsite.vol)
  iCH4ox <- ((VmaxCH4ox*(1/1000)*(SoilM/100)*microsite.vol))*(iCH4_l_conc_micromolperliter/(iCH4_l_conc_micromolperliter+KMCH4ox))*(iO2_l_micromolperliter_conc/(iO2_l_micromolperliter_conc+KmO2_CH4ox))*time.step 
  iCH4Flux_l_micromolpermicrosite <- iCH4prod-iCH4ox
  iCH4Flux_l_micromolperliter <- iCH4Flux_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol) 
  
  deltaO2 <- 9999.0+reltol 
  deltaCH4 <- 9999.0+reltol 
  iter <- 0 
  
  #Need to change here
  if(verbose>0) cat("iter ",iter," O2=",iO2," O2_l_micromolperliter=",iO2_l_micromolperliter,
                    " O2_l_micromolpermicrosite=",iO2_l_micromolpermicrosite," 
                    Rh_l_micromolpermicrosite=",iRh_l_micromolpermicrosite," Rh=",iRh,
                    " O2_l_micromolpermicrosite_conc=",iO2_l_micromolpermicrosite_conc,
                    " O2_l_micromolperliter_conc=",iO2_l_micromolperliter_conc,
                    " CH4=",iCH4,"CH4_l_micromolperliter=",iCH4_l_micromolperliter,
                    "CH4_l_micromolpermicrosite=",iCH4_l_micromolpermicrosite," CH4prod=", iCH4prod,
                    " CH4_l_conc_micromolpermicrosite=", iCH4_l_conc_micromolpermicrosite,
                    " CH4_l_conc_micromolperliter=", iCH4_l_conc_micromolperliter,
                    " CH4ox=",iCH4ox,"CH4Flux_l_micromolpermicrosite=",iCH4Flux_l_micromolpermicrosite, 
                    "CH4Flux_l_micromolperliter=",iCH4Flux_l_micromolperliter, "\n") 
  
  max_iter__ <- 99999
  max_iter_flag <- 0
  
  if(flagO2){
    ## If update O2, run O2 iteration until steady state
    ## run simulation until steady state 
    #while( (deltaCH4 > reltol) || (deltaO2 > reltol) )
    while( (deltaO2 > reltol))
    {
      if(iter > max_iter__){
        max_iter_flag <- 1
        break
      }

      #Need to change here
      O2_n_tmp <- O2airfrac*DO2*D_P #this is just to make it identical to 1st iteration value
      O2 <- O2_n_tmp
      O2_l_micromolperliter <- O2*Henry_constant_O2*(10^6)
      O2_l_micromolpermicrosite <- O2_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
      Rh_l_micromolpermicrosite <- VmaxSx*MMSx*(O2_l_micromolperliter/(KmO2_Rh+O2_l_micromolperliter))*(SoilM/100)*(1/1000)*(microsite.vol)*(time.step)
      Rh <- Rh_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol)
      O2_l_micromolpermicrosite_conc <- O2_l_micromolpermicrosite-(Rh_l_micromolpermicrosite*(1/RQ))
      O2_l_micromolperliter_conc <- O2_l_micromolpermicrosite_conc*1000*(100/SoilM)*(1/microsite.vol)     
      
       ## Update delta
      deltaO2 <- abs(iO2-O2)/iO2 
      ## update states

      iO2 <- O2
      iO2_l_micromolperliter <- O2_l_micromolperliter
      iO2_l_micromolpermicrosite <- O2_l_micromolpermicrosite
      iRh_l_micromolpermicrosite <- Rh_l_micromolpermicrosite
      iRh <- Rh
      iO2_l_micromolpermicrosite_conc <- O2_l_micromolpermicrosite_conc 
      iO2_l_micromolperliter_conc <- O2_l_micromolperliter_conc

      iter <- iter + 1 
      if(verbose>2) cat("iter ",iter," O2=",iO2," O2_l_micromolperliter=",iO2_l_micromolperliter,
                        " O2_l_micromolpermicrosite=",iO2_l_micromolpermicrosite,   
                        " Rh_l_micromolpermicrosite=",iRh_l_micromolpermicrosite," Rh=",iRh,
                        " O2_l_micromolpermicrosite_conc=",iO2_l_micromolpermicrosite_conc,
                        " O2_l_micromolperliter_conc=",iO2_l_micromolperliter_conc,
                        "(delta=",deltaO2,")","\n") 
    } ## End O2 iterations
    iter <- 0  ## re-set iterations for CH4
  } ## End update O2
  else{
    ## keep initial values
    O2 <- iO2
    O2_l_micromolperliter <- iO2_l_micromolperliter
    O2_l_micromolpermicrosite <-  iO2_l_micromolpermicrosite
    Rh_l_micromolpermicrosite <- iRh_l_micromolpermicrosite
    Rh <- iRh
    O2_l_micromolpermicrosite_conc <- iO2_l_micromolpermicrosite_conc 
    O2_l_micromolperliter_conc <- iO2_l_micromolperliter_conc

  } ## End not update O2
  if(verbose>0) cat("###### End O2 Iterations ##### \n iter = ",
                    iter," O2=",iO2,"(delta=",deltaO2,")","\n") 
  
  if(flagCH4){
    ## If update CH4, run CH4 iteration until steady state
    ## run simulation until steady state 
    #while( (deltaCH4 > reltol) || (deltaO2 > reltol) )
    while( (deltaCH4 > reltol))
    {
      if(iter > max_iter__){
        max_iter_flag <- 1
        break
      }

      #Need to change here
      CH4 <- CH4airfrac*DCH4*D_P #this is just to it identical to 1st iteration value
      CH4_l_micromolperliter <- CH4*Henry_constant_CH4
      CH4_l_micromolpermicrosite <- CH4_l_micromolperliter*(SoilM/100)*(1/1000)*microsite.vol
      CH4prod <- ((VmaxCH4prod*(1/1000)*(SoilM/100)*microsite.vol)/(1+(O2_l_micromolperliter_conc/kl_CH4prod)))*(Sx/(Sx+KMSx_CH4prod))*(time.step)
      CH4_l_conc_micromolpermicrosite <- CH4_l_micromolpermicrosite+CH4prod 
      CH4_l_conc_micromolperliter <- CH4_l_conc_micromolpermicrosite*1000*(100/SoilM)*(1/microsite.vol)
      CH4ox <- ((VmaxCH4ox*(1/1000)*(SoilM/100)*microsite.vol))*(CH4_l_conc_micromolperliter/(CH4_l_conc_micromolperliter+KMCH4ox))*(O2_l_micromolperliter_conc/(O2_l_micromolperliter_conc+KmO2_CH4ox))*time.step 
      CH4Flux_l_micromolpermicrosite <- CH4prod-CH4ox
      CH4Flux_l_micromolperliter <- CH4Flux_l_micromolpermicrosite*(100/SoilM)*(10^3)*(1/microsite.vol) 
      

      ##### Need to change here to end
      
      ## Update delta
      deltaCH4 <- abs(iCH4-CH4)/iCH4 

      iCH4 <- CH4
      iCH4_l_micromolperliter <- CH4_l_micromolperliter
      iCH4_l_micromolpermicrosite <- CH4_l_micromolpermicrosite 
      iCH4prod <- CH4prod
      iCH4_l_conc_micromolpermicrosite <- CH4_l_conc_micromolpermicrosite
      iCH4_l_conc_micromolperliter <- CH4_l_conc_micromolperliter
      iCH4ox <- CH4ox
      iCH4Flux_l_micromolpermicrosite <- CH4Flux_l_micromolpermicrosite
      iCH4Flux_l_micromolperliter <- CH4Flux_l_micromolperliter
        
      iter <- iter + 1 
      if(verbose>2) cat("iter ",iter," CH4=",iCH4," CH4_l_micromolperliter=", iCH4_l_micromolperliter,
                        " CH4_l_micromolpermicrosite=", iCH4_l_micromolpermicrosite," CH4prod=", iCH4prod,
                        " CH4_l_conc_micromolpermicrosite=", iCH4_l_conc_micromolpermicrosite,
                        " CH4_l_conc_micromolperliter=", iCH4_l_conc_micromolperliter,
                        " CH4ox=",iCH4ox,"CH4Flux_l_micromolpermicrosite=",iCH4Flux_l_micromolpermicrosite,
                        "CH4Flux_l_micromolperliter=",iCH4Flux_l_micromolperliter,
                        "(deltaCH4=",deltaCH4,")\n") 
    } ## End iterations
  } ## End updating CH4
  else{
    ## keep initial values
    CH4 <- iCH4
    CH4_l_micromolperliter <- iCH4_l_micromolperliter
    CH4_l_micromolpermicrosite <- iCH4_l_micromolpermicrosite 
    CH4prod <- iCH4prod
    CH4_l_conc_micromolpermicrosite <- iCH4_l_conc_micromolpermicrosite
    CH4_l_conc_micromolperliter <- iCH4_l_conc_micromolperliter
    CH4ox <- iCH4ox
    CH4Flux_l_micromolpermicrosite <- iCH4Flux_l_micromolpermicrosite
    CH4Flux_l_micromolperliter <- iCH4Flux_l_micromolperliter    
    
  } ## End not updating CH4
  if(verbose>0) cat("###### End CH4 Iterations ##### \n iter = ",
                    iter," CH4=",iCH4,"(delta=",deltaCH4,")","\n") 
  ## Return steady state
  damm_return <- rep(NA,18)
 
  
  names(damm_return) <- c("O2","O2_l_micromolperliter","O2_l_micromolpermicrosite","Rh_l_micromolpermicrosite","Rh",
                        "O2_l_micromolpermicrosite_conc","O2_l_micromolperliter_conc","CH4","CH4_l_micromolperliter",
                      "CH4_l_micromolpermicrosite","CH4prod","CH4_l_conc_micromolpermicrosite",
                     "CH4_l_conc_micromolperliter","CH4ox","CH4Flux_l_micromolpermicrosite",
                      "CH4Flux_l_micromolperliter","D_P","flag")
  damm_return[1] <- O2 
  damm_return[2] <- O2_l_micromolperliter  
  damm_return[3] <- O2_l_micromolpermicrosite 
  damm_return[4] <- Rh_l_micromolpermicrosite
  damm_return[5] <- Rh
  damm_return[6] <- O2_l_micromolpermicrosite_conc
  damm_return[7] <- O2_l_micromolperliter_conc
  damm_return[8] <- CH4 
  damm_return[9] <- CH4_l_micromolperliter 
  damm_return[10] <- CH4_l_micromolpermicrosite 
  damm_return[11] <- CH4prod 
  damm_return[12] <- CH4_l_conc_micromolpermicrosite
  damm_return[13] <- CH4_l_conc_micromolperliter
  damm_return[14] <- CH4ox
  damm_return[15] <- CH4Flux_l_micromolpermicrosite
  damm_return[16] <- CH4Flux_l_micromolperliter
  damm_return[17] <- D_P
  damm_return[18] <- max_iter_flag
  return(damm_return) 
}
test1 <- function()
{
  rm(list=ls()) 
  source("damm_11a-non-steady-state.R")
  microsite <- read.csv("microsite.csv")
  sink("./damm_liq_debug.txt")
  r <- damm_work( 
    R=0.0083, O2airfrac=0.209,  CH4airfrac=1.8,BD=1.009, PD=2.52,  
    Porosity=0.5220,p=0.00046, RQ=1, 
    SoilM.o=30, SoilM=22.5, SoilT=20, Sxtot=0.0010, 
    Depth_O2=10,Depth_CH4=10,Depth=10,  ###correct for diff depth
    total.microsite=10000, time.step=3600,
    iO2=0.209, iCH4=1.8, ## initial O2 and CH4 values
    flagO2=TRUE,  flagCH4=TRUE, ## whether run O2 and CH4 to steady states
    
    KmO2_Rh=245, alphaVmaxSx=230000000000,   EaVmaxSx=72,   KMSx=1000,   
    alphaVmaxCH4prod=3000000, EaVmaxCH4prod=100,  KMSx_CH4prod=350,
    alphaVmaxCH4ox=0.07, EaVmaxCH4ox=30,   KMCH4ox=0.01, KmO2_CH4ox=43, kl_CH4prod=3,
    
    reltol=sqrt(.Machine$double.eps),verbose=99
  )
  sink()
  r
 }

