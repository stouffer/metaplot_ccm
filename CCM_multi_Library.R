#setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")

#gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG      -fpic  -pipe  -g  -c SSR_predict_boot.c -o SSR_predict_boot.o
#gcc -std=gnu99 -shared -o SSR_predict_boot.so SSR_predict_boot.o -L/usr/lib/R/lib -lR

#To Do
# 1) randomize sp community

#program_init()
#ode_result<-make_comp_data(seednum=2000)
#plot_output(ode_result)
#ccm_out<-doCCM_environment(ode_result, "all")
#ssr_out<-ssr_data(ccm_out, predstepmax=10)
#plot_ccm(ccm_out)

################################################################
# Initialize the C functions for bootstrap CCM
################################################################

program_init_bootstrap<-function() {
  #########################
  # Load CCM function
  #########################
  file="CCM_bootstrap"
  if(is.loaded(file)) {dyn.unload(paste(file,".so",sep=""))}
  if(sum(grep(".so", dir()))>0) {
    system("rm *.so"); system("rm *.o")
  }
  system(paste("R CMD SHLIB ",file,".c",sep=""))
  dyn.load(paste(file,".so",sep=""))
  
  file2="SSR_predict_boot"
  if(is.loaded(file2)) {dyn.unload(paste(file2,".so",sep=""))}
  system(paste("R CMD SHLIB ",file2,".c",sep=""))
  dyn.load(paste(file2,".so",sep=""))
  
  file3="bmod"
  if(is.loaded(file3)) {dyn.unload(paste(file3,".so",sep=""))}
  system(paste("R CMD SHLIB ",file3,".c",sep=""))
  dyn.load(paste(file3,".so",sep=""))
}


program_init_bootstrap_short<-function() {
  #########################
  # Load CCM function
  #########################
  file="CCM_bootstrap"
  if(is.loaded(file)) {dyn.unload(paste(file,".so",sep=""))}
  dyn.load(paste(file,".so",sep=""))
  
  file2="SSR_predict_boot"
  if(is.loaded(file2)) {dyn.unload(paste(file2,".so",sep=""))}
  dyn.load(paste(file2,".so",sep=""))
  
  file3="bmod"
  if(is.loaded(file3)) {dyn.unload(paste(file3,".so",sep=""))}
  dyn.load(paste(file3,".so",sep=""))
}


################################################################
# Run bootstrap CCM
################################################################

CCM_boot<-function(A, B, E, tau=1, DesiredL=((tau*(E-1)+(E+1)):length(A)-E+2), iterations=100) {
  #Note - input must have "NA" between any gaps in the library (e.g., if you are "stacking" plots)
  ## Initialize parameters
  plengtht=length(A[!is.na(A)]); pLibLength=length(A[!is.na(A)]); Aest=rep(0, length(A)); rho=Aest; varrho=Aest
  
  # Check to make sure library is not too long
  if(plengtht>pLibLength) {plengtht=pLibLength}
  
  # Make list of acceptable starting positions
  gapdist<-tau*(E-1)
  acceptablelib<-as.numeric(is.finite(A))
  lA<-length(A)
  for(i in 1:gapdist) {
    acceptablelib<-acceptablelib*as.numeric(is.finite(c(rep(NA, i),A[-c((lA-i+1):lA)])))
  }
  acceptablelib<-which(acceptablelib>0)-1 #Convert into positions in C array
  acceptablelib<-acceptablelib[acceptablelib<((plengtht-1)-(tau))]
  lengthacceptablelib<-length(acceptablelib)
  
  DesiredL<-DesiredL+E-2 #Convert to positions in C array
  
  for(i in 1:length(DesiredL)) { #Load nearby points from acceptablelib vector
    DesiredL[i]<-acceptablelib[which.min(abs(acceptablelib-DesiredL[i]))] 
  }
  # Update input to match actual indices
  A[is.na(A)]<-0; B[is.na(B)]<-0 #Make vectors "C safe"
  
  if(tau*(E+1)>lengthacceptablelib) {
    print(paste("Error - too few records to test E =", E, "and tau =", tau))
    return(out=list(A=A, Aest=NA, B=B, rho=NA, varrho=NA, sdevrho=NA, Lobs=NA, E=out$E, tau=tau, FULLinfo=NA, rejectedL=NA))
  } else {
    
    out<-.C("CCM_bootstrap", A=as.double(A), Aest=as.double(Aest), B=as.double(B), rho=as.double(rho), E=as.integer(E), tau=as.integer(tau),
            plength=as.integer(plengtht), pLibLength=as.integer(pLibLength),DesiredL=as.integer(DesiredL), plengthDesL=as.integer(length(DesiredL)),
            piterations=as.integer(iterations), varrho=as.double(varrho), acceptablelib=as.integer(acceptablelib), plengthacceptablelib=as.integer(lengthacceptablelib))
    out$Aest[out$Aest==0]<-NA #Mark values that were not calculated  
    sdevrho<-sqrt(out$varrho[out$rho!=0]) #Calculate standard deviation
    return(list(A=out$A, Aest=out$Aest, B=out$B, rho=out$rho[out$rho!=0], varrho=out$varrho[out$rho!=0], sdevrho=sdevrho, Lobs=(1:length(A))[out$rho!=0]-E+1, E=out$E, tau=out$tau, FULLinfo=out))
  }
}

SSR_pred_boot<-function(A, B=A, E, tau=1, predstep=1, matchSugi=1) {
  repvec=as.numeric((sum(A[is.finite(A)]==B[is.finite(B)])==length(A[is.finite(A)]))&(length(A[is.finite(A)])==length(B[is.finite(B)])))
  
  #Predict elements of A using B
  #If A=B, uses cross-validation
  #matchSugi=1 removes only point i in cross validation - if 0, then removes all points within X(t-(E-1)):X(t+1)
  #repvec=0 if A and B are not the same vector (i.e, not using "leave one out cross validation")
  # Make list of acceptable starting positions
  gapdist<-tau*(E-1)+predstep
  acceptablelib<-as.numeric(is.finite(A))
  lA<-length(A)
  for(i in 1:gapdist) {
    acceptablelib<-acceptablelib*as.numeric(is.finite(c(rep(NA, i),A[-c((lA-i+1):lA)])))
  }
  acceptablelib<-which(acceptablelib>0)-1 #Convert into positions in C array
  lengthacceptablelib<-length(acceptablelib)
  
  if(tau*(E+1+predstep)>lengthacceptablelib) {
    print(paste("Error - too few records to test E =", E, "tau =", tau, "and prestep =", predstep))
    return(out=list(A=A, Aest=NA, B=B, E=E, tau=tau, pBlength=length(B), pAlength=length(A), predstep=predstep,
                    prepvec=repvec, pmatchSugi=matchSugi, acceptablelib=acceptablelib, plengthacceptablelib=lengthacceptablelib, rho=NA))
  } else { # Don't attempt to run algorithm using more lags than datapoints
    A[is.na(A)]<-0; B[is.na(B)]<-0 #Make vectors "C safe"  
    
    out<-.C("SSR_predict_boot", A=as.double(A), Aest=as.double(rep(0, length(A))), B=as.double(B), E=as.integer(E),
            tau=as.integer(tau),pBlength=as.integer(length(B)), pAlength=as.integer(length(A)), predstep=as.integer(predstep),
            prepvec=as.integer(repvec), pmatchSugi=as.integer(matchSugi), acceptablelib=as.integer(acceptablelib), plengthacceptablelib=as.integer(lengthacceptablelib))
    out$rho<-cor(out$A, out$Aest)
    out$Aest[out$Aest==0]<-NA
    return(out)
  }
}

################################################################
# Generate time series function, conventional
################################################################
make_comp_data_boot<-function(nspec=12,
                              a = 1,
                              w = 1,
                              K = 10,
                              m = 0.1,
                              Q = 2,
                              r = 1,
                              S = 100,
                              x_mean_sd = c(25, 1),
                              y_mean_sd = c(25, 1),
                              xstr=0,
                              B=c(rep(0.01, nspec)),
                              R=100,
                              times=seq(1, 100, by = 1),
                              OUrates=c(0.1, 0.1), #strength of return to mean
                              Wrates=c(1, 1), #strength of deviation from mean
                              seednum=1989,
                              number_of_chains=10 #Number of plots
) {
  
  #########################
  # Load diffeq function
  #########################
  times<-seq(1, max(times), by=1)
  pars=list(a = a,
            w = w,
            K = K,
            m = m,
            Q = Q,
            r = r,
            S = S,
            x_mean_sd = x_mean_sd,
            y_mean_sd = y_mean_sd,
            times=times)
  
  require(deSolve) # load differentail equation solver package
  if(seednum) {
    set.seed(seednum)
  }
  
  Bmod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      tpos<-which(times==floor(Time))
      xt<-xlist[tpos]
      yt<-ylist[tpos]
      
      gix<-r*exp(-0.5*((abs(xopt-xt)*(xstr)+abs(yopt-yt)*(1-xstr))/w)^2)
      
      dB<-(gix*(State[1]/(State[1]+K))-m)*State[-1]
      dR<-a*(S-State[1])-sum(Q*gix*(State[1]/(State[1]+K)*State[-1]))
      
      return(list(c(dR, dB)))
    })
  }
  
  #########################
  # make climate variables
  #########################
  xlist<- numeric(max(times))
  ylist<- numeric(max(times))
  
  xmid<-x_mean_sd[1]
  ymid<-y_mean_sd[1]
  
  xlist[1]<-runif(1, min=max(x_mean_sd[1]-2*x_mean_sd[2],0), max=x_mean_sd[1]+2*x_mean_sd[2])
  ylist[1]<-runif(1, min=max(y_mean_sd[1]-2*y_mean_sd[2],0), max=y_mean_sd[2]+2*y_mean_sd[2])
  Wxlist<-rnorm(max(times), 0, Wrates[1]) #White noise for x
  Wylist<-rnorm(max(times), 0, Wrates[2]) #White noise for y
  for(i in 2:max(times)) {
    xlist[i]<-xlist[i-1]+(xmid-xlist[i-1])*OUrates[1]+Wxlist[i-1]
    ylist[i]<-ylist[i-1]+(ymid-ylist[i-1])*OUrates[2]+Wylist[i-1]
    
  } #Check out speed of response (Smapping - which work, which don't? Fully stochastic vs. some signal)
  
  #########################
  # run diffeq function
  #########################
  pars$times<-times
  pars$xstr<-xstr
  
  pars$xlist<-xlist
  pars$ylist<-ylist
  
  pars$xopt<-rnorm(nspec, pars$x_mean_sd[1], pars$x_mean_sd[2])
  pars$yopt<-rnorm(nspec, pars$y_mean_sd[1], pars$y_mean_sd[2])
  
  out<-NULL
  for(plotiter in 1:number_of_chains) {
    Biter<-B*sample(c(1,0), length(B), rep=T)
    yini<-c(R, Biter)
    
    
    out_tmp   <- ode(yini, times, Bmod, pars)
    out<-rbind(out, as.matrix(out_tmp), NA)
  }
  
  return(list(out=out, pars=pars))
}

################################################################
# Generate time series function, C version
################################################################
make_comp_data_boot_Cfxn<-function(a = 1,
                                   w = 1,
                                   K = 10,
                                   m = 0.1,
                                   Q = 2,
                                   r = 1,
                                   S = 100,
                                   x_mean_sd = c(25, 1),
                                   y_mean_sd = c(25, 1),
                                   xstr=0.5,
                                   B=c(rep(1, 12)),
                                   R=100,
                                   times=seq(1, 1000, by = 1),
                                   OUrates=c(0.1, 0.1), #strength of return to mean
                                   Wrates=c(1, 1), #strength of deviation from mean
                                   seednum=1989,
                                   number_of_chains=10 #Number of plots
) {
  nspec=12 # must have 12 species
  
  #########################
  # Load diffeq function
  #########################
  require(deSolve) # load differentail equation solver package
  if(seednum) {
    set.seed(seednum)
  }
  
  #########################
  # run diffeq function
  #########################
  xopt<-rnorm(nspec, x_mean_sd[1], x_mean_sd[2])
  yopt<-rnorm(nspec, y_mean_sd[1], y_mean_sd[2])
  
  pars=c(OUrates[1], OUrates[2], Wrates[1],Wrates[2],x_mean_sd[1],  y_mean_sd[1],xstr,  Q,  r,  S,  a,  m,  K,  w,  xopt, yopt)
  
  xlist<- numeric(max(times))
  ylist<- numeric(max(times))
  xmid<-x_mean_sd[1]
  ymid<-y_mean_sd[1]
  
  xlist[1]<-runif(1, min=max(x_mean_sd[1]-2*x_mean_sd[2],0), max=x_mean_sd[1]+2*x_mean_sd[2])
  ylist[1]<-runif(1, min=max(y_mean_sd[1]-2*y_mean_sd[2],0), max=y_mean_sd[1]+2*y_mean_sd[2])
  Wxlist<-rnorm(max(times), 0, Wrates[1]) #White noise for x
  Wylist<-rnorm(max(times), 0, Wrates[2]) #White noise for y
  for(i in 2:max(times)) {
    xlist[i]<-xlist[i-1]+(xmid-xlist[i-1])*OUrates[1]+Wxlist[i-1]
    ylist[i]<-ylist[i-1]+(ymid-ylist[i-1])*OUrates[2]+Wylist[i-1]
  } #Check out speed of response (Smapping - which work, which don't? Fully stochastic vs. some signal)
  
  forcings <- list(x=cbind(times,xlist), y=cbind(times,ylist))
  #matplot(times, cbind(xlist, ylist), type="s", ylab="Climvar")
  
  out<-NULL
  for(plotiter in 1:number_of_chains) {
    Biter<-B*sample(c(1,0), length(B), rep=T)
    parms<-pars
    parms[which(Biter==0)+14]<-0
    
    yini<-c(R, Biter)
    out_tmp   <- ode(y=yini, times=times, func="derivs",
                     parms=parms, dllname="bmod", initforc="forcc",
                     forcings=forcings, initfun="parmsc", nout=14)#,
                     #outnames=NULL, jacfun="jac"), 
                     #method="ode23")
    
    #plot(out_tmp[,1], out_tmp[,2], type="l", ylab="R", xlab="time")
    #matplot(out_tmp[,1], out_tmp[,3:(3+12)], type="l", ylab="Sp", xlab="time")
    #matplot(out_tmp[,1], out_tmp[,which(Biter==0)+2], type="l", ylab="dead sp", xlab="time")
    out<-rbind(out, as.matrix(out_tmp), NA)
  }
  
  pars_save=list(a = a,
                 w = w,
                 K = K,
                 m = m,
                 Q = Q,
                 r = r,
                 S = S,
                 x_mean_sd = x_mean_sd,
                 y_mean_sd = y_mean_sd,
                 times=times,
                 xstr=xstr,
                 xlist=xlist,
                 ylist=ylist,
                 xopt=xopt,
                 yopt=yopt)
  
  return(list(out=out, pars=pars_save))
}

################################################################
# Make plots of results
################################################################
plot_output_boot<-function(ode_result) {
  ########################################
  ## Plotting
  ########################################
  out<-ode_result$out
  pars<-ode_result$pars
  
  matplot(out[ , 1], out[ , -c(1,2)], type = "l", xlab = "time", ylab = "B",
          main = "Bmod", lwd = 2)
  plot(out[,1], out[,2], type="l", xlab="time", ylab="R")
  plot(out[ , 1], rowSums(out[ , -c(1,2)]), type = "l", xlab = "time", ylab = "B",
       main = "Bmod", lwd = 2)
  
  timelist<-1:(min(which(is.na(out[,1])))-1)
  plot(out[timelist,1], pars$xlist, type="s", xlab="time", ylab="x")
  plot(out[timelist,1], pars$ylist, type="s", xlab="time", ylab="y")
  
  #Plot niche space
  plot(range(pars$xopt), c(0,1), xlab="x", ylab="f", type="n", main="x niche space")
  dlist<-seq(20,30,by=0.1)
  for(i in pars$xopt) {
    lines(dlist, pars$r*exp(-0.5*((abs(i-dlist))/pars$w)^2))
  }
  
  plot(range(pars$yopt), c(0,1), xlab="y", ylab="f", type="n", main="y niche space")
  dlist<-seq(20,30,by=0.1)
  for(i in pars$yopt) {
    lines(dlist, pars$r*exp(-0.5*((abs(i-dlist))/pars$w)^2))
  }
}

################################################################
# Make E plot
################################################################
#target_sp<-"all"; predstep=10; tau=1; maxE=20
makeEplot_environment_boot<-function(ode_result, target_sp, predstep=10, tau=1, maxE=20, matchSugi=1) {
  out<-ode_result[["out"]]
  pars<-ode_result[["pars"]]
  
  if(target_sp=="all") {
    Biomass_total<-rowSums(out[ , -c(1,2)])
  } else {
    Biomass_total<-out[ , target_sp+2]
  }
  temperature<-pars$xlist
  moisture<-pars$ylist
  
  Emat<-matrix(nrow=maxE-1, ncol=3); colnames(Emat)<-c("B", "T", "M")
  
  for(E in 2:maxE) {
    Emat[E-1,"B"]<-SSR_pred_boot(A=Biomass_total, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"T"]<-SSR_pred_boot(A=temperature, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"M"]<-SSR_pred_boot(A=moisture, E=E, predstep=predstep, tau=tau)$rho
  }
  matplot(2:(maxE), Emat, type="l", col=1:3, lty=1:3,
          xlab="E", ylab="rho", lwd=2,
          main=paste("pred. step = ", predstep, "; tau = ", tau, sep=""))
  legend("bottomleft", c("Biomass", "Temperature", "Moisture"), lty=1:3, col=1:3, lwd=2, bty="n")
  
  E_temp<-min(which(Emat[,2]>=(quantile(Emat[,2], 0.95, na.rm=T)-0.01)))+1
  E_moist<-min(which(Emat[,3]>=(quantile(Emat[,3], 0.95, na.rm=T)-0.01)))+1
  E_biom<-min(which(Emat[,1]>=(quantile(Emat[,1], 0.95, na.rm=T)-0.01)))+1
  
  return(list(Emat=Emat, Euse=c(E_temp=E_temp, E_moist=E_moist, E_biom=E_biom)))
}






################################################################
# Make E plot species
################################################################
makeEplot_species<-function(ode_result, target_sp_1, target_sp_2, predstep=10, tau=1, maxE=20) {
  out<-ode_result[["out"]]
  pars<-ode_result[["pars"]]
  #########################
  # Run CCM
  #########################
  
  Biomass_total_1<-out[ , target_sp_1+2]
  Biomass_total_2<-out[ , target_sp_2+2]
  
  DesiredL<-unique(round(c(seq(10, median(pars$times), length=20), seq(median(pars$times), max(pars$times)-10, length=10))),0)
  
  Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("B", "T")
  
  for(E in 2:maxE) {
    Emat[E-1,"B"]<-SSR_pred_boot(A=Biomass_total_1, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"T"]<-SSR_pred_boot(A=Biomass_total_2, E=E, predstep=predstep, tau=tau)$rho
  }
  matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,
          xlab="E", ylab="rho", lwd=2,
          main=paste("pred. step = ", predstep, "; tau = ", tau, sep=""))
  legend("bottomleft", c("Biomass 1", "Biomass 2"), lty=1:2, col=1:2, lwd=2, bty="n")
  
  E_biom_2<-min(which(Emat[,2]>=(quantile(Emat[,2], 0.95, na.rm=T)-0.01)))+1
  E_biom_1<-min(which(Emat[,1]>=(quantile(Emat[,1], 0.95, na.rm=T)-0.01)))+1
  
  return(list(Emat=Emat, Euse=c(E_biom_1=E_biom_1, E_biom_2=E_biom_2)))
}


################################################################
# CCM environment function
################################################################
doCCM_environment<-function(ode_result, target_sp, predstep=10, tau=1, maxE=20, iterations=100, twoway=TRUE, make_plot=TRUE) {
  predstep<-predstep
  out<-ode_result[["out"]]
  pars<-ode_result[["pars"]]
  #########################
  # Run CCM
  #########################
  
  if(target_sp=="all") {
    Biomass_total<-rowSums(out[ , -c(1,2)])
  } else {
    Biomass_total<-out[ , target_sp+2]
  }
  temperature<-pars$xlist
  moisture<-pars$ylist
  
  Emat<-matrix(nrow=maxE-1, ncol=3); colnames(Emat)<-c("B", "T", "M")
  
  for(E in 2:maxE) {
    Emat[E-1,"B"]<-SSR_pred_boot(A=Biomass_total, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"T"]<-SSR_pred_boot(A=temperature, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"M"]<-SSR_pred_boot(A=moisture, E=E, predstep=predstep, tau=tau)$rho
  }
  if(make_plot) {
  matplot(2:maxE, Emat, type="l", col=1:3, lty=1:3,
          xlab="E", ylab="rho", lwd=2,
          main=paste("pred. step = ", predstep, "; tau = ", tau, sep=""))
  legend("bottomleft", c("Biomass", "Temperature", "Moisture"), lty=1:3, col=1:3, lwd=2, bty="n")
  }
  
  ## Add on T and M libraries
  reptimes<-length(Biomass_total)/(length(temperature)+1)
  temperature<-rep(c(temperature, NA), reptimes)
  moisture<-rep(c(moisture, NA), reptimes)
  
  E_temp<-min(which(Emat[,2]>=(quantile(Emat[,2], 0.95, na.rm=T)-0.01)))+1
  E_moist<-min(which(Emat[,3]>=(quantile(Emat[,3], 0.95, na.rm=T)-0.01)))+1
  E_biom<-min(which(Emat[,1]>=(quantile(Emat[,1], 0.95, na.rm=T)-0.01)))+1
  
  
  maxE_use<-max(E_temp, E_moist, E_biom)
  DesiredL<-sort(unique(c(floor(seq((maxE_use)*tau, length(temperature)/4, length=10)),
                          floor(seq(length(temperature)/4, length(temperature)-tau, length=11)))))
  
  T_cause_B<-CCM_boot(A=temperature, B=Biomass_total, E=E_temp, DesiredL=DesiredL, iterations=iterations)
  M_cause_B<-CCM_boot(A=moisture, B=Biomass_total, E=E_moist, DesiredL=DesiredL, iterations=iterations)
  
  
  if(twoway) {
    B_cause_T<-CCM_boot(A=Biomass_total, B=temperature, E=E_biom, DesiredL=DesiredL, iterations=iterations)
    B_cause_M<-CCM_boot(A=Biomass_total, B=moisture, E=E_biom, DesiredL=DesiredL, iterations=iterations)
    
    return(list(T_cause_B=T_cause_B, M_cause_B=M_cause_B,
                B_cause_T=B_cause_T, B_cause_M=B_cause_M, Emat=Emat, Euse=c(E_temp=E_temp, E_moist=E_moist, E_biom=E_biom),
                datalist=list(Biomass_total=Biomass_total, temperature=temperature, moisture=moisture)))
  } else {
    return(list(T_cause_B=T_cause_B, M_cause_B=M_cause_B,
                B_cause_T=NULL, B_cause_M=NULL, Emat=Emat, Euse=c(E_temp=E_temp, E_moist=E_moist, E_biom=E_biom),
                datalist=list(Biomass_total=Biomass_total, temperature=temperature, moisture=moisture)))
  }
}

################################################################
# CCM species function
################################################################
doCCM_species<-function(ode_result, target_sp_1, target_sp_2, predstep=10, tau=1, maxE=20) {
  out<-ode_result[["out"]]
  pars<-ode_result[["pars"]]
  #########################
  # Run CCM
  #########################
  
  Biomass_total_1<-out[ , target_sp_1+2]
  Biomass_total_2<-out[ , target_sp_2+2]
  
  DesiredL<-unique(round(c(seq(10, median(pars$times), length=20), seq(median(pars$times), max(pars$times)-10, length=10))),0)
  if(max(DesiredL)>length(Biomass_total)) {DesiredL<-1:length(Biomass_total)}
  
  Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("B", "T")
  
  for(E in 2:maxE) {
    Emat[E-1,"B"]<-SSR_pred(A=Biomass_total_1, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"T"]<-SSR_pred(A=Biomass_total_2, E=E, predstep=predstep, tau=tau)$rho
  }
  matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,
          xlab="E", ylab="rho", lwd=2,
          main=paste("pred. step = ", predstep, "; tau = ", tau, sep=""))
  legend("bottomleft", c("Biomass 1", "Biomass 2"), lty=1:2, col=1:2, lwd=2, bty="n")
  
  E_biom_2<-min(which(Emat[,2]>=(quantile(Emat[,2], 0.95)-0.01)))+1
  E_biom_1<-min(which(Emat[,1]>=(quantile(Emat[,1], 0.95)-0.01)))+1
  
  B1_cause_B2<-CCM_varL(A=Biomass_total_1, B=Biomass_total_2, E=E_biom_1, DesiredL=DesiredL)
  B2_cause_B1<-CCM_varL(A=Biomass_total_2, B=Biomass_total_1, E=E_biom_2, DesiredL=DesiredL)
  
  return(list(B1_cause_B2, B2_cause_B1, Emat=Emat, Euse=c(E_biom_1=E_biom_1, E_biom_2=E_biom_2),
              datalist=list(Biomass_total_1=Biomass_total_1, Biomass_total_2=Biomass_total_2)))
}

################################################################
# Check for existence of manifold
################################################################
ssr_data<-function(ccm_out, predstepmax=10, tau=1) {
  
  pred_mat<-matrix(nrow=predstepmax, ncol=3)
  colnames(pred_mat)<-c("B", "T", "M")
  for(predstep in 1:predstepmax) {
    pred_mat[predstep,"B"]<-SSR_pred_boot(A=ccm_out$datalist$Biomass_total, E=ccm_out$Euse["E_biom"], predstep=predstep, tau=tau)$rho
    pred_mat[predstep,"T"]<-SSR_pred_boot(A=ccm_out$datalist$temperature, E=ccm_out$Euse["E_temp"], predstep=predstep, tau=tau)$rho
    pred_mat[predstep,"M"]<-SSR_pred_boot(A=ccm_out$datalist$moisture, E=ccm_out$Euse["E_moist"], predstep=predstep, tau=tau)$rho
  }
  
  matplot(1:predstepmax, pred_mat, type="l", col=1:3, lty=1:3,
          xlab="predstep", ylab="rho", lwd=2,
          main=paste("Data SSR"))
  legend("bottomleft", c("Biomass", "Temperature", "Moisture"), lty=1:3, col=1:3, lwd=2, bty="n")
  
  return(pred_mat)
}



################################################################
# Plot ccm
################################################################
plot_ccm<-function(ccm_out, ylimits=c(0, 1), twoway=TRUE) {
  #########################
  # Plot CCM
  #########################
  if(twoway) {
    
    with(ccm_out, {
      minL<-max(length(T_cause_B$Lobs), length(T_cause_B$rho), length(B_cause_T$rho), length(M_cause_B$rho), length(B_cause_M$rho))
      matplot(T_cause_B$Lobs[1:minL], cbind(T_cause_B$rho[1:minL],T_cause_B$rho[1:minL]+T_cause_B$sdevrho[1:minL],T_cause_B$rho[1:minL]-T_cause_B$sdevrho[1:minL],
                                            B_cause_T$rho[1:minL],B_cause_T$rho[1:minL]+B_cause_T$sdevrho[1:minL],B_cause_T$rho[1:minL]-B_cause_T$sdevrho[1:minL],
                                            M_cause_B$rho[1:minL],M_cause_B$rho[1:minL]+M_cause_B$sdevrho[1:minL],M_cause_B$rho[1:minL]-M_cause_B$sdevrho[1:minL],
                                            B_cause_M$rho[1:minL],B_cause_M$rho[1:minL]+B_cause_M$sdevrho[1:minL],B_cause_M$rho[1:minL]-B_cause_M$sdevrho[1:minL]),
              type="l", main="CCM results", xlab="Library length", ylab="rho", lwd=2,
              ylim=ylimits, lty=c(1,3,3,1,3,3,1,3,3,1,3,3), col=c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3)))
      legend("topleft", c("Temp. causes Biom.", "Biom. causes Temp.", "Moist. causes Biom.", "Biom. causes Moist."), col=c(1,2,3,4), lty=c(1), lwd=2, bty="n")
      abline(h=0, lty=3, col="darkgrey", lwd=2)
    })
    
  } else {
    with(ccm_out, {
      minL<-max(length(T_cause_B$Lobs), length(T_cause_B$rho), length(B_cause_T$rho), length(M_cause_B$rho), length(B_cause_M$rho))
      matplot(T_cause_B$Lobs[1:minL], cbind(T_cause_B$rho[1:minL],T_cause_B$rho[1:minL]+T_cause_B$sdevrho[1:minL],T_cause_B$rho[1:minL]-T_cause_B$sdevrho[1:minL],
                                            M_cause_B$rho[1:minL],M_cause_B$rho[1:minL]+M_cause_B$sdevrho[1:minL],M_cause_B$rho[1:minL]-M_cause_B$sdevrho[1:minL]),
              type="l", main="", xlab="Library length", ylab="rho", lwd=2,
              ylim=ylimits, lty=c(1,3,3,1,3,3,1,3,3), col=c(rep(1, 3), rep(2, 3)))
      abline(h=0, lty=3, col="darkgrey", lwd=2)
      par(xpd=TRUE)
      legend(mean(T_cause_B$Lobs[1:minL], na.rm=T), max(ylimits)+diff(ylimits)/3,
             c("Temp. causes Biom.", "Moist. causes Biom."), col=c(1,2), lty=c(1), lwd=2, bty="n",
             xjust=0.5, yjust=0.5)
      par(xpd=FALSE)
    })
  }
}



################################################################
# Check for existence of manifold SPECIES
################################################################
ssr_sp<-function(ccm_sp, predstepmax=10, tau=1) {
  
  pred_mat<-matrix(nrow=predstepmax, ncol=2)
  colnames(pred_mat)<-c("B", "T")
  for(predstep in 1:predstepmax) {
    pred_mat[predstep,"B"]<-SSR_pred_boot(A=ccm_sp$datalist$Biomass_total_1, E=ccm_sp$Euse["E_biom_1"], predstep=predstep, tau=tau)$rho
    pred_mat[predstep,"T"]<-SSR_pred_boot(A=ccm_sp$datalist$Biomass_total_2, E=ccm_sp$Euse["E_biom_2"], predstep=predstep, tau=tau)$rho
  }
  
  matplot(1:predstepmax, pred_mat, type="l", col=1:2, lty=1:2,
          xlab="predstep", ylab="rho", lwd=2,
          main=paste("Data SSR"))
  legend("bottomleft", c("Biomass 1", "Biomass 2"), lty=1:2, col=1:2, lwd=2, bty="n")
  
  return(pred_mat)
}


################################################################
# Plot ccm SPECIES
################################################################
plot_sp<-function(ccm_sp, ylimits=c(0, 1)) {
  #########################
  # Plot CCM
  #########################
  with(ccm_sp, {
    minL<-max(length(B1_cause_B2$Lobs), length(B1_cause_B2$rho), length(B2_cause_B1$rho))
    matplot(B1_cause_B2$Lobs[1:minL], cbind(B1_cause_B2$rho[1:minL], B2_cause_B1$rho[1:minL]), type="l", main="CCM results", xlab="Library length", ylab="rho", lwd=2,
            ylim=ylimits, lty=c(1,2), col=c("black", "red"))
    legend("topleft", c("B1. causes B2.", "B2. causes B1."), col=c("black", "red"), lty=c(1,2), lwd=2, bty="n")
    abline(h=0, lty=3, col="darkgrey", lwd=2)
  })
}


################################################################
# ODE out average
################################################################
odeout_average<-function(ode_result, env_avg=100) {
  ode_result$pars$times<-c(seq(1, max(ode_result$pars$times), by=env_avg), max(ode_result$pars$times))[-1]
  ode_result$out<-(ode_result$out)[ode_result$pars$times,]
  
  average_list<-cut(1:max(ode_result$pars$times), breaks=c(1, ode_result$pars$times))
  ode_result$pars$xlist<-aggregate(ode_result$pars$xlist, by=list(average_list), mean)$x
  ode_result$pars$ylist<-aggregate(ode_result$pars$ylist, by=list(average_list), mean)$x
  
  return(ode_result)
}



################################################################
# Test CCM_out
################################################################
testccm<-function(ccm_out) {
  lng_T<-length(ccm_out$T_cause_B$rho)
  
  lng_M<-length(ccm_out$M_cause_B$rho)
  
  #Test for increasing function, by estimating magnitudes of second derivative
  T_curve<-summary(lm(diff(diff(ccm_out$T_cause_B$rho))~ccm_out$T_cause_B$Lobs[c(-1, -lng_T)]))$coefficients[cbind(c(1,1), c(1,4))]
  M_curve<-summary(lm(diff(diff(ccm_out$M_cause_B$rho))~ccm_out$M_cause_B$Lobs[c(-1, -lng_M)]))$coefficients[cbind(c(1,1), c(1,4))]
  
  
  #Test whether final rho is significantly different from zero
  T_rho<-c(ccm_out$T_cause_B$rho[lng_T], pnorm(0, ccm_out$T_cause_B$rho[lng_T], ccm_out$T_cause_B$sdevrho[lng_T]))
  M_rho<-c(ccm_out$M_cause_B$rho[lng_M], pnorm(0, ccm_out$M_cause_B$rho[lng_M], ccm_out$M_cause_B$sdevrho[lng_M]))

  T_cause<-(T_curve[1]<0)&&(T_curve[2]<0.05)&&(T_rho[1]>0)&&(T_rho[2]<0.05)
  M_cause<-(M_curve[1]<0)&&(M_curve[2]<0.05)&&(M_rho[1]>0)&&(M_rho[2]<0.05)
    
  T_test<-c(T_curve, T_rho, T_cause)
    names(T_test)<-c("diff", "pdiff", "rho", "prho", "cause")
  M_test<-c(M_curve, M_rho, M_cause)
  names(M_test)<-c("diff", "pdiff", "rho", "prho", "cause")
  
  return(list(T_test=T_test, M_test=M_test))
}





