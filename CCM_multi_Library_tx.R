################################################################
# Initialize the C functions for bootstrap CCM
################################################################

program_init_bootstrap<-function() {
  #########################
  # Load CCM function
  #########################
  require(deSolve)
  
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
  
  file3="bmod_single"
  if(is.loaded(file3)) {dyn.unload(paste(file3,".so",sep=""))}
  system(paste("R CMD SHLIB ",file3,".c",sep=""))
  dyn.load(paste(file3,".so",sep=""))
}


program_init_bootstrap_short<-function() {
  #########################
  # Load CCM function
  #########################
  require(deSolve)
  
  file="CCM_bootstrap"
  if(!is.loaded(file)) {
    dyn.load(paste(file,".so",sep=""))
  }
  
  file2="SSR_predict_boot"
  if(!is.loaded(file2)) {
    dyn.load(paste(file2,".so",sep=""))
  }
  
  file3="bmod_single"
  if(!is.loaded(file3)) {
    dyn.load(paste(file3,".so",sep=""))
  }
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
  
  if(tau*(E+1)+predstep>=lengthacceptablelib) {
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
# Generate time series function, C version
################################################################
make_comp_data_boot<-function(a = 1,
                                   w = 1,
                                   K = 10,
                                   m = 0.1,
                                   Q = 2,
                                   r = 1,
                                   S = 100,
                                   x_mean_sd = c(0, 1),
                                   y_mean_sd = c(0, 1),
                                   B=c(rep(1, 12)),
                                   R=100,
                                   times=seq(1, 1000, by = 1),
                                   OUrates=c(0.1, 0.1), #strength of return to mean
                                   Wrates=c(1, 1), #strength of deviation from mean
                                   seednum=1989
) {
  nspec=12 # must have 12 species
  
  #########################
  # Load diffeq function
  #########################
  require(deSolve) # load differentail equation solver package
  if(seednum) {
    set.seed(seednum)
  }
  
  xstr<-NA
  #########################
  # run diffeq function
  #########################
  xopt<-seq(0, 1, length=nspec) #increase this?
  yopt<-seq(y_mean_sd[1]-y_mean_sd[2], y_mean_sd[1]+y_mean_sd[2], length=nspec)[c(6:1, 12:7)]
  
  pars=c(OUrates[1], OUrates[2], Wrates[1],Wrates[2],x_mean_sd[1],  y_mean_sd[1],xstr,  Q,  r,  S,  a,  m,  K,  w,  xopt, yopt)
  
  xlist<- numeric(max(times))
  ylist<- numeric(max(times))
  xmid<-x_mean_sd[1]
  ymid<-y_mean_sd[1]
  
  xlist[1]<-runif(1, min=x_mean_sd[1]-2*x_mean_sd[2], max=x_mean_sd[1]+2*x_mean_sd[2])
  ylist[1]<-runif(1, min=y_mean_sd[1]-2*y_mean_sd[2], max=y_mean_sd[1]+2*y_mean_sd[2])
  Wxlist<-rnorm(max(times), 0, Wrates[1]) #White noise for x
  Wylist<-rnorm(max(times), 0, Wrates[2]) #White noise for y
  for(i in 2:max(times)) {
    xlist[i]<-xlist[i-1]+(xmid-xlist[i-1])*OUrates[1]+Wxlist[i-1]
    ylist[i]<-ylist[i-1]+(ymid-ylist[i-1])*OUrates[2]+Wylist[i-1]
  } #Check out speed of response (Smapping - which work, which don't? Fully stochastic vs. some signal)
  
  xlist<-xlist+sin(1:max(times)/100*pi*2)
  ylist<-ylist+cos(1:max(times)/100*pi*2) 
  
  xlist<-scale(xlist)
  ylist<-scale(ylist)
  
  forcings <- list(x=cbind(times,xlist[(1:max(times))%in%times]))
  
  return(list(pars=pars, forcings=forcings, nspec=nspec, xopt=xopt, yopt=yopt, xlist=xlist, ylist=ylist))
}

make_comp_sim<-function(compsimout, number_of_chains=10,
                        a = 1,
                        w = 1,
                        K = 10,
                        m = 0.1,
                        Q = 2,
                        r = 1,
                        S = 100,
                        x_mean_sd = c(0, 1),
                        y_mean_sd = c(0, 1),
                        xstr=0.5,
                        B=c(rep(1, 12)),
                        R=100,
                        times=seq(1, 1000, by = 1),
                        OUrates=c(0.1, 0.1), #strength of return to mean
                        Wrates=c(1, 1), #strength of deviation from mean
                        seednum=1989) {
  require(deSolve)
  pars<-compsimout$pars
  pars[7]<-xstr
  forcings<-compsimout$forcings
  nspec<-compsimout$nspec
  R<-100
  xopt<-compsimout$xopt
  yopt<-compsimout$yopt
  
  xlist<-compsimout$xlist
  ylist<-compsimout$ylist
  
  out<-NULL
  for(plotiter in 1:number_of_chains) {

    Biter<-B
    parms<-pars
    parms[which(Biter==0)+14]<-0
    
    yini<-c(R, Biter)
    out_tmp   <- ode(y=yini, times=times, func="derivs",
                     parms=parms, dllname="bmod_single", initforc="forcc",
                     forcings=forcings, initfun="parmsc", nout=14,
                     outnames=NULL, jacfun="jac", 
                     method="ode23")
    
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
# CCM environment function
################################################################
#predstep<-1; tau<-1; make_plot<-TRUE; twoway<-FALSE; iterations<-100; maxE<-3; doscale=TRUE; target_sp<-"all"
doCCM_environment<-function(ode_result, target_sp, predstep=10, tau=1, maxE=20, iterations=100, twoway=TRUE, make_plot=TRUE, doscale=TRUE) {
  predstep<-predstep
  out<-ode_result[["out"]]
  pars<-ode_result[["pars"]]
  #########################
  # Run CCM
  #########################
  
  if(target_sp=="all") {
    Biomass_total<-rowSums(out[,3:14])
  } else {
    Biomass_total<-out[ , target_sp+2]
  }
  temperature<-pars$xlist
  moisture<-pars$ylist
  
  if(doscale==TRUE) {
    Biomass_total<-scale(Biomass_total)
    temperature<-scale(temperature)
    moisture<-scale(moisture)
  }
  
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
  reptimes<-length(Biomass_total)/(length(pars$times)+1)
  temperature<-rep(c(temperature, NA), reptimes)
  moisture<-rep(c(moisture, NA), reptimes)
  
  E_temp<-min(which(Emat[,2]>=(quantile(Emat[,2], 0.95, na.rm=T)-0.01)))+1
  E_moist<-min(which(Emat[,3]>=(quantile(Emat[,3], 0.95, na.rm=T)-0.01)))+1
  E_biom<-min(which(Emat[,1]>=(quantile(Emat[,1], 0.95, na.rm=T)-0.01)))+1
  
  
  maxE_use<-max(E_temp, E_moist, E_biom)
  #DesiredL<-sort(unique(c(floor(seq((maxE_use)*tau, length(temperature)/4, length=10)),
  #                        floor(seq(length(temperature)/4, length(temperature)-tau, length=11)))))
  
  T_cause_B<-CCM_boot(A=temperature, B=Biomass_total, E=E_temp, iterations=iterations)#, DesiredL=DesiredL)
  M_cause_B<-CCM_boot(A=moisture, B=Biomass_total, E=E_moist, iterations=iterations)#, DesiredL=DesiredL)
  
  
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
# Test CCM_out
################################################################
###Boot diff version
testccm_boot<-function(ccm_out, iter=1000) {
  lng_T<-length(ccm_out$T_cause_B$rho)
  lng_M<-length(ccm_out$M_cause_B$rho)
  
  
  T_dat<-matrix(nrow=iter, ncol=1)
  M_dat<-matrix(nrow=iter, ncol=1)
  Ttime<-ccm_out$T_cause_B$Lobs
  Mtime<-ccm_out$M_cause_B$Lobs
  Trho<-ccm_out$T_cause_B$rho
  Tsrho<-ccm_out$T_cause_B$sdevrho
  Mrho<-ccm_out$M_cause_B$rho
  Msrho<-ccm_out$M_cause_B$sdevrho
  
  for(i in 1:iter) {
    tmT<-sample((1:length(Ttime))[!is.na(Tsrho)], iter, rep=T)
    tmM<-sample((1:length(Mtime))[!is.na(Msrho)], iter, rep=T)
    
    subsT<-rnorm(n=iter, mean=Trho[tmT], sd=Tsrho[tmT])
    subsM<-rnorm(n=iter, mean=Mrho[tmM], sd=Msrho[tmM])
    
    tmT2<-sample((1:length(ccm_out$T_cause_B$Lobs))[!is.na(ccm_out$T_cause_B$sdevrho)], iter, rep=T)
    tmM2<-sample((1:length(ccm_out$M_cause_B$Lobs))[!is.na(ccm_out$M_cause_B$sdevrho)], iter, rep=T)
    
    subsT<-subsT-rnorm(n=iter, mean=Trho[tmT2], sd=ccm_out$T_cause_B$sdevrho[tmT2])
    subsM<-subsM-rnorm(n=iter, mean=Mrho[tmM2], sd=ccm_out$M_cause_B$sdevrho[tmM2])
    
    tmT<-tmT-tmT2
    tmM<-tmM-tmM2
    
    T_dat[i,1]<-summary(lm(subsT~tmT))$coefficients[2,1]
    M_dat[i,1]<-summary(lm(subsM~tmM))$coefficients[2,1] 
  }
  
  #hist(T_dat[,1])
  #hist(M_dat[,1])
  
  T_curve<-quantile(T_dat, c(0.025, 0.1586553, 0.5, 0.8413447, 0.975), na.rm=TRUE)
  M_curve<-quantile(M_dat, c(0.025, 0.1586553, 0.5, 0.8413447, 0.975), na.rm=TRUE)
  
  maxsd_T<-max(which(is.finite(ccm_out$T_cause_B$sdevrho)))
  maxsd_M<-max(which(is.finite(ccm_out$M_cause_B$sdevrho)))
  
  #Test whether final rho is significantly different from zero
  T_rho<-c(ccm_out$T_cause_B$rho[maxsd_T], ccm_out$T_cause_B$sdevrho[maxsd_T], pnorm(0, ccm_out$T_cause_B$rho[maxsd_T], ccm_out$T_cause_B$sdevrho[maxsd_T]))
  M_rho<-c(ccm_out$M_cause_B$rho[maxsd_M], ccm_out$M_cause_B$sdevrho[maxsd_M], pnorm(0, ccm_out$M_cause_B$rho[maxsd_M], ccm_out$M_cause_B$sdevrho[maxsd_M]))
  
  T_cause<-(T_curve[1]>0)&&(T_rho[1]>0)&&(T_rho[3]<0.05)
  M_cause<-(M_curve[1]>0)&&(M_rho[1]>0)&&(M_rho[3]<0.05)
  
  T_test<-c(T_curve, T_rho, T_cause)
  names(T_test)<-c("sl_025", "sl_160", "sl_500", "sl_840", "sl_975", "rho", "vrho", "prho", "cause")
  M_test<-c(M_curve, M_rho, M_cause)
  names(M_test)<-c("sl_025", "sl_160", "sl_500", "sl_840", "sl_975", "rho", "vrho", "prho", "cause")
  
  return(list(T_test=T_test, M_test=M_test))
}