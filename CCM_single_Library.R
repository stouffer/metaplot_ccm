#setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")

#program_init()
#ode_result<-make_comp_data(seednum=2000)
#plot_output(ode_result)
#ccm_out<-doCCM_environment(ode_result, "all")
#ssr_out<-ssr_data(ccm_out, predstepmax=10)
#plot_ccm(ccm_out)

################################################################
# Initialize the C functions for CCM
################################################################
program_init<-function() {
  #########################
  # Load CCM function
  #########################
  file="CCM_varL_130409"
  if(is.loaded(file)) {dyn.unload(paste(file,".so",sep=""))}
  if(sum(grep(".so", dir()))>0) {
    system("rm *.so"); system("rm *.o")
  }
  system(paste("R CMD SHLIB ",file,".c",sep=""))
  dyn.load(paste(file,".so",sep=""))
  
  
  file2="SSR_predict_130423"
  if(is.loaded(file2)) {dyn.unload(paste(file2,".so",sep=""))}
  system(paste("R CMD SHLIB ",file2,".c",sep=""))
  dyn.load(paste(file2,".so",sep=""))
}

CCM_varL<-function(A, B, E, tau=1, DesiredL=((tau*(E-1)+(E+1)):length(A)-E+2), plengtht=length(A), pLibLength=length(A), Aest=rep(0, length(A)), rho=rep(0,length(A))) {
  if(plengtht>pLibLength) {plengtht=pLibLength}
  DesiredL<-DesiredL+E-2
  out<-.C("CCM_varL_130409", A=as.double(A), Aest=as.double(Aest), B=as.double(B), rho=as.double(rho), E=as.integer(E), tau=as.integer(tau),
          plength=as.integer(plengtht), pLibLength=as.integer(pLibLength),DesiredL=as.integer(DesiredL), plengthDesL=as.integer(length(DesiredL)))
  out$Aest[out$Aest==0]<-NA #Mark values that were not calculated  
  return(list(A=out$A, Aest=out$Aest, B=out$B, rho=out$rho[out$rho!=0], Lobs=(1:length(A))[out$rho!=0]-E+1, E=out$E, tau=out$tau, FULLinfo=out))
}



SSR_pred<-function(A, B=A, E, tau=1, predstep=1, repvec=as.numeric((sum(A==B)==length(A))&(length(A)==length(B))), matchSugi=1) {
  #Predict elements of A using B
  #If A=B, uses cross-validation
  #matchSugi=1 removes only point i in cross validation - if 0, then removes all points within X(t-(E-1)):X(t+1)
  #repvec=0 if A and B are not the same vector (i.e, not using "leave one out cross validation")
  out<-.C("SSR_predict_130423", A=as.double(A), Aest=as.double(rep(0, length(A))), B=as.double(B), E=as.integer(E),
          tau=as.integer(tau),pBlength=as.integer(length(B)), pAlength=as.integer(length(A)), predstep=as.integer(predstep),
          prepvec=as.integer(repvec), pmatchSugi=as.integer(matchSugi))
  out$rho<-cor(out$A, out$Aest)
  out$Aest[out$Aest==0]<-NA
  return(out)
}



################################################################
# Generate time series function
################################################################
make_comp_data<-function(nspec=12,
                         a = 1,
                         w = 1,
                         K = 10,
                         m = 0.1,
                         Q = 2,
                         r = 1,
                         S = 100,
                         x_mean_sd = c(20, 30),
                         y_mean_sd = c(0, 10),
                         xstr=0,
                         B=c(rep(0.01, nspec)),
                         R=100,
                         times=seq(1, 200, by = 1),
                         OUrates=c(0.1, 0.1), #strength of return to mean
                         Wrates=c(1, 1), #strength of deviation from mean
                         seednum=1989
                         ) {
  
  #########################
  # Load diffeq function
  #########################
  
  pars=list(a = a,
            w = w,
            K = K,
            m = m,
            Q = Q,
            r = r,
            S = S,
            x_mean_sd = x_mean_sd,
            y_mean_sd = y_mean_sd)
  
  require(deSolve) # load differentail equation solver package
  set.seed(seednum)
  
  Bmod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      xt<-xlist[floor(Time)]
      yt<-ylist[floor(Time)]
      
      gix<-r*exp(-0.5*((abs(xopt-xt)*(xstr)+abs(yopt-yt)*(1-xstr))/w)^2)
      
      dB<-(gix*(State[1]/(State[1]+K))-m)*State[-1]
      dR<-a*(S-State[1])-sum(Q*gix*(State[1]/(State[1]+K)*State[-1]))
      
      return(list(c(dR, dB)))
    })
  }
  
  #########################
  # make climate variables
  #########################
  xlist<- numeric(length(times))
  ylist<- numeric(length(times))
  
  xmid<-pars$x_mean_sd[1]
  ymid<-pars$y_mean_sd[1]
    
  xlist[1]<-runif(1, min=pars$x_mean_sd[1]-2*pars$x_mean_sd[2], max=pars$x_mean_sd[1]+2*pars$x_mean_sd[2])
  ylist[1]<-runif(1, min=pars$y_mean_sd[1]-2*pars$y_mean_sd[2], max=pars$y_mean_sd[2]+2*pars$y_mean_sd[2])
  Wxlist<-rnorm(length(times), 0, Wrates[1]) #White noise for x
  Wylist<-rnorm(length(times), 0, Wrates[2]) #White noise for y
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
  yini<-c(R, B)
  
  
  out   <- ode(yini, times, Bmod, pars)
  
  return(list(out=out, pars=pars))
}


################################################################
# Make plots of results
################################################################
plot_output<-function(ode_result) {
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
  plot(out[,1], pars$xlist, type="s", xlab="time", ylab="x")
  plot(out[,1], pars$ylist, type="s", xlab="time", ylab="y")
  
  #Plot niche space
  plot(range(pars$xopt), c(0,1), xlab="x", ylab="f", type="n", main="x niche space")
  dlist<-seq(20,30,by=0.1)
  for(i in pars$xopt) {
    lines(dlist, pars$r*exp(-0.5*((abs(i-dlist))/pars$w)^2))
  }
  
  plot(range(pars$yopt), c(0,1), xlab="y", ylab="f", type="n", main="y niche space")
  dlist<-seq(0,10,by=0.1)
  for(i in pars$yopt) {
    lines(dlist, pars$r*exp(-0.5*((abs(i-dlist))/pars$w)^2))
  }
}

################################################################
# Make E plot
################################################################
makeEplot_environment<-function(ode_result, target_sp, predstep=10, tau=1, maxE=20) {
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
  
  for(E in 2:(maxE)) {
    Emat[E-1,"B"]<-SSR_pred(A=Biomass_total, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"T"]<-SSR_pred(A=temperature, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"M"]<-SSR_pred(A=moisture, E=E, predstep=predstep, tau=tau)$rho
  }
  matplot(2:(maxE), Emat, type="l", col=1:3, lty=1:3,
          xlab="E", ylab="rho", lwd=2,
          main=paste("pred. step = ", predstep, "; tau = ", tau, sep=""))
  legend("bottomleft", c("Biomass", "Temperature", "Moisture"), lty=1:3, col=1:3, lwd=2, bty="n")
  
  E_temp<-min(which(Emat[,2]>=(quantile(Emat[,2], 0.95)-0.01)))+1
  E_moist<-min(which(Emat[,3]>=(quantile(Emat[,3], 0.95)-0.01)))+1
  E_biom<-min(which(Emat[,1]>=(quantile(Emat[,1], 0.95)-0.01)))+1
  
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
    Emat[E-1,"B"]<-SSR_pred(A=Biomass_total_1, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"T"]<-SSR_pred(A=Biomass_total_2, E=E, predstep=predstep, tau=tau)$rho
  }
  matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,
          xlab="E", ylab="rho", lwd=2,
          main=paste("pred. step = ", predstep, "; tau = ", tau, sep=""))
  legend("bottomleft", c("Biomass 1", "Biomass 2"), lty=1:2, col=1:2, lwd=2, bty="n")
  
  E_biom_2<-min(which(Emat[,2]>=(quantile(Emat[,2], 0.95)-0.01)))+1
  E_biom_1<-min(which(Emat[,1]>=(quantile(Emat[,1], 0.95)-0.01)))+1
  
  return(list(Emat=Emat, Euse=c(E_biom_1=E_biom_1, E_biom_2=E_biom_2)))
}


################################################################
# CCM environment function
################################################################
doCCM_environment<-function(ode_result, target_sp, predstep=10, tau=1, maxE=20) {
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
    
  DesiredL<-unique(round(c(seq(10, median(pars$times), length=20), seq(median(pars$times), max(pars$times)-10, length=10))),0)
  if(max(DesiredL)>length(Biomass_total)) {DesiredL<-1:length(Biomass_total)}
  
  Emat<-matrix(nrow=maxE-1, ncol=3); colnames(Emat)<-c("B", "T", "M")
  
  for(E in 2:maxE) {
    Emat[E-1,"B"]<-SSR_pred(A=Biomass_total, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"T"]<-SSR_pred(A=temperature, E=E, predstep=predstep, tau=tau)$rho
    Emat[E-1,"M"]<-SSR_pred(A=moisture, E=E, predstep=predstep, tau=tau)$rho
  }
  matplot(2:maxE, Emat, type="l", col=1:3, lty=1:3,
          xlab="E", ylab="rho", lwd=2,
          main=paste("pred. step = ", predstep, "; tau = ", tau, sep=""))
  legend("bottomleft", c("Biomass", "Temperature", "Moisture"), lty=1:3, col=1:3, lwd=2, bty="n")
  
  E_temp<-min(which(Emat[,2]>=(quantile(Emat[,2], 0.95)-0.01)))+1
  E_moist<-min(which(Emat[,3]>=(quantile(Emat[,3], 0.95)-0.01)))+1
  E_biom<-min(which(Emat[,1]>=(quantile(Emat[,1], 0.95)-0.01)))+1
  
  T_cause_B<-CCM_varL(A=temperature, B=Biomass_total, E=E_temp, DesiredL=DesiredL)
  M_cause_B<-CCM_varL(A=moisture, B=Biomass_total, E=E_moist, DesiredL=DesiredL)
  
  B_cause_T<-CCM_varL(A=Biomass_total, B=temperature, E=E_biom, DesiredL=DesiredL)
  B_cause_M<-CCM_varL(A=Biomass_total, B=moisture, E=E_biom, DesiredL=DesiredL)
  
  return(list(T_cause_B=T_cause_B, M_cause_B=M_cause_B,
                B_cause_T=B_cause_T, B_cause_M=B_cause_M, Emat=Emat, Euse=c(E_temp=E_temp, E_moist=E_moist, E_biom=E_biom),
              datalist=list(Biomass_total=Biomass_total, temperature=temperature, moisture=moisture)))
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
    pred_mat[predstep,"B"]<-SSR_pred(A=ccm_out$datalist$Biomass_total, E=ccm_out$Euse["E_biom"], predstep=predstep, tau=tau)$rho
    pred_mat[predstep,"T"]<-SSR_pred(A=ccm_out$datalist$temperature, E=ccm_out$Euse["E_temp"], predstep=predstep, tau=tau)$rho
    pred_mat[predstep,"M"]<-SSR_pred(A=ccm_out$datalist$moisture, E=ccm_out$Euse["E_moist"], predstep=predstep, tau=tau)$rho
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
plot_ccm<-function(ccm_out, ylimits=c(0, 1)) {
  #########################
  # Plot CCM
  #########################
  with(ccm_out, {
  minL<-max(length(T_cause_B$Lobs), length(T_cause_B$rho), length(B_cause_T$rho), length(M_cause_B$rho), length(B_cause_M$rho))
  matplot(T_cause_B$Lobs[1:minL], cbind(T_cause_B$rho[1:minL], B_cause_T$rho[1:minL], M_cause_B$rho[1:minL], B_cause_M$rho[1:minL]), type="l", main="CCM results", xlab="Library length", ylab="rho", lwd=2,
          ylim=ylimits, lty=c(1,2,1,2), col=c("black", "black", "red", "red"))
  legend("topleft", c("Temp. causes Biom.", "Biom. causes Temp.", "Moist. causes Biom.", "Biom. causes Moist."), col=c("black", "black", "red", "red"), lty=c(1,2,1,2), lwd=2, bty="n")
  abline(h=0, lty=3, col="darkgrey", lwd=2)
  })
}


################################################################
# Check for existence of manifold SPECIES
################################################################
ssr_sp<-function(ccm_sp, predstepmax=10, tau=1) {
  
  pred_mat<-matrix(nrow=predstepmax, ncol=2)
  colnames(pred_mat)<-c("B", "T")
  for(predstep in 1:predstepmax) {
    pred_mat[predstep,"B"]<-SSR_pred(A=ccm_sp$datalist$Biomass_total_1, E=ccm_sp$Euse["E_biom_1"], predstep=predstep, tau=tau)$rho
    pred_mat[predstep,"T"]<-SSR_pred(A=ccm_sp$datalist$Biomass_total_2, E=ccm_sp$Euse["E_biom_2"], predstep=predstep, tau=tau)$rho
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





