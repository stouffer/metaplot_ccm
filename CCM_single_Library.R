setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")

################################################################
# Initialize the C functions for CCM
################################################################
program_init()

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
}

CCM_varL<-function(A, B, E, tau=1, DesiredL=((tau*(E-1)+(E+1)):length(A)-E+2), plengtht=length(A), pLibLength=length(A), Aest=rep(0, length(A)), rho=rep(0,length(A))) {
  if(plengtht>pLibLength) {plengtht=pLibLength}
  DesiredL<-DesiredL+E-2
  out<-.C("CCM_varL_130409", A=as.double(A), Aest=as.double(Aest), B=as.double(B), rho=as.double(rho), E=as.integer(E), tau=as.integer(tau),
          plength=as.integer(plengtht), pLibLength=as.integer(pLibLength),DesiredL=as.integer(DesiredL), plengthDesL=as.integer(length(DesiredL)))
  out$Aest[out$Aest==0]<-NA #Mark values that were not calculated  
  return(list(A=out$A, Aest=out$Aest, B=out$B, rho=out$rho[out$rho!=0], Lobs=(1:length(A))[out$rho!=0]-E+1, E=out$E, tau=out$tau, FULLinfo=out))
}

################################################################
# Generate time series function
################################################################
ode_result<-make_comp_data(seednum=2000)

make_comp_data<-function(nspec=12,
                         pars=list(a = 1,
                                       w = 1,
                                       K = 10,
                                       m = 0.1,
                                       Q = 2,
                                       r = 1,
                                       S = 100,
                                       xrange = c(20, 30),
                                       yrange = c(0, 10),
                                       xstr=0),
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
  
  xmid<-mean(pars$xrange)
  ymid<-mean(pars$yrange)
    
  xlist[1]<-runif(1, min=pars$xrange[1], max=pars$xrange[2])
  ylist[1]<-runif(1, min=pars$yrange[1], max=pars$yrange[2])
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
  
  pars$xlist<-xlist
  pars$ylist<-ylist
  
  pars$xopt<-runif(nspec, min=pars$xrange[1], max=pars$xrange[2])
  pars$yopt<-runif(nspec, min=pars$yrange[1], max=pars$yrange[2])
  yini<-c(R, B)
  
  
  out   <- ode(yini, times, Bmod, pars)
  
  return(list(out=out, pars=pars))
}


################################################################
# Make plots of results
################################################################
plot_output(ode_result)

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
  plot(pars$xrange, c(0,1), xlab="x", ylab="f", type="n", main="x niche space")
  dlist<-seq(20,30,by=0.1)
  for(i in pars$xopt) {
    lines(dlist, pars$r*exp(-0.5*((abs(i-dlist))/pars$w)^2))
  }
  
  plot(pars$yrange, c(0,1), xlab="y", ylab="f", type="n", main="y niche space")
  dlist<-seq(0,10,by=0.1)
  for(i in pars$yopt) {
    lines(dlist, pars$r*exp(-0.5*((abs(i-dlist))/pars$w)^2))
  }
}

################################################################
# CCM environment function
################################################################
ccm_out<-doCCM_environment(ode_result, "all")

doCCM_environment<-function(ode_result, target_sp) {
  out<-ode_result[["out"]]
  pars<-ode_result[["pars"]]
  #########################
  # Run CCM
  #########################
    
  if(target_sp=="all") {
    Biomass_total<-rowSums(out[ , -c(1,2)])
  } else {
    Biomass_total<-out[ , target_sp]
  }
  temperature<-pars$xlist
  moisture<-pars$ylist
    
  DesiredL<-unique(round(c(seq(10, median(pars$times), length=20), seq(median(pars$times), max(pars$times)-10, length=10))),0)
  
  source("SSR_functions.R")
  Emat<-matrix(nrow=20, ncol=3); colnames(Emat)<-c("B", "T", "M")
  
  for(E in 2:21) {
    Emat[E-1,"B"]<-univariate_SSR(Biomass_total, c(1,200), c(1,200), E, tau = 1, tp = 2, b = E+1)$stats$rho
    Emat[E-1,"T"]<-univariate_SSR(temperature, c(1,200), c(1,200), E, tau = 1, tp = 1, b = E+1)$stats$rho
    Emat[E-1,"M"]<-univariate_SSR(moisture, c(1,200), c(1,200), E, tau = 1, tp = 1, b = E+1)$stats$rho
  }
  matplot(1:20, Emat, type="l", col=1:3, lty=1:3)
  
  
  (E<-min(which(Emat[,2]>=(quantile(Emat[,2], 0.95)-0.01)))+1)
  E<-3
  T_cause_B<-CCM_varL(A=temperature, B=Biomass_total, E=E, DesiredL=DesiredL)
  (E<-min(which(Emat[,3]>=(quantile(Emat[,3], 0.95)-0.01)))+1)
  M_cause_B<-CCM_varL(A=moisture, B=Biomass_total, E=E, DesiredL=DesiredL)
  
  (E<-min(which(Emat[,1]>=(quantile(Emat[,1], 0.95)-0.01)))+1)
  B_cause_T<-CCM_varL(A=Biomass_total, B=temperature, E=E, DesiredL=DesiredL)
  B_cause_M<-CCM_varL(A=Biomass_total, B=moisture, E=E, DesiredL=DesiredL)
  
  #########################
  # Plot CCM
  #########################
  plot(T_cause_B$Lobs, T_cause_B$rho, type="l", main="T_cause_B", ylab="Library length", xlab="rho", lwd=2)
  plot(M_cause_B$Lobs, M_cause_B$rho, type="l", main="M_cause_B", ylab="Library length", xlab="rho", lwd=2)
  
  plot(B_cause_T$Lobs, B_cause_T$rho, type="l", main="B_cause_T", ylab="Library length", xlab="rho", lwd=2)
  plot(B_cause_M$Lobs, B_cause_M$rho, type="l", main="B_cause_M", ylab="Library length", xlab="rho", lwd=2)
    
  return(list(T_cause_B=T_cause_B, M_cause_B=M_cause_B,
                B_cause_T=B_cause_T, B_cause_M=B_cause_M))
}



