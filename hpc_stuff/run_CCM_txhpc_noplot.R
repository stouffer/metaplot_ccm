#Load functions
#setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")

#Set simulation parameters
fname<-"CCM_noplot"

times=seq(1, 1201, by = 100)
number_of_chains=5
clnumber<-2 #Number of CPU nodes to use
unique_signal=TRUE #Should each plot get its own signal?
rare_signal=FALSE #Are plots viewed less than 1 per timestep?
rare_by=100 #How often are they viewed?

#Set constant parameters
seednum=1114
x_mean_sd = c(0, 1)
Wrates=c(1)
OUrates=c(0.1)
a = 1; w = 1; K = 10; m = 0.1; Q = 2; r = 1; S = 100; R=100; nspec<-12
B=c(rep(0.01, nspec))
xopt<-seq(22, 28, length=nspec)
datout<-NULL
#Set interaction strength
Sstrlist<-c(0, 0.001, 0.01, 0.1, 1, 2, 4, 8, 16, 32, 64, 128)
minreading<-200

#Function for diffeq
Bmod <- function(Time, State, Pars, input) {
  with(as.list(c(State, Pars)), {
    
    #Calculate Shannon diversity for current assemblage
    State[-c(1)][State[-c(1)]<=0]<-0 #State[1] is resource abundance
    pS<-State[-c(1)]/sum(State[-c(1)])
    try(Shannon_tmp<-(-pS*log(pS)), silent=T)
    Shannon_tmp[!is.finite(Shannon_tmp)]<-0
    Shannon<-exp(sum(Shannon_tmp))
        
    gix<-r*exp(-0.5*((xopt-input(Time))/w)^2)
    dB<-(gix*(State[1]/(State[1]+K))-m)*State[-1]
    
    diffsupply<-Shannon*Sstr #Change supply rate proportionally to Shannon diversity
    dR<-(a+diffsupply)*(S-State[1])-sum(Q*gix*(State[1]/(State[1]+K)*State[-1]))
    return(list(c(dR, dB)))
  })
}

require(deSolve)
require(doParallel)
cl <- makeCluster(clnumber)
registerDoParallel(cl)

datout<-foreach(Sstr=Sstrlist, .combine=rbind) %dopar% { #Run simulation for each strength size
  #Run diffeq function
  yini<-c(R, B)
  require(deSolve)
  source("CCM_multi_Library_tx.R")
  program_init_bootstrap_short()
  
  pars=list(a = a,
            w = w,
            K = K,
            m = m,
            Q = Q,
            r = r,
            xopt=xopt,
            x_mean_sd = x_mean_sd,
            times=times,
            Sstr=Sstr)
  
  eHccm<-NULL
  B_totccm<-NULL
  time_ccm<-NULL
  
  if(unique_signal==FALSE) {
    signal <- data.frame(times = times,
                         import = rep(0, length(times)))
    signal$import<-runif(nrow(signal), min=20, max=30)
    sigimp <- approxfun(signal$times, signal$import, rule = 2)
  }
  
  
  for(i in 1:number_of_chains) { #Build plots
    if(unique_signal==TRUE) {
      signal <- data.frame(times = times,
                           import = rep(0, length(times)))
      signal$import<-runif(nrow(signal), min=20, max=30)
      sigimp <- approxfun(signal$times, signal$import, rule = 2)
    }
    
    out_tmp   <- ode(yini, times, Bmod, pars, input=sigimp)
    
    #Plot output from diffeq
    #matplot(out_tmp[,1], out_tmp[,2], type="l", xlab="time", ylab="R")
    #matplot(out_tmp[,1], out_tmp[,3:14], type="l", xlab="time", ylab="B")
    
    #Calculate Shannon
    out_tmp[,3:14][out_tmp[,3:14]<=0|!is.finite(out_tmp[,3:14])]<-0
    B_tot<-rowSums(out_tmp[,3:14], na.rm=T)
    #plot(out_tmp[,1], B_tot, type="l", xlab="time", ylab="B")
    pB_tot<-out_tmp[,3:14]/B_tot
    lB_tot<-log(pB_tot)
    lB_tot[!is.finite(lB_tot)]<-0
    sh_div<--(pB_tot)*lB_tot
    eH<-exp(rowSums(sh_div))
    
    #plot(out_tmp[,1], eH, type="l")
    
    #Save variables
    eHccm<-c(eHccm, NA, eH)
    B_totccm<-c(B_totccm, NA, B_tot)
    time_ccm<-c(time_ccm, NA, out_tmp[,1])
    
  } # End plot building
  
  eHccm<-eHccm[time_ccm>minreading] #Remove first years of establishment
  B_totccm<-B_totccm[time_ccm>minreading] #Remove first years of establishment
  time_ccm<-time_ccm[time_ccm>minreading]
  
  if(rare_signal) { #Is signal only observed occasionally?
    rareseq<-seq(min(time_ccm, na.rm=T), max(time_ccm, na.rm=T), by=rare_by)
    
    eHccm<-eHccm[(time_ccm%in%rareseq)|(is.na(time_ccm))]
    B_totccm<-B_totccm[(time_ccm%in%rareseq)|(is.na(time_ccm))]
    time_ccm<-time_ccm[(time_ccm%in%rareseq)|(is.na(time_ccm))]
  }
  
  #Scale variables
  eHccm<-scale(eHccm)
  B_totccm<-scale(B_totccm)
  
  #Calculate optimal E
  maxE<-10
  Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("B", "eH")
  for(E in 2:maxE) {
    Emat[E-1,"B"]<-SSR_pred_boot(A=B_totccm, E=E, predstep=1, tau=1)$rho
    Emat[E-1,"eH"]<-SSR_pred_boot(A=eHccm, E=E, predstep=1, tau=1)$rho
  }
  #matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,
  #          xlab="E", ylab="rho", lwd=2)
  #legend("bottomleft", c("Biomass", "eH"), lty=1:3, col=1:3, lwd=2, bty="n")
  
  E_Bt<-which.max(Emat[,1])+1
  E_eH<-which.max(Emat[,2])+1
  
  #Do CCM
  CCM_boot_eHccm<-CCM_boot(eHccm, B_totccm, E_eH, tau=1, iterations=100)
  CCM_boot_Btccm<-CCM_boot(B_totccm, eHccm, E_Bt, tau=1, iterations=100)
  
  #Save outputs
  max_lmax<-1:(max(length(CCM_boot_eHccm$Lobs), length(CCM_boot_Btccm$Lobs)))
  cbind(Sstr=Sstr, eHL=CCM_boot_eHccm$Lobs[max_lmax], eHrm=CCM_boot_eHccm$rho[max_lmax], eHrsd=CCM_boot_eHccm$sdevrho[max_lmax],
        BL=CCM_boot_Btccm$Lobs[max_lmax], Brm=CCM_boot_Btccm$rho[max_lmax], Brsd=CCM_boot_Btccm$sdevrho[max_lmax])
  
} #End run simulations

#Write outputs
write.csv(datout, paste(fname, ".csv", sep=""), row.names=FALSE)

#Plot outputs
datout<-data.frame(datout)
plot_limits_rho<-c(min=min(c(datout$eHrm-datout$ehrsd, datout$Brm-datout$Brsd), na.rm=T),
          max=max(c(datout$eHrm+datout$ehrsd, datout$Brm+datout$Brsd), na.rm=T))
plot_limits_L<-c(min=min(c(datout$eHL, datout$BL), na.rm=T),
                   max=max(c(datout$eHL, datout$BL), na.rm=T))

pdf(paste(fname, ".pdf", sep=""), width=8.5, height=11)
par(mfrow=c(1,2))
collst<-rev(grey.colors(length(Sstrlist)))
collst_red<-rev(heat.colors(length(Sstrlist)))

plot(0, 0,type="n", xlab="L", ylab="rho",
        ylim=plot_limits_rho, xlim=plot_limits_L, main="Diversity 'causes' Productivity")
for(i in 1:length(Sstrlist)) {
  L<-datout[datout$Sstr==Sstrlist[i],]$eHL
  ym<-datout[datout$Sstr==Sstrlist[i],]$eHrm
  ysdL<-ym-datout[datout$Sstr==Sstrlist[i],]$eHrsd
  ysdU<-ym+datout[datout$Sstr==Sstrlist[i],]$eHrsd
  
  matlines(L, cbind(ym, ysdL, ysdU),
          type="l", lty=c(1,2,2), col=collst[i], lwd=c(2,1,1))
}

plot(0, 0,type="n", xlab="L", ylab="rho",
     ylim=plot_limits_rho, xlim=plot_limits_L, main="Productivity 'causes' Diversity")
for(i in 1:length(Sstrlist)) {
  L<-datout[datout$Sstr==Sstrlist[i],]$BL
  ym<-datout[datout$Sstr==Sstrlist[i],]$Brm
  ysdL<-ym-datout[datout$Sstr==Sstrlist[i],]$Brsd
  ysdU<-ym+datout[datout$Sstr==Sstrlist[i],]$Brsd
  
  matlines(L, cbind(ym, ysdL, ysdU),
           type="l", lty=c(1,2,2), col=collst_red[i], lwd=c(2,1,1))
}

dev.off()


