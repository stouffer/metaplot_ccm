#Load functions
setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")

#Set simulation parameters
times=seq(1, 1001, by = 1)
number_of_chains=1

#Set constant parameters
seednum=1114
x_mean_sd = c(0, 1)
Wrates=c(1)
OUrates=c(0.1)
a = 1
w = 1
K = 10
m = 0.1
Q = 2
r = 1
S = 100
nspec<-12
B=c(rep(0.01, nspec))
R=100
seednum=1989
xopt<-seq(22, 28, length=nspec)

#Function for diffeq
Bmod <- function(Time, State, Pars, input) {
  with(as.list(c(State, Pars)), {
    
    State[-c(1)][State[-c(1)]<=0]<-0
    pS<-State[-c(1)]/sum(State[-c(1)])
    try(Shannon_tmp<-(-pS*log(pS)), silent=T)
    Shannon_tmp[!is.finite(Shannon_tmp)]<-0
    Shannon<-exp(sum(Shannon_tmp))
        
    gix<-r*exp(-0.5*((xopt-input(Time))/w)^2)
    dB<-(gix*(State[1]/(State[1]+K))-m)*State[-1]
    
    diffsupply<-Shannon*Sstr
    dR<-(a+diffsupply)*(S-State[1])-sum(Q*gix*(State[1]/(State[1]+K)*State[-1]))
    return(list(c(dR, dB)))
  })
}

#External signal with rectangle impulse
signal <- data.frame(times = times,
                     import = rep(0, length(times)))

signal$import<-runif(nrow(signal), min=20, max=30)

sigimp <- approxfun(signal$times, signal$import, rule = 2)

#Set interaction strength
Sstr=1

#Run diffeq function
yini<-c(R, B)
require(deSolve)
source("CCM_multi_Library_tx.R")
program_init_bootstrap_short()

times<-seq(1, max(times), by=1)
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

out_tmp   <- ode(yini, times, Bmod, pars, input=sigimp)

#Plot output from diffeq
#matplot(out_tmp[,1], out_tmp[,2], type="l", xlab="time", ylab="R")
#matplot(out_tmp[,1], out_tmp[,3:14], type="l", xlab="time", ylab="B")

#Calculate Shannon
out_tmp[,3:14][out_tmp[,3:14]<=0]<-0
B_tot<-rowSums(out_tmp[,3:14])
plot(out_tmp[,1], B_tot, type="l", xlab="time", ylab="B")
pB_tot<-out_tmp[,3:14]/B_tot
lB_tot<-log(pB_tot)
lB_tot[!is.finite(lB_tot)]<-0
sh_div<--(pB_tot)*lB_tot
eH<-exp(rowSums(sh_div))

#plot(out_tmp[,1], eH, type="l")

#Scale variables
eHccm<-scale(eH[-c(1:200)])
B_totccm<-scale(B_tot[-c(1:200)])

#Calculate optimal E
maxE<-10
Emat<-matrix(nrow=maxE-1, ncol=2); colnames(Emat)<-c("B", "eH")
for(E in 2:maxE) {
  Emat[E-1,"B"]<-SSR_pred_boot(A=B_totccm, E=E, predstep=1, tau=1)$rho
  Emat[E-1,"eH"]<-SSR_pred_boot(A=eHccm, E=E, predstep=1, tau=1)$rho
}
matplot(2:maxE, Emat, type="l", col=1:2, lty=1:2,
          xlab="E", ylab="rho", lwd=2)
legend("bottomleft", c("Biomass", "eH"), lty=1:3, col=1:3, lwd=2, bty="n")

E_Bt<-which.max(Emat[,1])+1
E_eH<-which.max(Emat[,2])+1

#Do CCM
CCM_boot_eHccm<-CCM_boot(eHccm, B_totccm, E_eH, tau=1, iterations=100)
CCM_boot_Btccm<-CCM_boot(B_totccm, eHccm, E_Bt, tau=1, iterations=100)

#Save outputs

#Plot outputs
lng<-length(CCM_boot_eHccm$Lobs)
matplot(CCM_boot_eHccm$Lobs[-lng], cbind(
  CCM_boot_eHccm$rho[-lng],
  CCM_boot_eHccm$rho[-lng]-CCM_boot_eHccm$sdevrho[-lng],
  CCM_boot_eHccm$rho[-lng]+CCM_boot_eHccm$sdevrho[-lng]),
        type="l", lty=c(1,2,2), col=1, lwd=c(2,1,1))

lng<-length(CCM_boot_Btccm$Lobs)
matplot(CCM_boot_Btccm$Lobs[-lng], cbind(
  CCM_boot_Btccm$rho[-lng],
  CCM_boot_Btccm$rho[-lng]-CCM_boot_Btccm$sdevrho[-lng],
  CCM_boot_Btccm$rho[-lng]+CCM_boot_Btccm$sdevrho[-lng]),
        type="l", lty=c(1,2,2), col=1, lwd=c(2,1,1))


