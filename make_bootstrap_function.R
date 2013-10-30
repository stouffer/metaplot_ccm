setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")
source("CCM_single_Library.R")

  program_init()
  ode_result<-make_comp_data(seednum=1000, xstr=0.5,  times=seq(1, 200, by = 1), x_mean_sd = c(25, 1), y_mean_sd = c(5, 1))

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
}

CCM_boot<-function(A, B, E, tau=1, DesiredL=((tau*(E-1)+(E+1)):length(A)-E+2), iterations=1000) {
  ## Initialize parameters
  plengtht=length(A); pLibLength=length(A); Aest=rep(0, length(A)); rho=Aest; varrho=Aest
  
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
  lengthacceptablelib<-length(acceptablelib)
  
  #Remove desired libraries that are too long
  rejectedL<-DesiredL[DesiredL>lengthacceptablelib]
  DesiredL<-DesiredL[DesiredL<=lengthacceptablelib]
  
  # Update input to match actual indices
  DesiredL<-DesiredL+E-2
  A[is.na(A)]<-0; B[is.na(B)]<-0 #Make vectors "C safe"

  out<-.C("CCM_bootstrap", A=as.double(A), Aest=as.double(Aest), B=as.double(B), rho=as.double(rho), E=as.integer(E), tau=as.integer(tau),
          plength=as.integer(plengtht), pLibLength=as.integer(pLibLength),DesiredL=as.integer(DesiredL), plengthDesL=as.integer(length(DesiredL)),
          piterations=as.integer(iterations), varrho=as.double(varrho), acceptablelib=as.integer(acceptablelib), plengthacceptablelib=as.integer(lengthacceptablelib))
  out$Aest[out$Aest==0]<-NA #Mark values that were not calculated  
  sdevrho<-sqrt(out$varrho[out$rho!=0]) #Calculate standard deviation
  return(list(A=out$A, Aest=out$Aest, B=out$B, rho=out$rho[out$rho!=0], varrho=out$varrho[out$rho!=0], sdevrho=sdevrho, Lobs=(1:length(A))[out$rho!=0]-E+1, E=out$E, tau=out$tau, FULLinfo=out, rejectedL=rejectedL))
}


#### Make fake data to analyze
program_init_bootstrap()

x<-ode_result$out[,2]
y<-ode_result$out[,3]
system.time(ccm_old_out<-CCM_varL(x[1:150], y[1:150], E=2, tau=1))
  
  
x[seq(10, 200, by=10)]<-NA
y[seq(10, 200, by=10)]<-NA


# Run functions!!
system.time(ccm_out<-CCM_boot(x, y, E=2, tau=1, iterations=1000)) #30 seconds for L=157 with 1000 boots


#Try with NON correct functions
x<-ode_result$out[,2]
y<-runif(length(x))
system.time(ccm_old_out_2<-CCM_varL(x[1:150], y[1:150], E=2, tau=1))

x[seq(10, 200, by=10)]<-NA
y[seq(10, 200, by=10)]<-NA


# Run functions!!
system.time(ccm_out_2<-CCM_boot(x, y, E=2, tau=1, iterations=1000)) #30 seconds for L=157 with 1000 boots


#pdf("bootstrap_v_wrapping_equilized.pdf", width=8.5, height=11)
par(mfrow=c(2,1))
plot(ccm_old_out$Lobs, ccm_old_out$rho, type="l", lwd=2, xlab="Library Length", ylab="rho", main="Yes causation", xlim=c(0, 170), ylim=c(0, 0.75))
legend(100, 0.2, c("Wrapping", "Boot Strapping"), col=c(1,2), lwd=2, bty="n", cex=1.2);
#ccm_out$varrho<-sqrt(ccm_out$varrho)*(1/sqrt(1000))
matlines(ccm_out$Lobs, cbind(ccm_out$rho, ccm_out$rho-ccm_out$sdevrho, ccm_out$rho+ccm_out$sdevrho),
         col=2, lty=c(1,2,2), lwd=2)

plot(ccm_old_out_2$Lobs, ccm_old_out_2$rho, type="l", lwd=2, xlab="Library Length", ylab="rho", main="No causation", ylim=c(0, 0.25), xlim=c(0, 170))
#ccm_out_2$varrho<-sqrt(ccm_out_2$varrho)
matlines(ccm_out_2$Lobs[-196], cbind(ccm_out_2$rho, ccm_out_2$rho-ccm_out_2$sdevrho, ccm_out_2$rho+ccm_out_2$sdevrho)[-196,],
         col=2, lty=c(1,2,2), lwd=2)
#dev.off()