setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")

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
  return(list(A=out$A, Aest=out$Aest, B=out$B, rho=out$rho[out$rho!=0], varrho=out$varrho[out$rho!=0], Lobs=(1:length(A))[out$rho!=0]-E+1, E=out$E, tau=out$tau, FULLinfo=out, rejectedL=rejectedL))
}


#### Make fake data to analyze
x<-ode_result$out[,2]
y<-ode_result$out[,3]

x[seq(10, 200, by=10)]<-NA
y[seq(10, 200, by=10)]<-NA


# Run functions!!
program_init_bootstrap()
system.time(ccm_out<-CCM_boot(x, y, E=2, tau=1, iterations=1000)) #30 seconds for L=157 with 1000 boots

matplot(ccm_out$Lobs, cbind(ccm_out$rho, ccm_out$rho-ccm_out$varrho, ccm_out$rho+ccm_out$varrho), type="l",
        col=1, lty=c(1,2,2), lwd=2,
        xlab="Library Length", ylab="rho")

