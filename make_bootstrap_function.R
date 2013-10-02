setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")

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


#### Make fake data to analyze
x<-ode_result$out[,2]
y<-ode_result$out[,3]

x[seq(10, 200, by=10)]<-NA
y[seq(10, 200, by=10)]<-NA

time<-rep(c(1:9, NA), 20)

#### Make bootstrap lists
