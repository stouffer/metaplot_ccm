################################################################
# Initialize the C functions for bootstrap CCM
################################################################
#R -d valgrind --vanilla < testccm.R > log.txt 2>&1

setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")

  file2="SSR_predict_boot"
  if(is.loaded(file2)) {dyn.unload(paste(file2,".so",sep=""))}
  dyn.load(paste(file2,".so",sep=""))

################################################################
# Run bootstrap CCM
################################################################

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
  acceptablelib[(length(A)-predstep*tau):length(A)]<-0 #make sure we have enough future spaces for prediction
  acceptablelib<-which(acceptablelib>0)-1 #Convert into positions in C array
  lengthacceptablelib<-length(acceptablelib)
  
  if(E+1>lengthacceptablelib) {
    print(paste("Error - too few records to test E =", E))
  } else {
  #Remove desired libraries that are too long
  A[is.na(A)]<-0; B[is.na(B)]<-0 #Make vectors "C safe"  
  
  out<-.C("SSR_predict_boot", A=as.double(A), Aest=as.double(rep(0, length(A))), B=as.double(B), E=as.integer(E),
          tau=as.integer(tau),pBlength=as.integer(length(B)), pAlength=as.integer(length(A)), predstep=as.integer(predstep),
          prepvec=as.integer(repvec), pmatchSugi=as.integer(matchSugi), acceptablelib=as.integer(acceptablelib), plengthacceptablelib=as.integer(lengthacceptablelib))
  out$rho<-cor(out$A, out$Aest)
  out$Aest[out$Aest==0]<-NA
  return(out)
  }
}



dat<-read.csv("dat_sim.csv")[,-1]
A=dat; B=dat; E=20; tau=1; predstep=10; matchSugi=1
SSR_pred_boot(A=dat, B=dat, E=20, tau=1, predstep=10, matchSugi=1)