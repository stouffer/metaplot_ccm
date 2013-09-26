# load SSR functions
path="~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/"
setwd(path)

file="CCM_varL_130409"
if(is.loaded(file)) {dyn.unload(paste(file,".so",sep=""))}
system(paste("R CMD SHLIB ",file,".c",sep=""))
dyn.load(paste(file,".so",sep=""))

CCM_varL<-function(A, B, E, tau=1, DesiredL=((tau*(E-1)+(E+1)):length(A)-E+2), plengtht=length(A), pLibLength=length(A), Aest=rep(0, length(A)), rho=rep(0,length(A))) {
  if(plengtht>pLibLength) {plengtht=pLibLength}
  DesiredL<-DesiredL+E-2
  out<-.C("CCM_varL_130409", A=as.double(A), Aest=as.double(Aest), B=as.double(B), rho=as.double(rho), E=as.integer(E), tau=as.integer(tau),
          plength=as.integer(plengtht), pLibLength=as.integer(pLibLength),DesiredL=as.integer(DesiredL), plengthDesL=as.integer(length(DesiredL)))
  out$Aest[out$Aest==0]<-NA #Mark values that were not calculated  
  return(list(A=out$A, Aest=out$Aest, B=out$B, rho=out$rho[out$rho!=0], Lobs=(1:length(A))[out$rho!=0]-E+1, E=out$E, tau=out$tau, FULLinfo=out))
}


tau<-1
E<-2
L<-78
plengtht<-78


tlist<-1:1000
y<-sin(tlist/10)+runif(1000)

x<-numeric(1000); x[1]<-0.1
for(i in 2:1000) {
  x[i]<-x[i-1]*(3.8-3.8*x[i-1] - 0.1*y[i-1])
  if(x[i]<0) {x[i]<-0.1}
  if(x[i]>1) {x[i]<-0.9}
}

pdf("ccm_est_1.pdf", width=8, height=6)
plot(tlist, y, type="l", main="A (black) causes B (red)", xlab="Library Length", ylab="Value")
lines(tlist, x, col="red")
dev.off()

Lcalc<-c(seq(10,100,by=5), seq(200, 950, by=50))

my_1<-CCM_varL(A=y, B=x, E=E, DesiredL=Lcalc)
my_2<-CCM_varL(A=x, B=y, E=E, DesiredL=Lcalc)


pdf("ccm_est_2.pdf", width=8, height=6)
plot(Lcalc, my_1$rho[my_1$rho!=0], type="l", xlim=c(0, 1000), ylim=c(-0.2,1), lwd=4,
     xlab="Library Length", ylab="rho", main="CCM estimate of causality")
lines(Lcalc, my_2$rho[my_1$rho!=0], lty=2, lwd=4, col="red")
legend(400, 0.6, c("A causes B (correct)", "B causes A (correct)"), lty=c(1,2), lwd=4, col=c(1,2), cex=1.4, bty="n")
dev.off()



x<-numeric(1000); x[1]<-0.1
y<-numeric(1000); y[1]<-0.3
for(i in 2:1000) {
  x[i]<-x[i-1]*(3.8-3.8*x[i-1] - 0.8*y[i-1])
  y[i]<-y[i-1]*(3.2-3.2*y[i-1] - 0.0*y[i-1])
  if(x[i]<0) {x[i]<-0.1}
  if(x[i]>1) {x[i]<-0.9}
  if(y[i]<0) {y[i]<-0.1}
  if(y[i]>1) {y[i]<-0.9}
}

pdf("ccm_est_3.pdf", width=8, height=6)
plot(tlist, y, type="l", main="A (black) causes B (red)", xlab="Library Length", ylab="Value",
     ylim=c(0,1))
lines(tlist, x, col="red")
dev.off()


my_1_2<-CCM_varL(A=y, B=x, E=E, DesiredL=Lcalc)
my_2_2<-CCM_varL(A=x, B=y, E=E, DesiredL=Lcalc)

pdf("ccm_est_4.pdf", width=8, height=6)
plot(Lcalc, my_1_2$rho[my_1$rho!=0], type="l", xlim=c(0, 1000), ylim=c(-0.2,1), lwd=4,
     xlab="Library Length", ylab="rho", main="CCM estimate of causality")
lines(Lcalc, my_2_2$rho[my_1$rho!=0], lty=2, lwd=4, col="red")
legend(400, 0.6, c("A causes B (correct)", "B causes A (incorrect)"), lty=c(1,2), lwd=4, col=c(1,2), cex=1.4, bty="n")

dev.off()
