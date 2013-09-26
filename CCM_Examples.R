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
plot(tlist, y, type="l")
x<-numeric(1000); x[1]<-0.1
for(i in 2:1000) {
  x[i]<-x[i-1]*(3.8-3.8*x[i-1] - 0.1*y[i-1])
  if(x[i]<0) {x[i]<-0.1}
  if(x[i]>1) {x[i]<-0.9}
}
lines(tlist, x, col="red")

my_1<-CCM_varL(A=y, B=x, E=E)
my_2<-CCM_varL(A=x, B=y, E=E)
