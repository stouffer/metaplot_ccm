diag(amat)<-1


mmat<-rep(0.0001, 16)



####### Run simulation in one grid

Ntmat<-matrix(nrow=20,ncol=16);
Ntmat[1,]<-0.0

j<-1
Kmat[j]<-1

Ntmat<-matrix(nrow=20, ncol=nspec)
Ntmat[1,]<-0.01

for(i in 2:20) {
  Ntmat[i,]<-Ntmat[i-1,]+Ntmat[i-1,]*rmat*(1-Ntmat[i-1,]%*%amat)/Kmat[j]
}

matplot(1:20, Ntmat, type="l")

CCmod <- function(Time, Ntvec, Pars) {
  with(as.list(c(Ntvec, Pars)), {
    
    dNtvec<-Ntmat*rmat*(1-Ntmat%*%amat)/Kmat[j]
    
    return(list(c(dNtvec)))
  })
}

nspec<-6
pars  <- list(rmat=sort(runif(nspec, min=0.01, max=1)),
              Kmat<-rep(1, nspec),
              amat<-matrix(nrow=nspec, ncol=nspec, data=runif(nspec^2)))

yini  <- c(rep(0.01, nspec))
times <- seq(0, 100, by = 1)

CCmod(1, yini, pars)
require(deSolve)
out   <- ode(yini, times, CCmod, pars)



matplot(out[,1], out[,-1], type="l")


prodmat<-rowSums(Ntmat)
divmat<-rowSums(Ntmat>0.001)

plot(1:1000, prodmat, type='l')
plot(1:1000, divmat, type='l')


#######
require(deSolve)

CCmod <- function(Time, Ntvec, Pars) {
  with(as.list(c(Ntvec, Pars)), {
    veclen<-length(Ntvec)
    
    dNtvec<-cvec*Ntvec*(1-D-sum(Ntvec)) - mvec*Ntvec
    dNtvec[-1]<-dNtvec[-1] - cumsum(Ntvec[-veclen]*cvec[-veclen])*Ntvec[-1]
    dNtvec[-veclen]<-dNtvec[-veclen] + rev(cumsum(rev(Ntvec)[-veclen]))*cvec[-veclen]*Ntvec[-veclen]
    
    dNtvec[Ntvec<0]<--Ntvec[Ntvec<0]
    dNtvec[Ntvec>1]<--Ntvec[Ntvec>1]
    
    return(list(c(dNtvec)))
  })
}

nspec<-6
pars  <- list(cvec = c(1,2,4,16,32,64),
           mvec = rep(0.2, nspec),
           D = 0.25)

yini  <- c(rep(0.01, nspec))
times <- seq(0, 20, by = 1)

CCmod(1, yini, pars)

out   <- ode(yini, times, CCmod, pars)

## User specified plotting
matplot(out[ , 1], out[ , -1], type = "l", xlab = "time", ylab = "pcover",
        main = "CCmod", lwd = 2)
legend("topright", as.character(pars$cvec), col = 1:nspec, lty = 1:nspec, cex=0.8, ncol=3)

