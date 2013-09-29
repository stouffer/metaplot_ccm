####### Lehman biomass model
a = 1 # resource supply rate
w = 1 # width of niche
K = 10 # half-saturation rate
m = 0.1 # mortality/respiration rate
Q = 2 # biomass conversion coefficient
r = 1 # maximum gross productivity
S = 100 # resource supply point
xrange = c(20, 30) # range of x environmental variable
yrange = c(0, 10) # range of y environmental variable

require(deSolve) # load differentail equation solver package
set.seed(1989)

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

nspec<-12
pars  <- list(a = 1,
              w = 1,
              K = 10,
              m = 0.1,
              Q = 2,
              r = 1,
              S = 100,
              xrange = c(20, 30),
              yrange = c(0, 10))


B <- c(rep(0.01, nspec))
R <- 100
yini<-c(R, B)
times <- seq(1, 200, by = 1)

pars$xlist<-runif(max(times), min=pars$xrange[1], max=pars$xrange[2])
pars$ylist<-runif(max(times), min=pars$yrange[1], max=pars$yrange[2])

pars$xstr<-0.5

pars$xopt<-runif(nspec, min=pars$xrange[1], max=pars$xrange[2])
pars$yopt<-runif(nspec, min=pars$yrange[1], max=pars$yrange[2])

out   <- ode(yini, times, Bmod, pars)
########################################
## Plotting
########################################

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
  lines(dlist, r*exp(-0.5*((abs(i-dlist))/w)^2))
}

plot(pars$yrange, c(0,1), xlab="y", ylab="f", type="n", main="y niche space")
dlist<-seq(0,10,by=0.1)
for(i in pars$yopt) {
  lines(dlist, r*exp(-0.5*((abs(i-dlist))/w)^2))
}




