setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/results/")
dat<-read.csv("datout2.csv")

head(dat)

mod_T<-glm(Tcause~xstr, data=dat)
mod_M<-glm(Mcause~xstr, data=dat)

summary(mod_T)
summary(mod_M)


summary(mod_rho_T<-lm(Trho~xstr, data=dat))
summary(mod_rho_M<-lm(Mrho~xstr, data=dat))

modtmp<-lm(Trho~log(xstr), data=dat)

xstr<-seq(0, 1, by=0.001)
Tfit<-predict(mod_rho_T, newdata=data.frame(xstr=xstr), interval="confidence")
Mfit<-predict(mod_rho_M, newdata=data.frame(xstr=xstr), interval="confidence")

#Rho plots
matplot(xstr, cbind(Tfit, Mfit), lty=c(1,2,2), lwd=c(2,1,1), col=c(1,1,1,2,2,2), xlab="xstr", ylab="rho", type="l")

##Prob plots
plot(c(0, 1), c(-0.3,1), type="n", xlab="xstr", ylab="Pr. causation predicted", main="Env T or M")
subs<-sample(1:nrow(dat), 500, rep=F)
with(dat, rug(xstr[Tcause==0][subs]))
with(dat, rug(xstr[Tcause==1][subs], side=3))
Tcfit<-predict(mod_T, newdat=data.frame(xstr=xstr), se=TRUE)
matlines(xstr, cbind(Tcfit$fit, Tcfit$fit-Tcfit$se.fit, Tcfit$fit+Tcfit$se.fit), lwd=c(2,1,1), lty=c(1,2,2), col=1)

sd_levels<-cut(dat$xstr, 5)
sdseq<-c(0.1,0.3,0.5,0.7,0.9)
plotdat<-t(matrix(nrow=3, data=unlist(tapply(dat$Tcause, sd_levels, function(x) cbind(mean(x), mean(x)-sd(x), mean(x)+sd(x))))))

segments(sdseq-0.01, plotdat[,2], sdseq-0.01, plotdat[,3])


subs<-sample(1:nrow(dat), 500, rep=F)
with(dat, rug(xstr[Mcause==0][subs], col=2))
with(dat, rug(xstr[Mcause==1][subs], side=3, col=2))
Mcfit<-predict(mod_M, newdat=data.frame(xstr=xstr), se=TRUE)
matlines(xstr, cbind(Mcfit$fit, Mcfit$fit-Mcfit$se.fit, Mcfit$fit+Mcfit$se.fit), lwd=c(2,1,1), lty=c(1,2,2), col=2)

sd_levels<-cut(dat$xstr, 5)
sdseq<-c(0.1,0.3,0.5,0.7,0.9)
plotdat<-t(matrix(nrow=3, data=unlist(tapply(dat$Mcause, sd_levels, function(x) cbind(mean(x), mean(x)-sd(x), mean(x)+sd(x))))))

segments(sdseq+0.01, plotdat[,2], sdseq+0.01, plotdat[,3], col=2)
par(xpd=TRUE)
legend(0.1, 1.6, c("T", "M"), lty=1, col=c(1,2), lwd=2)

