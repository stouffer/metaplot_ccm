setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")
source("CCM_multi_Library.R")


system.time({
  #Initialize and make data
  program_init_bootstrap_short()
  ode_result<-make_comp_data_boot_Cfxn(seednum=123, xstr=0,  times=seq(1, 100, by = 1), number_of_chains=5)
  plot_output_boot(ode_result)
  
  #Compare species to environment
  eplot_out_boot<-makeEplot_environment_boot(ode_result, "all", tau=1, predstep=1, maxE=5)
  
  ccm_out<-doCCM_environment(ode_result=ode_result, target_sp="all", predstep=1, tau=1, maxE=5, iterations=100, twoway=FALSE)
  ssr_out<-ssr_data(ccm_out, predstepmax=10, tau=1)
  plot_ccm(ccm_out,  ylimits=c(-0.2, 1), twoway=FALSE)
  
  test_out<-testccm(ccm_out)
  test_out
})




#Simulate:
#dat_out<-matrix(nrow=(iter_use<-100), ncol=11)
colnames(dat_out)<-c("Tdiff", "Tpdiff", "Trho", "Tprho", "Tcause", "Mdiff", "Mpdiff", "Mrho", "Mprho", "Mcause", "xstr")
require(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
iter_max<-100
system.time({
dat_out<-foreach(z_iter_counter=1:iter_max, .combine=rbind) %dopar% {
  program_init_bootstrap_short()
  #Initialize and make data
    xstr<-runif(1, 0, 1)
  
    ode_result<-FALSE
    
  while(length(ode_result)<2) { #catch strange deSolve error
  try({
      ode_result<-make_comp_data_boot_Cfxn(seednum=FALSE, xstr=xstr,  times=seq(1, 100, by = 1), number_of_chains=5)
    })
  
  if(length(ode_result)==1) {
    detach("package:deSolve", unload=TRUE)
    require(deSolve)
  }
  }
  
  #plot_output_boot(ode_result)
  
  #Compare species to environment
  #eplot_out_boot<-makeEplot_environment_boot(ode_result, "all", tau=1, predstep=1, maxE=5)
  
  ccm_out<-doCCM_environment(ode_result=ode_result, target_sp="all", predstep=1, tau=1, maxE=5, iterations=100, twoway=FALSE, make_plot=FALSE)
  #ssr_out<-ssr_data(ccm_out, predstepmax=10, tau=1)
  #plot_ccm(ccm_out,  ylimits=c(-0.2, 1), twoway=FALSE)
  
  test_out<-testccm(ccm_out)
  c(test_out$T_test, test_out$M_test, xstr)
}
})

colnames(dat_out)<-c("Tdiff", "Tpdiff", "Trho", "Tprho", "Tcause", "Mdiff", "Mpdiff", "Mrho", "Mprho", "Mcause", "xstr")
stopImplicitCluster()

#Plot results
matplot(dat_out[,"xstr"], cbind(dat_out[,"Trho"], dat_out[,"Mrho"]), type="p", pch=c(1,2),
        xlab="xstr", ylab="rho")

matplot(dat_out[,"xstr"], cbind(dat_out[,"Tcause"], dat_out[,"Mcause"]), type="p", pch=c(1,2),
        xlab="xstr", ylab="rho")

mod1<-glm(dat_out[,"Tcause"]~dat_out[,"xstr"], family="binomial")
summary(mod1)

mod2<-glm(dat_out[,"Mcause"]~dat_out[,"xstr"], family="binomial")
summary(mod2)
####
The below still need to be fixed
##


system.time({
  #Initialize and make data
  program_init_bootstrap()
  ode_result<-make_comp_data_boot_Cfxn(seednum=123, xstr=0,  times=seq(1, 100, by = 1), number_of_chains=5)
  plot_output_boot(ode_result)
  
  #Compare species to environment
  eplot_out_boot<-makeEplot_environment_boot(ode_result, 5, tau=1, predstep=1, maxE=5)
  
  ccm_out<-doCCM_environment(ode_result=ode_result, target_sp=1, predstep=1, tau=1, maxE=5, iterations=100, twoway=FALSE)
  ssr_out<-ssr_data(ccm_out, predstepmax=10, tau=1)
  plot_ccm(ccm_out,  ylimits=c(-0.2, 1), twoway=FALSE)
})


#Try different levels of averaging
ode_result<-make_comp_data(seednum=1000, xstr=0.5,  times=seq(1, 2000, by = 1), x_mean_sd = c(25, 1), y_mean_sd = c(5, 1))
ode_agg<-odeout_average(ode_result, env_avg=10)

eplot_out<-makeEplot_environment(ode_agg, "all", tau=1, predstep=1, maxE=5)
ccm_out<-doCCM_environment(ode_agg, "all", tau=1, predstep=1, maxE=5)
ssr_out<-ssr_data(ccm_out, predstepmax=5, tau=1)
plot_ccm(ccm_out,  ylimits=c(-0.2, 1))