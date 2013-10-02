setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")
source("CCM_single_Library.R")

system.time({
program_init()
ode_result<-make_comp_data(seednum=1000, xstr=0.5,  times=seq(1, 200, by = 1), x_mean_sd = c(25, 1), y_mean_sd = c(5, 1))
plot_output(ode_result)

#Compare species to environment
eplot_out<-makeEplot_environment(ode_result, "all", tau=2, predstep=tau)
ccm_out<-doCCM_environment(ode_result, "all", tau=2)
ssr_out<-ssr_data(ccm_out, predstepmax=10, tau=5)
plot_ccm(ccm_out,  ylimits=c(-0.2, 1))

#Compare species
eplot_sp_out<-makeEplot_species(ode_result, 1, 2, tau=2, predstep=tau)
ccm_sp<-doCCM_species(ode_result, 1, 2, predstep=5, tau=1)
ssr_sp_out<-ssr_sp(ccm_sp, predstepmax=10, tau=5)
plot_sp(ccm_sp,  ylimits=c(-0.2, 1))
})


#Try different levels of averaging
ode_result<-make_comp_data(seednum=1000, xstr=0.5,  times=seq(1, 2000, by = 1), x_mean_sd = c(25, 1), y_mean_sd = c(5, 1))
ode_agg<-odeout_average(ode_result, env_avg=10)

eplot_out<-makeEplot_environment(ode_agg, "all", tau=1, predstep=1, maxE=5)
ccm_out<-doCCM_environment(ode_agg, "all", tau=1, predstep=1, maxE=5)
ssr_out<-ssr_data(ccm_out, predstepmax=5, tau=1)
plot_ccm(ccm_out,  ylimits=c(-0.2, 1))

