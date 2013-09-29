setwd("~/Dropbox/Active Work/Sugihara/Metacommunity_CCM/")
source("CCM_single_Library.R")

system.time({
program_init()
ode_result<-make_comp_data(seednum=2000, xstr=1,  times=seq(1, 1000, by = 1))
plot_output(ode_result)
ccm_out<-doCCM_environment(ode_result, "all")
ssr_out<-ssr_data(ccm_out, predstepmax=10)
plot_ccm(ccm_out,  ylimits=c(-0.2, 0.5))
})
