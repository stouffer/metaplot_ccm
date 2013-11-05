#!/usr/bin/env Rscript
source("CCM_multi_Library.R")

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
    
    #while(length(ode_result)<2) { #catch strange deSolve error
     #try({
        ode_result<-make_comp_data_boot_Cfxn(seednum=FALSE, xstr=xstr,  times=seq(1, 100, by = 1), number_of_chains=5)
      #})
      
      #if(length(ode_result)==1) {
        #detach("package:deSolve", unload=TRUE)
        #require(deSolve)
      #}
    #}
    
    #plot_output_boot(ode_result)
    
    #Compare species to environment
    #eplot_out_boot<-makeEplot_environment_boot(ode_result, "all", tau=1, predstep=1, maxE=5)
    
    ccm_out<-doCCM_environment(ode_result=ode_result, target_sp="all", predstep=1, tau=1, maxE=5, iterations=100, twoway=FALSE, make_plot=FALSE)
    #ssr_out<-ssr_data(ccm_out, predstepmax=10, tau=1)
    #plot_ccm(ccm_out,  ylimits=c(-0.2, 1), twoway=FALSE)
    
    test_out<-testccm(ccm_out)
    c(test_out$T_test, test_out$M_test, xstr)
  }
  colnames(dat_out)<-c("Tdiff", "Tpdiff", "Trho", "Tprho", "Tcause", "Mdiff", "Mpdiff", "Mrho", "Mprho", "Mcause", "xstr")
  write.csv(dat_out, "datout.csv", row.names=FALSE)
})