source("CCM_multi_Library.R")

require(doParallel)
require(deSolve)

cl <- makeCluster(4)
registerDoParallel(cl)
iter_max<-100
system.time({
  dat_out<-foreach(z_iter_counter=1:iter_max, .combine=rbind) %dopar% {
    program_init_bootstrap_short()
    #Initialize and make data
    xstr<-runif(1, 0, 1)
    
    ode_result<-FALSE
    
    ode_result<-make_comp_data_boot_Cfxn(seednum=FALSE, xstr=xstr,  times=seq(1, 100, by = 1), number_of_chains=5)
    
    ccm_out<-doCCM_environment(ode_result=ode_result, target_sp="all", predstep=1, tau=1, maxE=5, iterations=100, twoway=FALSE, make_plot=FALSE)
    
    test_out<-testccm(ccm_out)
    c(test_out$T_test, test_out$M_test, xstr)
  }
  colnames(dat_out)<-c("Tdiff", "Tpdiff", "Trho", "Tprho", "Tcause", "Mdiff", "Mpdiff", "Mrho", "Mprho", "Mcause", "xstr")
  write.csv(dat_out, "datout.csv", row.names=FALSE)
})